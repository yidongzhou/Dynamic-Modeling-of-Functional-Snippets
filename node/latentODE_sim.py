# Import necessary libraries
import os
import argparse
import time
import numpy as np
import matplotlib
matplotlib.use('agg') # Use 'agg' backend for saving plots without displaying
import torch
import torch.nn as nn
import torch.optim as optim
import concurrent.futures

# Set up argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('--adjoint', type=eval, default=False)  # Use adjoint method for memory efficiency
parser.add_argument('--niters', type=int, default=2000)  # Number of training iterations
parser.add_argument('--lr', type=float, default=0.01)  # Learning rate for optimizer
parser.add_argument('--gpu', type=int, default=0)  # GPU index (if available)
parser.add_argument('--train_dir', type=str, default=None)  # Directory for saving checkpoints
parser.add_argument('--type', type=str, default='OU')  # Type of process to simulate ('OU', 'HL', or 'GBM')
parser.add_argument('--x0', type=float, default=0)  # Initial value of the process
parser.add_argument('--m', type=int, default=2)  # Number of sampled data points per trajectory
args = parser.parse_args()


# Import the appropriate ODE solver
if args.adjoint:
    from torchdiffeq import odeint_adjoint as odeint  # Adjoint method for memory efficiency
else:
    from torchdiffeq import odeint


# Functions for different stochastic processes
def ou_step(theta, sigma, delta, x0, z1):
    """Generate the next step in an Ornstein-Uhlenbeck (OU) process."""
    mean = x0 * np.exp(-theta * delta)
    std_dev = np.sqrt(sigma**2 * (1 - np.exp(-2 * theta * delta)) / (2 * theta))
    return mean + std_dev * z1


def hl_step(sigma, t0, t1, x0, z1):
    """Generate the next step in a Ho-Lee (HL) process."""
    mean = x0 + np.sin(t1) - np.sin(t0)
    std_dev = sigma * np.sqrt(t1 - t0)
    return mean + std_dev * z1


def gbm_step(theta, sigma, delta, x0, z1):
    """Generate the next step in a Geometric Brownian Motion (GBM) process."""
    return x0 * np.exp(delta * (theta - sigma**2 / 2) + sigma * np.sqrt(delta) * z1)


def generate_data(n=50, m=2, delta=0.05, nu=0.0, x0=0.0, ptype='OU'):
    """
    Generate n true trajectories and sample m noisy points from each trajectory.

    Parameters:
    - n: Number of trajectories
    - m: Number of sampled data points per trajectory
    - delta: Time interval between points in the true trajectory
    - nu: Standard deviation of noise added to the sampled points
    - x0: Initial value of the process
    - ptype: Type of process ('OU', 'HL', or 'GBM')

    Returns:
    - orig_trajs: Array of true trajectories (n x K)
    - tgrid: Time points used for the trajectories
    - samp_trajs: Sampled noisy points for each trajectory
    - samp_ts: Corresponding time points for the sampled points
    """
    tgrid = np.arange(0, 1 + delta, delta)  # Time grid
    K = len(tgrid)
    orig_trajs = np.zeros((n, K))  # Initialize true trajectories
    orig_trajs[:, 0] = x0  # Set initial values

    # Generate true trajectories based on the specified process type
    for i in range(n):
        for j in range(1, K):
            z1 = np.random.normal()  # Standard normal random variable
            if ptype == "OU":
                orig_trajs[i, j] = ou_step(1, 1, delta, orig_trajs[i, j - 1], z1)
            elif ptype == "HL":
                orig_trajs[i, j] = hl_step(1, tgrid[j - 1], tgrid[j], orig_trajs[i, j - 1], z1)
            elif ptype == "GBM":
                orig_trajs[i, j] = gbm_step(0.3, 0.5, delta, orig_trajs[i, j - 1], z1)
            else:
                raise ValueError("Invalid process type. Choose 'OU', 'HL', or 'GBM'.")
    
    # Add Gaussian noise to the true trajectories
    noisy_trajs = orig_trajs + np.random.normal(scale=nu, size=(n, K))
    
    # Sample m consecutive points for each trajectory
    samp_trajs = []
    samp_ts = []
    for _ in range(n):
        start_idx = np.random.randint(0, K - m)
        samp_ts.append(tgrid[start_idx:start_idx + m])
        samp_trajs.append(noisy_trajs[i, start_idx:start_idx + m])
    
    # batching for sample trajectories is good for RNN; batching for original
    # trajectories only for ease of indexing
    orig_trajs = np.stack(orig_trajs, axis=0)
    samp_trajs = np.stack(samp_trajs, axis=0)
    samp_ts = np.stack(samp_ts, axis=0)    
    return orig_trajs, tgrid, samp_trajs, samp_ts


# Define the Latent ODE function
class LatentODEfunc(nn.Module):
    """
    Neural network representing the ODE function for the latent space dynamics.
    """
    def __init__(self, latent_dim=4, nhidden=20):
        super(LatentODEfunc, self).__init__()
        self.elu = nn.ELU(inplace=True)
        self.fc1 = nn.Linear(latent_dim, nhidden)
        self.fc2 = nn.Linear(nhidden, nhidden)
        self.fc3 = nn.Linear(nhidden, latent_dim)
        self.nfe = 0

    def forward(self, t, x):
        self.nfe += 1
        out = self.fc1(x)
        out = self.elu(out)
        out = self.fc2(out)
        out = self.elu(out)
        out = self.fc3(out)
        return out


# Recognition RNN for inferring initial conditions
class RecognitionRNN(nn.Module):
    """
    Recurrent neural network for recognizing the initial latent state from observations.
    """
    def __init__(self, latent_dim=4, obs_dim=2, nhidden=25, nbatch=1):
        super(RecognitionRNN, self).__init__()
        self.nhidden = nhidden
        self.nbatch = nbatch
        self.i2h = nn.Linear(obs_dim + nhidden, nhidden)
        self.h2o = nn.Linear(nhidden, latent_dim * 2)

    def forward(self, x, h):
        combined = torch.cat((x, h), dim=1)
        h = torch.tanh(self.i2h(combined))
        out = self.h2o(h)
        return out, h

    def initHidden(self):
        return torch.zeros(self.nbatch, self.nhidden)


# Decoder for mapping latent states to observations
class Decoder(nn.Module):
    """
    Neural network for decoding the latent states back to observed space.
    """
    def __init__(self, latent_dim=4, obs_dim=2, nhidden=20):
        super(Decoder, self).__init__()
        self.relu = nn.ReLU(inplace=True)
        self.fc1 = nn.Linear(latent_dim, nhidden)
        self.fc2 = nn.Linear(nhidden, obs_dim)

    def forward(self, z):
        out = self.fc1(z)
        out = self.relu(out)
        out = self.fc2(out)
        return out
    
    
# Running average meter for tracking metrics
class RunningAverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self, momentum=0.99):
        self.momentum = momentum
        self.reset()

    def reset(self):
        self.val = None
        self.avg = 0

    def update(self, val):
        if self.val is None:
            self.avg = val
        else:
            self.avg = self.avg * self.momentum + val * (1 - self.momentum)
        self.val = val


def log_normal_pdf(x, mean, logvar):
    """
    Compute the log probability density of a normal distribution.

    Parameters:
    - x: Input tensor (observed values)
    - mean: Mean of the normal distribution
    - logvar: Log variance of the normal distribution

    Returns:
    - Log probability density of the normal distribution
    """
    const = torch.from_numpy(np.array([2. * np.pi])).float().to(x.device)
    const = torch.log(const)
    return -.5 * (const + logvar + (x - mean) ** 2. / torch.exp(logvar))


def normal_kl(mu1, lv1, mu2, lv2):
    """
    Compute the KL divergence between two normal distributions.

    Parameters:
    - mu1, lv1: Mean and log variance of the first normal distribution
    - mu2, lv2: Mean and log variance of the second normal distribution

    Returns:
    - KL divergence between the two distributions
    """
    v1 = torch.exp(lv1)
    v2 = torch.exp(lv2)
    lstd1 = lv1 / 2.
    lstd2 = lv2 / 2.

    kl = lstd2 - lstd1 + ((v1 + (mu1 - mu2) ** 2.) / (2. * v2)) - .5
    return kl


def run_simulation(n, nu):
    """
    Run a single simulation of the latent ODE model.

    Parameters:
    - n: Number of trajectories (sample size)
    - nu: Noise level for the observations

    Returns:
    - RMSE (root-mean-square error) of the predicted trajectories
    """
    latent_dim = 4
    nhidden = 20
    rnn_nhidden = 25
    obs_dim = 1  # Observed dimension
    delta = 0.05
    device = torch.device('cuda:' + str(args.gpu) if torch.cuda.is_available() else 'cpu')

   # Generate data using the specified process type
    orig_trajs, tgrid, samp_trajs, samp_ts = generate_data(
        n=n, m=args.m, delta=delta, nu=nu, x0=args.x0, ptype = args.type
    )
    
    # Convert generated data to tensors
    orig_trajs = torch.from_numpy(orig_trajs).float().to(device)
    orig_ts = tgrid
    samp_trajs = torch.from_numpy(samp_trajs).float().to(device)
    samp_ts = torch.from_numpy(samp_ts).float().to(device)


    # model
    func = LatentODEfunc(latent_dim, nhidden).to(device)
    rec = RecognitionRNN(latent_dim, obs_dim, rnn_nhidden, n).to(device)
    dec = Decoder(latent_dim, obs_dim, nhidden).to(device)
    params = (list(func.parameters()) + list(dec.parameters()) + list(rec.parameters()))
    optimizer = optim.Adam(params, lr=args.lr)
    loss_meter = RunningAverageMeter()
    
    # Load checkpoint if available
    if args.train_dir is not None:
        if not os.path.exists(args.train_dir):
            os.makedirs(args.train_dir)
        ckpt_path = os.path.join(args.train_dir, 'ckpt.pth')
        if os.path.exists(ckpt_path):
            checkpoint = torch.load(ckpt_path)
            func.load_state_dict(checkpoint['func_state_dict'])
            rec.load_state_dict(checkpoint['rec_state_dict'])
            dec.load_state_dict(checkpoint['dec_state_dict'])
            optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            orig_trajs = checkpoint['orig_trajs']
            samp_trajs = checkpoint['samp_trajs']
            orig_ts = checkpoint['orig_ts']
            samp_ts = checkpoint['samp_ts']
            print('Loaded ckpt from {}'.format(ckpt_path))
            
    previous_loss = None  # To store the loss from the previous iteration
    threshold = 0.01  # Stop if loss decrease is less than 1%
    patience = 5  # Number of iterations to wait for improvement
    patience_counter = 0 
    
    # Training loop
    try:
        for itr in range(1, args.niters + 1):
            optimizer.zero_grad()
            
            # Infer initial latent state using RNN (backward pass)
            h = rec.initHidden().to(device)
            for t in reversed(range(samp_trajs.size(1))):
                obs = samp_trajs[:, t].unsqueeze(1)
                out, h = rec.forward(obs, h)
            qz0_mean, qz0_logvar = out[:, :latent_dim], out[:, latent_dim:]
            epsilon = torch.randn(qz0_mean.size()).to(device)
            z0 = epsilon * torch.exp(.5 * qz0_logvar) + qz0_mean

            # forward in time and solve ode for reconstructions
            results = []
            for i in range(samp_ts.size(0)):
                pred_z_i = odeint(func, z0[i], samp_ts[i])  # Result shape: [m, 4]
                results.append(pred_z_i)
            
            # Stack the results along a new dimension
            pred_z = torch.stack(results)  # Shape: [n, m, 4]

            pred_x = dec(pred_z).squeeze(-1)

            # Compute loss using negative ELBO
            noise_std_ = torch.zeros(pred_x.size()).to(device) + 0.1
            noise_logvar = 2. * torch.log(noise_std_).to(device)
            logpx = log_normal_pdf(
                samp_trajs, pred_x, noise_logvar).sum(-1).sum(-1)
            pz0_mean = pz0_logvar = torch.zeros(z0.size()).to(device)
            analytic_kl = normal_kl(qz0_mean, qz0_logvar,
                                    pz0_mean, pz0_logvar).sum(-1)
            loss = torch.mean(-logpx + analytic_kl, dim=0)
            loss.backward()
            optimizer.step()
            loss_meter.update(loss.item())
            print(f'n = {n}, nu = {nu}, Iter: {itr}, running avg elbo: {-loss_meter.avg:.4f}')
                
            # Early stopping based on loss improvement
            if previous_loss is not None:
                loss_diff = abs(previous_loss - loss_meter.avg)
                relative_decrease = loss_diff / max(abs(previous_loss), 1e-10)
                if relative_decrease < threshold:
                    patience_counter += 1
                else:
                    patience_counter = 0
                    
            # Stop training if patience limit is reached
            if patience_counter >= patience:
                print(f"n = {n}, nu = {nu}, Stopping early at iteration {itr} due to minimal running avg elbo decrease.")
                break        

            # Update the previous loss
            previous_loss = loss_meter.avg

    except KeyboardInterrupt:
        if args.train_dir is not None:
            ckpt_path = os.path.join(args.train_dir, 'ckpt.pth')
            torch.save({
                'func_state_dict': func.state_dict(),
                'rec_state_dict': rec.state_dict(),
                'dec_state_dict': dec.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'orig_trajs': orig_trajs,
                'samp_trajs': samp_trajs,
                'orig_ts': orig_ts,
                'samp_ts': samp_ts,
            }, ckpt_path)
            print('Stored ckpt at {}'.format(ckpt_path))

    # Evaluate the model
    with torch.no_grad():
        h = rec.initHidden().to(device)
        for t in reversed(range(samp_trajs.size(1))):
            obs = samp_trajs[:, t].unsqueeze(1)
            out, h = rec.forward(obs, h)
        qz0_mean, qz0_logvar = out[:, :latent_dim], out[:, latent_dim:]
        epsilon = torch.randn(qz0_mean.size()).to(device)
        z0 = epsilon * torch.exp(.5 * qz0_logvar) + qz0_mean
        orig_ts = torch.from_numpy(orig_ts).float().to(device)

        zs = odeint(func, z0, orig_ts).permute(1, 0, 2)
        pred_trajs = dec(zs).squeeze(-1)

    # Compute RMSE of the predictions
    orig_trajs = orig_trajs.cpu().numpy()  # Move to CPU and convert to numpy array
    pred_trajs = pred_trajs.cpu().numpy()
    K = len(tgrid)
    rmse = np.sqrt(np.mean((orig_trajs[:, K - 1] - pred_trajs[:, K - 1]) ** 2)).item()
    return rmse

def run_Q_sequential(n, nu, Q):
    """
    Run Q simulations for a given sample size (n) and noise level (nu) sequentially.

    Parameters:
    - n: Sample size (number of trajectories)
    - nu: Noise level for the observations
    - Q: Number of simulations to run

    Returns:
    - n: The input sample size
    - nu: The input noise level
    - mean_rmse: Mean RMSE (Root Mean Square Error) across the Q simulations
    - std_rmse: Standard deviation of the RMSE across the Q simulations
    """
    rmses = [run_simulation(n, nu) for _ in range(Q)]
    mean_rmse = np.mean(rmses)
    std_rmse = np.std(rmses)
    return n, nu, mean_rmse, std_rmse

def run_Q_parallel(n, nu, Q, num_cores):
    """
    Run Q simulations for a given sample size (n) and noise level (nu) in parallel.

    Parameters:
    - n: Sample size (number of trajectories)
    - nu: Noise level for the observations
    - Q: Number of simulations to run
    - num_cores: Number of CPU cores to use for parallel processing

    Returns:
    - n: The input sample size
    - nu: The input noise level
    - mean_rmse: Mean RMSE (Root Mean Square Error) across the Q simulations
    - std_rmse: Standard deviation of the RMSE across the Q simulations
    """
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [executor.submit(run_simulation, n, nu) for _ in range(Q)]
        rmses = [future.result() for future in concurrent.futures.as_completed(futures)]
    mean_rmse = np.mean(rmses)
    std_rmse = np.std(rmses)
    return n, nu, mean_rmse, std_rmse


# Main script block
if __name__ == '__main__':
    """
    Main script for running latent ODE simulations with varying sample sizes and noise levels.
    
    The script runs Q simulations for each combination of sample size (n) and noise level (nu)
    and saves the results to a text file. It can utilize parallel processing to speed up the simulations.
    """
    # Simulation parameters
    Q = 500  # Number of simulations to run for each (n, nu) combination
    nvec = [50, 200, 1000]  # List of sample sizes
    nuvec = [0, 0.01, 0.1]  # List of noise levels
    num_cores = 8  # Number of CPU cores to use for parallel processing
    
    # Start timer to measure total runtime
    start_time = time.time()

    # Open the results file once for writing
    filename = f'rmse_{args.type}.txt'
    with open(filename, 'w') as f:
        # Write the header line
        f.write("Sample Size\tNoise Level\tMean RMSE\tSD RMSE\n")
        
        # Run the simulations for each combination of n and nu in sequence
        for n in nvec:
            for nu in nuvec:
                # Run Q simulations in parallel for the current (n, nu) combination
                n, nu, mean_rmse, std_rmse = run_Q_parallel(n, nu, Q, num_cores)
                
                # Write the results for this (n, nu) pair to the file immediately
                f.write(f"{n}\t{nu}\t{mean_rmse:.4f}\t{std_rmse:.4f}\n")
    
    # End timer and calculate total runtime
    end_time = time.time()
    total_runtime = (end_time - start_time) / 60

    # Append the total runtime at the end of the file
    with open(filename, 'a') as f:
        f.write(f"\nTotal Runtime: {total_runtime:.2f} minutes\n")
    
    print(f"Results saved to {filename}")


# =============================================================================
# n = 1000
# # Start timer
# start_time = time.time()
# run_Q_sequential(n, m, Q)
# # End timer and print total runtime
# end_time = time.time()
# (end_time - start_time) / Q
# =============================================================================
