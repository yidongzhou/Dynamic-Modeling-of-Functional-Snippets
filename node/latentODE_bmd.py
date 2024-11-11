# Import necessary libraries
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('agg') # Use 'agg' backend for saving plots without displaying
import torch
import torch.nn as nn
import torch.optim as optim


# Set up argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('--adjoint', type=eval, default=False)  # Use adjoint method for memory efficiency
parser.add_argument('--niters', type=int, default=500)  # Number of training iterations
parser.add_argument('--lr', type=float, default=0.01)  # Learning rate for optimizer
parser.add_argument('--gpu', type=int, default=0)  # GPU index (if available)
parser.add_argument('--train_dir', type=str, default=None)  # Directory for saving checkpoints
args = parser.parse_args()


# Import the appropriate ODE solver
if args.adjoint:
    from torchdiffeq import odeint_adjoint as odeint  # Adjoint method for memory efficiency
else:
    from torchdiffeq import odeint

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


# Main script block
if __name__ == '__main__':
    latent_dim = 4
    nhidden = 20
    rnn_nhidden = 25
    obs_dim = 1
    device = torch.device('cuda:' + str(args.gpu) if torch.cuda.is_available() else 'cpu')
    
    # Load the .txt file for female
    yf = np.loadtxt("bmd/yf.txt")
    tf = np.loadtxt("bmd/tf.txt")
    n = yf.shape[0]
    M = n
    
    tgrid = np.concatenate(([10.1], np.arange(11, 25)))
    x0 = 0.778
    
    # Convert data to tensors
    samp_trajs = torch.from_numpy(yf).float().to(device)
    samp_ts = torch.from_numpy(tf).float().to(device)
    tgrid = torch.from_numpy(tgrid).float().to(device)
    
    # model
    func = LatentODEfunc(latent_dim, nhidden).to(device)
    rec = RecognitionRNN(latent_dim, obs_dim, rnn_nhidden, n).to(device)
    dec = Decoder(latent_dim, obs_dim, nhidden).to(device)
    params = (list(func.parameters()) + list(dec.parameters()) + list(rec.parameters()))
    optimizer = optim.Adam(params, lr=args.lr)
    loss_meter = RunningAverageMeter()
    
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
            samp_trajs = checkpoint['samp_trajs']
            samp_ts = checkpoint['samp_ts']
            print('Loaded ckpt from {}'.format(ckpt_path))
            
    previous_loss = None  # To store the loss from the previous iteration
    threshold = 0.01  # Stop if loss decrease is less than 1%
    patience = 5  # Number of iterations to wait for improvement
    patience_counter = 0
       
    try:
        for itr in range(1, args.niters + 1):
            optimizer.zero_grad()
            # backward in time to infer q(z_0)
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
    
            # compute loss
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
            print(f'Iter: {itr}, running avg elbo: {-loss_meter.avg:.4f}')
    
    except KeyboardInterrupt:
        if args.train_dir is not None:
            ckpt_path = os.path.join(args.train_dir, 'ckpt.pth')
            torch.save({
                'func_state_dict': func.state_dict(),
                'rec_state_dict': rec.state_dict(),
                'dec_state_dict': dec.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'samp_trajs': samp_trajs,
                'samp_ts': samp_ts,
            }, ckpt_path)
            print('Stored ckpt at {}'.format(ckpt_path))
    
    with torch.no_grad():
        # sample from a given initial value
        h = rec.initHidden().to(device)
        obs = torch.from_numpy(np.full(M, x0)).float().to(device).unsqueeze(1)
        out, h = rec.forward(obs, h)
        qz0_mean, qz0_logvar = out[:, :latent_dim], out[:, latent_dim:]
        epsilon = torch.randn(qz0_mean.size()).to(device)
        z0 = epsilon * torch.exp(.5 * qz0_logvar) + qz0_mean
    
        zs = odeint(func, z0, tgrid).permute(1, 0, 2)
        pred_trajs = dec(zs).squeeze(-1)
        
    pred_trajs = pred_trajs.cpu().numpy()  # Move to CPU and convert to numpy array
    np.savetxt('bmd/bmdf.txt', pred_trajs, delimiter='\t', fmt='%f')
        
    # Load the .txt file for male
    ym = np.loadtxt("bmd/ym.txt")
    tm = np.loadtxt("bmd/tm.txt")
    n = ym.shape[0]
    M = n
    
    tgrid = np.arange(9, 25)
    x0 = 0.642
    
    # Convert data to tensors
    samp_trajs = torch.from_numpy(ym).float().to(device)
    samp_ts = torch.from_numpy(tm).float().to(device)
    tgrid = torch.from_numpy(tgrid).float().to(device)
    
    # model
    func = LatentODEfunc(latent_dim, nhidden).to(device)
    rec = RecognitionRNN(latent_dim, obs_dim, rnn_nhidden, n).to(device)
    dec = Decoder(latent_dim, obs_dim, nhidden).to(device)
    params = (list(func.parameters()) + list(dec.parameters()) + list(rec.parameters()))
    optimizer = optim.Adam(params, lr=args.lr)
    loss_meter = RunningAverageMeter()
    
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
            samp_trajs = checkpoint['samp_trajs']
            samp_ts = checkpoint['samp_ts']
            print('Loaded ckpt from {}'.format(ckpt_path))
            
    previous_loss = None  # To store the loss from the previous iteration
    threshold = 0.01  # Stop if loss decrease is less than 1%
    patience = 5  # Number of iterations to wait for improvement
    patience_counter = 0
        
    try:
        for itr in range(1, args.niters + 1):
            optimizer.zero_grad()
            # backward in time to infer q(z_0)
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
    
            # compute loss
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
            print(f'Iter: {itr}, running avg elbo: {-loss_meter.avg:.4f}')
    
    except KeyboardInterrupt:
        if args.train_dir is not None:
            ckpt_path = os.path.join(args.train_dir, 'ckpt.pth')
            torch.save({
                'func_state_dict': func.state_dict(),
                'rec_state_dict': rec.state_dict(),
                'dec_state_dict': dec.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'samp_trajs': samp_trajs,
                'samp_ts': samp_ts,
            }, ckpt_path)
            print('Stored ckpt at {}'.format(ckpt_path))
    
    with torch.no_grad():
        # sample from a given initial value
        h = rec.initHidden().to(device)
        obs = torch.from_numpy(np.full(M, x0)).float().to(device).unsqueeze(1)
        out, h = rec.forward(obs, h)
        qz0_mean, qz0_logvar = out[:, :latent_dim], out[:, latent_dim:]
        epsilon = torch.randn(qz0_mean.size()).to(device)
        z0 = epsilon * torch.exp(.5 * qz0_logvar) + qz0_mean
    
        zs = odeint(func, z0, tgrid).permute(1, 0, 2)
        pred_trajs = dec(zs).squeeze(-1)
    
    pred_trajs = pred_trajs.cpu().numpy()  # Move to CPU and convert to numpy array
    np.savetxt('bmd/bmdm.txt', pred_trajs, delimiter='\t', fmt='%f')