# Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets With Stochastic Differential Equations

This Github repo contains the data and codes necessary to replicate
**Zhou and Müller (2024)**: “Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets With Stochastic Differential Equations.” Journal of the Royal Statistical Society: Series B, just accepted

## Folder Structure

The folder structure of this repo is as follows:

| Folder      | Usage                                                                                                                                              |
|:------------|:---------------------------------------------------------------------------------------------------------------------------------------------------|
| code        | R scripts for the proposed approach                                                                                                                |
| data        | Data files                                                                                                                                         |
| mcfda       | R scripts for Lin and Wang (2022 JASA), adapted from https://github.com/linulysses/mcfda                                                           |
| node        | Python scripts for neural ordinary differential equations, adapted from https://github.com/rtqichen/torchdiffeq/blob/master/examples/latent_ode.py |
| sim         | R scripts for simulations                                                                                                                          |
| Figure1.R   | R script to replicate Figure 1                                                                                                                     |
| Figure2.R   | R script to replicate Figure 2                                                                                                                     |
| Figure3&5.R | R script to replicate Figure 3 and Figure 5                                                                                                        |
| Figure4.R   | R script to replicate Figure 4                                                                                                                     |

## code

R scripts implementing the proposed dynamic modeling approach.

| Data files | Usage                                                                   |
|:-----------|:------------------------------------------------------------------------|
| gdm.R      | Dynamic modeling of functional snippets with multiple linear regression |
| kerFctn.R  | Kernel function                                                         |
| ldm.R      | Dynamic modeling of functional snippets with local linear regression    |

## data

Zhou and Müller (2024) uses the following datasets, which are based on West Jr et al. (1997), Bachrach et al. (1999), and Tuddenham and Snyder (1954).

| Data files    | Details                                                                                       |
|:--------------|:----------------------------------------------------------------------------------------------|
| bgd.RData     | Berkeley growth study data (Tuddenham and Snyder 1954), also available in the R package `fda` |
| bgdplot.RData | Data to replicate Figure 4                                                                    |
| bmd.RData     | Spinal bone mineral density data (Bachrach et al. 1999)                                       |
| bmdplot.RData | Data to replicate Figure 3 and Figure 5                                                       |
| ngd.RData     | Nepal growth study data (West Jr et al. 1997)                                                 |
| ngdplot.RData | Data to replicate Figure 2                                                                    |
| ouplot.RData  | Data to replicate Figure 1                                                                    |

## node

Python scripts to replicate simulation results in subsection S.5.4 of the Supplementary Material.

| Data files       | Details                                                                                              |
|:-----------------|:-----------------------------------------------------------------------------------------------------|
| latentODE_bmd.py | Spinal bone mineral density data (Bachrach et al. 1999) using neural ordinary differential equations |
| latentODE_sim.py | Simulation using neural ordinary differential equations                                              |

## sim

R scripts to replicate simulation results in subsection 5.2 of the main text and Section S.5 of the Supplementary Material.

| Data files    | Details                                                                                |
|:--------------|:---------------------------------------------------------------------------------------|
| arma.R | Simulation using the autoregressive integrated moving average model                           |
| gbm.R  | Simulation for the geometric Brownian motion                                                  |
| hl.R   | Simulation for the Ho-Lee model                                                               |
| ou.R   | Simulation for the Ornstein-Uhlenbeck process                                                 |
| oucb.R | Simulation for uncertainty quantification using the Ornstein-Uhlenbeck process                |
| our.R  | Simulation for the case of five measurements per subject using the Ornstein-Uhlenbeck process |

## Report Errors

To report errors, please contact <ydzhou@ucdavis.edu>. Comments and suggestions are welcome.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-imbensrubinsacerdote" class="csl-entry">

Bachrach, L. K., Hastie, T., Wang, M.-C., Narasimhan, B. and Marcus, R. (1999) Bone mineral acquisition in healthy Asian, Hispanic, Black, and Caucasian youth: a longitudinal study. *Journal of Clinical Endocrinology & Metabolism*, **84**, 4702–4712.

</div>

<div id="ref-calonico2017women" class="csl-entry">

Tuddenham, R. D. and Snyder, M. M. (1954) Physical growth of California boys and girls from birth to eighteen years. *University of California Publications in Child Development*, **1**, 183–364.

</div>

<div id="ref-dehejiawahba" class="csl-entry">

West Jr, K. P., LeClerq, S. C., Shrestha, S. R., Wu, L. S.-F., Pradhan, E. K., Khatry, S. K., Katz, J., Adhikari, R. and Sommer, A. (1997) Effects of vitamin A on growth of vitamin A-deficient children: field studies in Nepal. *Journal of Nutrition*, **127**, 1957–1965.

</div>

<div id="ref-imbensxu" class="csl-entry">

Zhou, Y. and Müller, H.G., 2024. Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets With Stochastic Differential Equations. Journal of the Royal Statistical Society: Series B, just accepted.

</div>

</div>

