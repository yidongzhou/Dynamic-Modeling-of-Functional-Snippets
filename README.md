# Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets With Stochastic Differential Equations
================

This Github repo contains the data and codes necessary to replicate
**Zhou and Müller (2024+)**: “Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets With Stochastic Differential Equations.” [\[arXiv\]](https://arxiv.org/pdf/2306.10221)

## Folder Structure

The folder structure of this repo is as follows:

| folder | usage                                                                                    |
|:-------|:-----------------------------------------------------------------------------------------|
| code   | R scripts for the proposed approach                                                      |
| data   | Data files                                                                               |
| mcfda  | R scripts for Lin and Wang (2022 JASA), adapted from https://github.com/linulysses/mcfda |
| sim    | R scripts for simulations                                                                |

## Data Files

Zhou and Müller (2024+) uses the following datasets, which are based on
West Jr et al. (1997), Bachrach et al. (1999), and Tuddenham and Snyder (1954).

| Data.files    | Details                                                          | File_Type | Experimental |
|:--------------|:-----------------------------------------------------------------|:----------|:-------------|
| bgd.RData     | NSW experimental data, used in LaLonde (1986)                    | Stata     | Yes          |
| bgdplot.RData | Subset of NSW experimental data, used in Dehejia & Wahba (1999)  | Stata     | Yes          |
| bmd.RData     | CPS-SSA-1 controls, used in both papers                          | Stata     | No           |
| bmdplot.RData | PSID-1 controls, used in both papers                             | Stata     | No           |
| ngd.RData     | Data of lottery winners, used in Imbens, Rubin & Sacerdote (201) | R         | No           |
| ngdplot.RData | Reconstructed NSW AFDC female samples                            | Stata     | Both         |
| ouplot.RData  | Reconstructed NSW AFDC female samples                            | Stata     | Both         |

## R Scripts

To replicate all findings, set the directory to the root folder of the
replications files and execute `master.R`. Below are explanations for
the usage of each R script:

| Data.files          | Usage                                                  |
|:--------------------|:-------------------------------------------------------|
| master.R            | Install necessary R packages and excute all R scripts. |
| functions_est.R     | Store functions for estimation                         |
| functions_plot.R    | Store functions for making plots                       |
| lalonde1_prepare.R  | Preprocess LaLonde datasets                            |
| lalonde2_trim.R     | Trim datasets to improve overlap                       |
| lalonde3_overlap.R  | Visualize overlap in propensity scores                 |
| lalonde4_estimate.R | Estimate the ATT                                       |
| lalonde5_catt.R     | Estimate and visualize CATT                            |
| lalonde6_qte.R      | Estimate and visualize quantile treatment effects      |
| laldone7_sens.R     | Conduct sensitivity analyses                           |
| lalonde8_lcs.R      | Analyze the female samples by Calonico & Smith (2017)  |
| irs1_est.R          | Estimate the ATT using the IRS data                    |
| irs2_big.R          | Additional analyses for winning big prizes             |
| irs3_small.R        | Additional analyses for winning small prizes           |

## Install R Packages

For successful replication, the following R packages need to be
installed:

``` r
# required packages
packages <- c("haven", "labelled", "Matching", "grf", "sensemakr", "qte",
    "estimatr", "CBPS", "hbal", "DoubleML", "mlr3learners", "fixest", "ggplot2")

# install packages
install_all <- function(packages) {
  installed_pkgs <- installed.packages()[, "Package"]
  for (pkg in packages) {
    if (!pkg %in% installed_pkgs) {
      install.packages(pkg)
    }
  }
}
install_all(packages)
```

## Report Errors

To report errors, please contact <yiqingxu@stanford.edu>. Comments and
suggestions are welcome.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-calonico2017women" class="csl-entry">

Calónico, Sebastian, and Jeffrey Smith. 2017. “The Women of the National
Supported Work Demonstration.” *Journal of Labor Economics* 35 (S1):
S65–97.

</div>

<div id="ref-dehejiawahba" class="csl-entry">

Dehejia, Rajeev H, and Sadek Wahba. 1999. “Causal Effects in
Nonexperimental Studies: Reevaluating the Evaluation of Training
Programs.” *Journal of the American Statistical Association* 94 (448):
1053–62.

</div>

<div id="ref-imbensrubinsacerdote" class="csl-entry">

Imbens, Guido W, Donald B Rubin, and Bruce I Sacerdote. 2001.
“Estimating the Effect of Unearned Income on Labor Earnings, Savings,
and Consumption: Evidence from a Survey of Lottery Players.” *American
Economic Review* 91 (4): 778–94.

</div>

<div id="ref-imbensxu" class="csl-entry">

Imbens, Guido W, and Yiqing Xu. 2024. “LaLonde (1986) After Nearly Four
Decades: Lessons Learned.” arXiv:2406.00827.

</div>

<div id="ref-LaLonde" class="csl-entry">

LaLonde, Robert J. 1986. “Evaluating the Econometric Evaluations of
Training Programs with Experimental Data.” *The American Economic
Review* 76 (4): 604–20.

</div>

</div>

