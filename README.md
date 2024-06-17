# Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets With Stochastic Differential Equations

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

| Data.files    | Details                                                          |
|:--------------|:-----------------------------------------------------------------|
| bgd.RData     | NSW experimental data, used in LaLonde (1986)                    |
| bgdplot.RData | Subset of NSW experimental data, used in Dehejia & Wahba (1999)  |
| bmd.RData     | CPS-SSA-1 controls, used in both papers                          |
| bmdplot.RData | PSID-1 controls, used in both papers                             |
| ngd.RData     | Data of lottery winners, used in Imbens, Rubin & Sacerdote (201) |
| ngdplot.RData | Reconstructed NSW AFDC female samples                            |
| ouplot.RData  | Reconstructed NSW AFDC female samples                            |

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

To report errors, please contact <ydzhou@ucdavis.edu>. Comments and
suggestions are welcome.

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

Zhou, Y. and Müller, H.G., 2023. Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets With Stochastic Differential Equations. arXiv preprint arXiv:2306.10221.

</div>

</div>

