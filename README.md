# Correcting for Missing Data When Evaluating Surrogate Markers in a Clinical Trial

This repository contains `R` code and simulation data to reproduce results from the manuscript by [Lotspeich, Nguyen, and Parast (2026+)](https://arxiv.org/abs/). 

These simulations rely on the `missSurrogate` package, which implements the various estimators from the paper. The package can be found in its own repo [here](https://github.com/sarahlotspeich/missSurrogate) and installed in `R` as follows:

``` r
# Run once
# install.packages("devtools")
devtools::install_github("sarahlotspeich/missSurrogate", ref = "main")
```

## Tables of Simulation Results

**Table 1.** Simulation results for $\widehat{R}_S$, the proportion of treatment effect explained (PTE) in Setting 1, where the surrogate marker $S$ was missing completely at random for $32-39\%$ of patients in the sample of $n = 2000$. With $R_S = 0.5$, the surrogate $S$ was of moderate strength.

  - [Script (Run Simulations)](simulations/sett1_mcar.R)
  - [Script (Make Table)](tables/sett1_mcar.R)
  - [Data (Simulation Results)](simulations/sett1_mcar.csv)

**Table 2 (top).** Simulation results for $\widehat{R}_S$, the proportion of treatment effect explained (PTE) in Settings 2 and 3, when the surrogate marker $S$ was missing at random depending on treatment $Z$ (top portion) and depending on the primary outcome $Y$ (bottom portion) in the sample of $n = 2000$. With $R_S = 0.5$, the surrogate $S$ was of moderate strength.

  - [Script (Run Simulations)](simulations/sett2_mar_givY.R)
  - [Script (Make Table)](tables/sett2_mar_givY.R)
  - [Data (Simulation Results)](simulations/sett2_mar_givY.csv)

**Table 2 (bottom).** Simulation results for $\widehat{R}_S$, the proportion of treatment effect explained (PTE) in Settings 2 and 3, when the surrogate marker $S$ was missing at random depending on treatment $Z$ (top portion) and depending on the primary outcome $Y$ (bottom portion) in the sample of $n = 2000$. With $R_S = 0.5$, the surrogate $S$ was of moderate strength.

  - [Script (Run Simulations)](simulations/sett3_mar_givZ.R)
  - [Script (Make Table)](tables/sett3_mar_givZ.R)
  - [Data (Simulation Results)](simulations/sett3_mar_givZ.csv)

**Table 3.** Simulation results for $\widehat{R}_S$, the proportion of treatment effect explained (PTE) in Setting 5, when the distributions of the surrogate markers $S$ do not overlap between the treatment and control groups. 

  - [Script (Run Simulations)](simulations/sett5_non_overlap.R)
  - [Script (Make Table)](tables/sett5_non_overlap.R)
  - [Data (Simulation Results)](simulations/sett5_non_overlap.csv)

## Figures of Simulation Results
**Figure 1.** The nonparametric estimator $\widehat{R}_S$, of the proportion of treatment effect explained (PTE), assumes that the distributions of the surrogate marker $S$ in the treatment and control groups overlap, as in **A)**. In Setting 5, we consider the setting where these distributions do not fully overlap, as in **B)**.

  - [Figure](figures/surrogate_density_overlap.pdf)
  - [Script (Make Figure)](figures/surrogate_density_overlap.R)

**Figure 2.** Simulation results for $\widehat{R}_S$, the proportion of treatment effect explained (PTE) in Setting 4, when the surrogate marker $S$ was missing at random given primary outcome $Y$, treatment group $Z$, and their interaction $Y \times Z$. The weights model for the IPW approaches is misspecified when it includes $Y$ only or $Z$ only. 

  - [Figure](figures/sett5_non_overlap_boxplot.pdf)
  - [Script (Make Figure)](figures/sett5_non_overlap_boxplot.R)
  - [Data (Simulation Results)](simulations/sett5_non_overlap.csv)
