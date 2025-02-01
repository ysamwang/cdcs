
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The cdcs package

<!-- badges: start -->
<!-- badges: end -->

This R package contains source code and replication files for the paper
[Confidence Sets for Causal Orderings](https://arxiv.org/abs/2305.14506)
by Wang, Kolar, and Drton. Specifically, the project seeks to quantify
uncertainty for causal discovery tasks. The primary task is to create a
confidence set for causal orderings—a total ordering of the variables in
which ancestors precede descendants. This set can be post-processed to
form a lower/upper envelope of ancestral relations and confidence
intervals for causal effects which incorporate uncertainty about the
causal structure.

## Installation

You can install the R package (the acronym is for Causal Discovery
Confidence Sets) from it’s [GitHub
repository](https://github.com/ysamwang/cdcs) with the following code:

``` r
# install.packages("pak")
pak::pak("ysamwang/cdcs")
```

## Replicating the numerical results in the paper

The simulations provided in the paper are quite extensive and would
require a high performance cluster to replicate the exact settings used
in the paper. Nonetheless, we provide the code used to generate
numerical results and under settings where the problem size and number
of replications can run on most personal machines.

- To replicate the simulations used to create Table 1 use the script
  table1_replication.R
- To replicate the simulations used to create Figure 2 use the script
  figure2_replication.R
- To replicate the simulations used to create Table 2 use the script
  table3_replication.R
- To replicate the simulations used to create Figure 3 use the script
  figure3_replication.R
- To replicate the real data analysis use the script
  dataAnalysis_replication.R. The data is the 12 portfolio average value
  weighted data for 2019-2023 from [Kenneth French’s data
  library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/12_Industry_Portfolios_daily_CSV.zip)
  and is located in the data folder. It was accessed in Mar 2024. Note:
  this may take roughly 24 hrs to run on a personal computer.
- To replicate the simulations in the appendix use the script
  appendix_replication.R
