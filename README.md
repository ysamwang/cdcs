
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The cdcs package

<!-- badges: start -->
<!-- badges: end -->

This R package contains source code and replication files for the paper
[Confidence Sets for Causal Orderings](https://arxiv.org/abs/2305.14506)
by Wang, Kolar, and Drton.

## Installation

You can install cdcs from [GitHub](https://github.com/) with:

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

- To replicate Table 1 use the script table1_replication.R
- To replicate Figure 2 use the script figure2_replication.R
- To replicate Table 2 use the script table3_replication.R
- To replicate Figure 3 use the script figure3_replication.R
- To replicate the real data analysis
