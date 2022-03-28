
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spaSim

<!-- badges: start -->
<!-- badges: end -->

The goal of spaSim is to facilitate tissue image simulations! It
simulates cells with 2D locations and cell types in a tissue. The
available patterns include background cells, cell clusters, immune cell
rings and vessels. It also enables simulations that generate a set of
images in one run!

Simulations from spaSim can be applied to test and benchmark metrics
developed for spatial image analysis.

## Installation

You can install the development version of spaSim like so:

``` r
# enable this later
# install.packages("devtools")
# devtools::install_github("TrigosTeam/spaSim")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(spaSim)
set.seed(610)
mix_background <- TIS(n_cells = 5000, width = 2000, height = 2000, min_d = 10,
                      names_of_bg_cells = c("Tumour","Immune","Others"), 
                      proportions_of_bg_cells = c(0.1, 0.2, 0.7),
                      plot_image = TRUE)
```

<img src="man/figures/README-example-1.png" width="100%" /><img src="man/figures/README-example-2.png" width="100%" />

## Vignette

The vignette with an overview of the package can be accessed from the
top Menu under Articles [here](https://trigosteam.github.io/spaSim/).
