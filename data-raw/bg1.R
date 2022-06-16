## code to prepare `bg1` dataset goes here
n_cells <- 5000
width <- 2000
height <- 2000
min_d <- 10

# set seed for this background image simulation for reproducibility
set.seed(610)
bg1 <- simulate_background_cells(n_cells = n_cells, width = width, height = height,
                                 min_d = min_d, oversampling_rate = 1.5)

usethis::use_data(bg1, overwrite = T)
