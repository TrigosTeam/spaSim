## code to prepare `bg1` dataset goes here
n_cells <- 5000
width <- 2000
height <- 2000
min_d <- 10
names_of_bg_cells <- c("Tumour", "Immune", "Others")
proportions_of_bg_cells <- c(0.4, 0.4, 0.2)
bg1 <- simulate_background_cells(n_cells = n_cells, width = width, height = height, min_d = min_d, oversample = 1.5)
