## code to prepare `D_shape1` dataset goes here
# D_shape1 -  ####
size <- 300
centre_loc_offset <- data.frame("x" = 1000, "y" = 1000)
infiltration_types <- c("Immune1", "Immune2", "Others")
infiltration_proportions <- c(0.15, 0.05, 0.05)
D_shape1 <- list(
  D1 = list(
    name_of_cluster_cell = "Tumour",
    size = size,
    shape = "Circle",
    centre_loc = data.frame("x" = 1000, "y" = 1000),
    infiltration_types = infiltration_types,
    infiltration_proportions = infiltration_proportions,
    name_of_ring_cell = "Immune1",
    immune_ring_width = 150,
    immune_ring_infiltration_types = c("Others"),
    immune_ring_infiltration_proportions = c(0.15),
    name_of_double_ring_cell = "Immune2",
    double_ring_width = 100,
    double_ring_infiltration_types = c("Others"),
    double_ring_infiltration_proportions = c(0.15)
  ),
  D2 = list(
    name_of_cluster_cell = "Tumour",
    size = size,
    shape = "Oval",
    centre_loc = data.frame("x" = 1200, "y" = 1200),
    infiltration_types = infiltration_types,
    infiltration_proportions = infiltration_proportions,
    name_of_ring_cell = "Immune1",
    immune_ring_width = 150,
    immune_ring_infiltration_types = c("Others"),
    immune_ring_infiltration_proportions = c(0.15),
    name_of_double_ring_cell = "Immune2",
    double_ring_width = 100,
    double_ring_infiltration_types = c("Others"),
    double_ring_infiltration_proportions = c(0.15)
  )
)

