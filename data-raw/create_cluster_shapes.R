## code to prepare `C_shape1`, `C_shape2`, `C_shape3` dataset goes here
# C_shape1 - Tumour cluster ####
centre_loc_offset <- data.frame("x" = 1000, "y" = 1000)
infiltration_types <- c("Immune", "Others")
infiltration_proportions <- c(0.1, 0.1)
C_shape1 <- list(
  C1 = list(
    name_of_cluster_cell = "Tumour",
    size = 100,
    shape = "Circle",
    centre_loc = centre_loc_offset,
    infiltration_types = infiltration_types,
    infiltration_proportions = infiltration_proportions
  ),
  C2 = list(
    name_of_cluster_cell = "Tumour",
    size = 0,
    shape = "Oval",
    centre_loc = centre_loc_offset + c(450, 250),
    infiltration_types = infiltration_types,
    infiltration_proportions = infiltration_proportions
  ),
  C3 = list(
    name_of_cluster_cell = "Tumour",
    size = 100,
    shape = "Oval",
    centre_loc = centre_loc_offset + c(-250, -300),
    infiltration_types = infiltration_types,
    infiltration_proportions = infiltration_proportions
  )
)


# C_shape2 - Tumour cluster ####
centre_loc_offset <- data.frame("x" = 1000, "y" = 1000)
infiltration_types <- c("Immune", "Others")
infiltration_proportions <- c(0.1, 0.1)
C_shape2 <- list(
  C1 = list(
    name_of_cluster_cell = "Tumour",
    size = 100,
    shape = "Circle",
    centre_loc = centre_loc_offset + c(0, -200),
    infiltration_types = infiltration_types,
    infiltration_proportions = infiltration_proportions
  ),
  C2 = list(
    name_of_cluster_cell = "Tumour",
    size = 50,
    shape = "Oval",
    centre_loc = centre_loc_offset + c(100, 50),
    infiltration_types = infiltration_types,
    infiltration_proportions = infiltration_proportions
  )
)

# C_shape3 - Immune cluster ####
centre_loc_offset <- data.frame("x" = 1000, "y" = 1000)
infiltration_types <- c("Immune", "Others")
infiltration_proportions <- c(0.1, 0.1)
C_shape3 <- list(
  C1 = list(
    name_of_cluster_cell = "Immune",
    size = 100,
    shape = "Irregular",
    centre_loc = centre_loc_offset,
    infiltration_types = infiltration_types,
    infiltration_proportions = infiltration_proportions
  ))
