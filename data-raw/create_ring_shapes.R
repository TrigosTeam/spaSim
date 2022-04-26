## code to prepare `R_shape1`, `R_shape2` and `R_shape3` dataset goes here
# R_shape1 - Tumour cluster and immune ring ####
centre_loc_offset <- data.frame("x" = 1000, "y" = 1000)
infiltration_types <- c("Immune", "Others")
infiltration_proportions <- c(0, 0.1)
R_shape1 <- list(
    I1 = list(
        name_of_cluster_cell = "Tumour",
        size = 100,
        shape = "Circle",
        centre_loc = centre_loc_offset,
        infiltration_types = infiltration_types,
        infiltration_proportions = infiltration_proportions,
        name_of_ring_cell = "Immune",
        immune_ring_width = 50,
        immune_ring_infiltration_types = c("Others"),
        immune_ring_infiltration_proportions = c(0.15)
    ),
    I2 = list(
        name_of_cluster_cell = "Tumour",
        size = 0,
        shape = "Oval",
        centre_loc = centre_loc_offset + c(450, 250),
        infiltration_types = infiltration_types,
        infiltration_proportions = infiltration_proportions,
        name_of_ring_cell = "Immune",
        immune_ring_width = 50,
        immune_ring_infiltration_types = c("Others"),
        immune_ring_infiltration_proportions = c(0.15)
    ),
    I3 = list(
        name_of_cluster_cell = "Tumour",
        size = 100,
        shape = "Oval",
        centre_loc = centre_loc_offset + c(-250, -300),
        infiltration_types = infiltration_types,
        infiltration_proportions = infiltration_proportions,
        name_of_ring_cell = "Immune",
        immune_ring_width = 50,
        immune_ring_infiltration_types = c("Others"),
        immune_ring_infiltration_proportions = c(0.15)
    )
)

# R_shape2 - Tumour cluster and immune ring ####
centre_loc_offset <- data.frame("x" = 1000, "y" = 1000)
infiltration_types <- c("Immune", "Others")
infiltration_proportions <- c(0, 0.1)

R_shape2 = list(
    I1 = list(
        name_of_cluster_cell = "Tumour",
        size = 100,
        shape = "Circle",
        centre_loc = centre_loc_offset + c(0, -200),
        infiltration_types = infiltration_types,
        infiltration_proportions = infiltration_proportions,
        name_of_ring_cell = "Immune",
        immune_ring_width = 50,
        immune_ring_infiltration_types = c("Others"),
        immune_ring_infiltration_proportions = c(0.15)
    ),
    I2 = list(
        name_of_cluster_cell = "Tumour",
        size = 50,
        shape = "Oval",
        centre_loc = centre_loc_offset + c(100, 50),
        infiltration_types = infiltration_types,
        infiltration_proportions = infiltration_proportions,
        name_of_ring_cell = "Immune",
        immune_ring_width = 50,
        immune_ring_infiltration_types = c("Others"),
        immune_ring_infiltration_proportions = c(0.15)
    )
)


# R_shape3 - Tumour cluster and immune ring ####
centre_loc_offset <- data.frame("x" = 1000, "y" = 1000)
infiltration_types <- c("Immune", "Others")
infiltration_proportions <- c(0, 0.1)

R_shape3 = list(
    I1 = list(
        name_of_cluster_cell = "Tumour",
        size = 100,
        shape = "Circle",
        centre_loc = centre_loc_offset + c(0, -200),
        infiltration_types = infiltration_types,
        infiltration_proportions = infiltration_proportions,
        name_of_ring_cell = "Immune",
        immune_ring_width = 50,
        immune_ring_infiltration_types = c("Tumour", "Others"),
        immune_ring_infiltration_proportions = c(0.3, 0.15)
    ),
    I2 = list(
        name_of_cluster_cell = "Tumour",
        size = 50,
        shape = "Oval",
        centre_loc = centre_loc_offset + c(100, 50),
        infiltration_types = infiltration_types,
        infiltration_proportions = infiltration_proportions,
        name_of_ring_cell = "Immune",
        immune_ring_width = 50,
        immune_ring_infiltration_types = c("Tumour","Others"),
        immune_ring_infiltration_proportions = c(0.3, 0.15)
    )
)
