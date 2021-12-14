simulate_immune_rings <- function(background_sample,
                                  bg_type = "Others",
                                  n_immune_rings = 2,
                                  win = NULL,
                                  properties_of_immune_rings = list(
                                    I1 = list(
                                      name_of_cluster_cell = "Tumour",
                                      size = 700,
                                      shape = "Circle",
                                      centre_loc = data.frame("x" = 1800, "y" = 1500),
                                      infiltration_types = c("Immune1", "Immune2", "Others"),
                                      infiltration_proportions = c(0.15, 0.05, 0.05),
                                      name_of_ring_cell = "Immune1",
                                      immune_ring_width = 150,
                                      immune_ring_infiltration_types = c("Others"),
                                      immune_ring_infiltration_proportions = c(0.15)
                                    ),
                                    I2 = list(
                                      name_of_cluster_cell = "Tumour",
                                      size = 700,
                                      shape = "Circle",
                                      centre_loc = data.frame("x" = 1800, "y" = 1500),
                                      infiltration_types = c("Immune1", "Immune2", "Others"),
                                      infiltration_proportions = c(0.15, 0.05, 0.05),
                                      name_of_ring_cell = "Immune1",
                                      immune_ring_width = 150,
                                      immune_ring_infiltration_types = c("Others"),
                                      immune_ring_infiltration_proportions = c(0.15)
                                    )
                                  )

) {

  ## CHECK
  # check if the specified cluster properties match n_immune_rings
  if (as.numeric(length(properties_of_immune_rings)) != n_immune_rings){
    stop("`n_immune_rings` does not match the length of `properties_of_immune_rings`!")
  }

  # Get the window
  # if window is specified, use the specified window
  # otherwise, use the window of the background sample
  if (is.null(win)) {
    X <- max(background_sample$Cell.X.Position)
    Y <- max(background_sample$Cell.Y.Position)
    win <- owin(c(0, X), c(0,Y))
  }

  # Default phenotype is specified by bg_type
  # (when background sample does not have Phenotype)
  if (is.null(background_sample$Phenotype)){
    background_sample[, "Phenotype"] <- bg_type
  }

  n_cells <- dim(background_sample)[1]

  for (k in 1:n_immune_rings) { # for each cluster
    # get the arguments

    cluster_cell_type <- properties_of_immune_rings[[k]]$name_of_cluster_cell
    size <- properties_of_immune_rings[[k]]$size
    shape <- properties_of_immune_rings[[k]]$shape
    centre_loc <- properties_of_immune_rings[[k]]$centre_loc
    infiltration_types <- properties_of_immune_rings[[k]]$infiltration_types
    infiltration_proportions <- properties_of_immune_rings[[k]]$infiltration_proportions
    ring_cell_type = properties_of_immune_rings[[k]]$name_of_ring_cell
    ring_width = properties_of_immune_rings[[k]]$immune_ring_width
    ring_infiltration_types = properties_of_immune_rings[[k]]$immune_ring_infiltration_types
    ring_infiltration_proportions = properties_of_immune_rings[[k]]$immune_ring_infiltration_proportions

    # generate a location as the centre of the cluster
    if (is.null(centre_loc)){
      seed_point <- runifpoint(1, win=win)}
    else seed_point <- centre_loc
    a <- seed_point$x
    b <- seed_point$y

    # cluster size is the radius of the cluster
    r <- size
    R <- r^2

    # cluster shape
    shape <- shape
    Circle <- (shape == "Circle")
    Oval <- (shape == "Oval")

    # immune ring radius
    I_R <- (r+ring_width)^2

    # determine if each cell is in the cluster or in the immune ring or neither
    for (i in 1:n_cells){
      x <- background_sample[i, "Cell.X.Position"]
      y <- background_sample[i, "Cell.Y.Position"]
      pheno <- background_sample[i, "Phenotype"]

      # squared distance to the cluster centre
      A <- (x - a)^2
      B <- (y - b)^2
      AB <- (x-a)*(y-b)
      D <- Circle*(A + B) + Oval*(A + AB + B)

      # determine which region the point falls in
      if (D < R){
        # generate random number to decide the phenotype
        random <- runif(1)
        n_infiltration_types <- length(infiltration_types)
        pheno <- cluster_cell_type
        n <- 1
        current_p <- 0
        while (n <= n_infiltration_types){
          current_p <- current_p + infiltration_proportions[n]
          if (random <= current_p) {
            pheno <- infiltration_types[n]
            break
          }
          n <- n+1
        }
      }

      else if(D < I_R){
        # generate random number to decide the phenotype
        random <- runif(1)
        n_ring_infiltration_types <- length(ring_infiltration_types)
        pheno <- ring_cell_type
        n <- 1
        current_p <- 0
        while (n <= n_ring_infiltration_types){
          current_p <- current_p + ring_infiltration_proportions[n]
          if (random <= current_p) {
            pheno <- ring_infiltration_types[n]
            break
          }
          n <- n+1
        }
      }

      if (background_sample[i, "Phenotype"] != cluster_cell_type){
        background_sample[i, "Phenotype"] <- pheno }
    }
  }

  return(background_sample)
}

