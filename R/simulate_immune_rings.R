#' simulate_immune_rings
#'
#' @description Based on an existing background image, simulate rings of immune
#'   cells that surround tumour clusters. The tumour clusters and immune rings
#'   are simulated at the same time.
#'
#' @param background_sample (OPTIONAL) Data.Frame. The image that the immune
#'   rings are simulated on. By default use the internal background image `bg1`.
#' @param bg_type (OPTIONAL) String The name of the background cell type. By
#'   default is "Others".
#' @param n_immune_rings Number of immune rings. This must match the
#'   `length(properties_of_immune_rings)`.
#' @param win (OPTIONAL) owin object output from spatstat.geom::owin function.
#'   By default is the window of the background image
#' @param properties_of_immune_rings List of properties of the immune rings
#'
#' @family simulate pattern functions
#' @seealso   \code{\link{simulate_background_cells}} for all cell simulation,
#'   \code{\link{simulate_mixing}} for mixed background simulation,
#'   \code{\link{simulate_clusters}} for cluster simulation,
#'   \code{\link{simulate_double_rings}} for double immune ring simulation, and
#'   \code{\link{simulate_stripes}} for vessel simulation.
#'
#' @return A data.frame of the simulated image
#' @export
#'
#' @examples
#' set.seed(610)
#' # manually define the properties of the immune ring
#' properties_of_immune_rings <- list(I1 = list(name_of_cluster_cell = "Tumour",
#' size = 600,shape = "Circle",centre_loc = data.frame("x" = 930, "y" = 1000),
#' infiltration_types = c("Immune1", "Immune2", "Others"), infiltration_proportions
#' = c(0.15, 0.05, 0.05), name_of_ring_cell = "Immune1", immune_ring_width = 150,
#' immune_ring_infiltration_types = c("Others"), immune_ring_infiltration_proportions = c(0.15)))
#' # simulate immune rings (`n_immune_rings` should match the length of `properties_of_immune_rings`)
#' immune_ring_image <- simulate_immune_rings(background_sample = bg1,
#' n_immune_rings = 1, properties_of_immune_rings = properties_of_immune_rings)
#' # library(SPIAT)
#' # plot_cell_categories(immune_ring_image, c("Tumour","Immune1"),c("red","blue"),
#' # "Phenotype")
#'
simulate_immune_rings <- function(background_sample = bg1,
                                  bg_type = "Others",
                                  n_immune_rings = 2,
                                  win = NULL,
                                  properties_of_immune_rings = list(
                                    I1 = list(
                                      name_of_cluster_cell = "Tumour",
                                      size = 600,
                                      shape = "Circle",
                                      centre_loc = data.frame("x" = 930, "y" = 1000),
                                      infiltration_types = c("Immune1", "Immune2", "Others"),
                                      infiltration_proportions = c(0.15, 0.05, 0.05),
                                      name_of_ring_cell = "Immune1",
                                      immune_ring_width = 150,
                                      immune_ring_infiltration_types = c("Others"),
                                      immune_ring_infiltration_proportions = c(0.15)
                                    ),
                                    I2 = list(
                                      name_of_cluster_cell = "Tumour",
                                      size = 500,
                                      shape = "Oval",
                                      centre_loc = data.frame("x" = 1330, "y" = 1100),
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
  # is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SummarizedExperiment::colData(background_sample))}

  # check if the specified cluster properties match n_immune_rings
  if (as.numeric(length(properties_of_immune_rings)) != n_immune_rings){
    stop("`n_immune_rings` does not match the length of `properties_of_immune_rings`!")
  }

  # add a new column to store the position label for each cell (0 for core cluster,
  # 1 for first ring, 2 for background cells)
  background_sample$lab <- 2

  ## Get the window
  # if window is specified, use the specified window
  # otherwise, use the window of the background sample
  if (is.null(win)) {
    X <- max(background_sample$Cell.X.Position)
    Y <- max(background_sample$Cell.Y.Position)
    win <- spatstat.geom::owin(c(0, X), c(0,Y))
  }

  ## Default phenotype is specified by bg_type
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

    # if the location of the cluster is not specified,
    # generate a location as the centre of the cluster
    if (is.null(centre_loc)){
      seed_point <- spatstat.core::runifpoint(1, win=win)}
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
        # assign the primary label of the cell
        background_sample[i, "lab"] <- 0
        # generate random number to decide the phenotype
        random <- stats::runif(1)
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
        # determine the primary label of the cell, if the primary label is lower
        # than 2, keep the primary label, skip out of the conditional
        background_sample[i, "lab"] <- min(1, background_sample[i, "lab"])
        if (background_sample[i , "lab"] == 1){
          # generate random number to decide the phenotype
          random <- stats::runif(1)
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
      }
      background_sample[i, "Phenotype"] <- pheno
    }
  }

  return(background_sample)
}

