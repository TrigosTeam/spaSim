#' simulate_immune_rings
#'
#' @description Based on an existing background image, simulate rings of immune
#'   cells that surround tumour clusters. The tumour clusters and immune rings
#'   are simulated at the same time.
#'
#' @param bg_sample (OPTIONAL) A data.frame or SingleCellExperiment class object
#'   with locations of points representing background cells. Further cell types
#'   will be simulated based on this background sample. The data.frame or the
#'   metadata of the SCE object should have colnames including
#'   "Cell.X.Positions" and "Cell.Y.Positions". By default use the internal
#'   \code{\link{bg1}} background image.
#' @param bg_type (OPTIONAL) String The name of the background cell type. By
#'   default is "Others".
#' @param n_ir Number of immune rings. This must match the arg
#'   `length(ir_properties)`.
#' @param win (OPTIONAL) owin object output from spatstat.geom::owin function.
#'   By default is the window of the background image.
#' @param ir_properties List of properties of the immune rings. Please refer to the examples for the structure of `ir_properties`.
#' @param plot_image Boolean. Whether the simulated image is plotted.
#' @param plot_categories String Vector specifying the order of the cell
#'   categories to be plotted.
#' @param plot_colours String Vector specifying the order of the colours that
#'   correspond to the `plot_categories` arg.
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
#' ir_properties <- list(I1 = list(name_of_cluster_cell = "Tumour",
#' size = 600,shape = "Circle",centre_loc = data.frame("x" = 930, "y" = 1000),
#' infiltration_types = c("Immune1", "Immune2", "Others"), infiltration_proportions
#' = c(0.15, 0.05, 0.05), name_of_ring_cell = "Immune1", immune_ring_width = 150,
#' immune_ring_infiltration_types = c("Others"), immune_ring_infiltration_proportions = c(0.15)))
#' # simulate immune rings (`n_ir` should match the length of `ir_properties`)
#' immune_ring_image <- simulate_immune_rings(bg_sample = bg1,
#' n_ir = 1, ir_properties = ir_properties)
#' 
simulate_immune_rings <- function(bg_sample = bg1,
                                  bg_type = "Others",
                                  n_ir = 2,
                                  win = NULL,
                                  ir_properties = list(
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
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = NULL,
                                  plot_colours = NULL
) {

  ## CHECK
  # is the background sample a dataframe?
  if (!is.data.frame(bg_sample)) {
    bg_sample <- data.frame(SummarizedExperiment::colData(bg_sample))}

  # check if the specified cluster properties match n_ir
  if (as.numeric(length(ir_properties)) != n_ir){
    stop("`n_ir` does not match the length of `ir_properties`!")
  }

  # add a new column to store the position label for each cell (0 for core cluster,
  # 1 for first ring, 2 for background cells)
  bg_sample$lab <- 2

  ## Get the window
  # if window is specified, use the specified window
  # otherwise, use the window of the background sample
  if (is.null(win)) {
    X <- max(bg_sample$Cell.X.Position)
    Y <- max(bg_sample$Cell.Y.Position)
    win <- spatstat.geom::owin(c(0, X), c(0,Y))
  }

  ## Default phenotype is specified by bg_type
  # (when background sample does not have Phenotype)
  if (is.null(bg_sample$Phenotype)){
    bg_sample[, "Phenotype"] <- bg_type
  }

  n_cells <- dim(bg_sample)[1]

  for (k in seq_len(n_ir)) { # for each cluster
    # get the arguments
    cluster_cell_type <- ir_properties[[k]]$name_of_cluster_cell
    size <- ir_properties[[k]]$size
    shape <- ir_properties[[k]]$shape
    centre_loc <- ir_properties[[k]]$centre_loc
    infiltration_types <- ir_properties[[k]]$infiltration_types
    infiltration_proportions <- ir_properties[[k]]$infiltration_proportions
    ring_cell_type = ir_properties[[k]]$name_of_ring_cell
    ring_width = ir_properties[[k]]$immune_ring_width
    ring_infiltration_types = ir_properties[[k]]$immune_ring_infiltration_types
    ring_infiltration_proportions = ir_properties[[k]]$immune_ring_infiltration_proportions

    # if the location of the cluster is not specified,
    # generate a location as the centre of the cluster
    if (is.null(centre_loc)){
      seed_point <- spatstat.random::runifpoint(1, win=win)}
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
    for (i in seq_len(n_cells)){
      x <- bg_sample[i, "Cell.X.Position"]
      y <- bg_sample[i, "Cell.Y.Position"]
      pheno <- bg_sample[i, "Phenotype"]

      # squared distance to the cluster centre
      A <- (x - a)^2
      B <- (y - b)^2
      AB <- (x-a)*(y-b)
      D <- Circle*(A + B) + Oval*(A + AB + B)

      # determine which region the point falls in
      if (D < R){
        # assign the primary label of the cell
        bg_sample[i, "lab"] <- 0
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
        bg_sample[i, "lab"] <- min(1, bg_sample[i, "lab"])
        if (bg_sample[i , "lab"] == 1){
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
      bg_sample[i, "Phenotype"] <- pheno
    }
  }
  if (plot_image){
    if(is.null(plot_categories)) plot_categories <- unique(bg_sample$Phenotype)
    if (is.null(plot_colours)){
      plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown", "purple", "lightblue",
                        "lightgreen", "yellow", "black", "pink")
    }
    phenos <- plot_categories
    plot_cells(bg_sample, phenos, plot_colours[seq_len(length(phenos))], "Phenotype")
  }
  return(bg_sample)
}

