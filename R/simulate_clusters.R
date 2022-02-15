#' Simulate clusters
#' @description Based on an existing background image, simulate clusters of
#'   cells where the same type of cells aggregate.
#' @param background_sample (OPTIONAL) Data.Frame. The background image that the
#'   clusters are generated on. Default is the internal background image `bg1`.
#' @param n_clusters Numeric. Number of clusters.
#' @param bg_type (OPTIONAL) String. The name of the background cell type if the
#'   background sample does not have a "Phenotype" column. By default is
#'   "Others".
#' @param win (OPTIONAL) `owin` object output from spatstat.geom::owin function.
#'   By default is the window of the background image.
#' @param properties_of_clusters List of properties of the clusters. See
#'   examples for the format of this arg.
#'
#' @family simulate pattern functions
#' @seealso   \code{\link{simulate_background_cells}} for all cell simulation,
#'   \code{\link{simulate_mixing}} for mixed background simulation,
#'   \code{\link{simulate_immune_rings}}/\code{\link{simulate_double_rings}} for
#'   immune ring simulation, and \code{\link{simulate_stripes}} for vessel
#'   simulation.
#'
#' @return A data.frame of the simulated image
#' @export
#' @examples
#' set.seed(610)
#' cluster_image <- simulate_clusters(background_sample = bg1,
#' n_clusters = 1, properties_of_clusters= list(C1 = list(
#' name_of_cluster_cell = "Tumour", size = 300, shape = "Oval", centre_loc =
#' data.frame("x" = 500, "y" = 500), infiltration_types = c("Immune1", "Others"),
#' infiltration_proportions = c(0.1, 0.05))))

simulate_clusters <- function(background_sample = bg1,
                              n_clusters = 2,
                              bg_type = "Others",
                              win = NULL,
                              properties_of_clusters = list(
                                C1 = list(
                                  name_of_cluster_cell = "Tumour",
                                  size = 300,
                                  shape = "Oval",
                                  centre_loc = data.frame("x" = 500, "y" = 500),
                                  infiltration_types = c("Immune1", "Others"),
                                  infiltration_proportions = c(0.1, 0.05)),
                                C2 = list(
                                  name_of_cluster_cell = "Immune1",
                                  size = 500,
                                  shape = "Irregular",
                                  centre_loc = data.frame("x" = 1500, "y" = 500),
                                  infiltration_types = c("Immune2", "Others"),
                                  infiltration_proportions = c(0.1, 0.05))
                              )
){

  # CHECK is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SingleCellExperiment::colData(background_sample))}

  # Get the window
  # if window is specified, use the specified window
  # otherwise, use the window of the background sample
  if (is.null(win)) {
    X <- max(background_sample$Cell.X.Position)
    Y <- max(background_sample$Cell.Y.Position)
    win <- spatstat.geom::owin(c(0, X), c(0,Y))
  }

  # Default phenotype is specified by bg_type
  # (when background sample does not have Phenotype)
  if (is.null(background_sample$Phenotype)){
    background_sample[, "Phenotype"] <- bg_type
  }

  n_cells <- dim(background_sample)[1]

  for (k in 1:n_clusters) { # for each cluster
    # get the arguments

    cell_type <- properties_of_clusters[[k]]$name_of_cluster_cell
    size <- properties_of_clusters[[k]]$size
    shape <- properties_of_clusters[[k]]$shape
    centre_loc <- properties_of_clusters[[k]]$centre_loc
    infiltration_types <- properties_of_clusters[[k]]$infiltration_types
    infiltration_proportions <- properties_of_clusters[[k]]$infiltration_proportions

    # generate a location as the centre of the cluster
    if (is.null(centre_loc)){
      seed_point <- spatstat.core::runifpoint(1, win=win)}
    else seed_point <- centre_loc
    a <- seed_point$x
    b <- seed_point$y

    r <- size
    R <- r^2
    shape <- shape
    #r_theta <- stats::runif(1, min = -2 , max = 1) # for the irregular shape
    # `r_theta` not random
    r_theta <- 0.5

    Circle <- (shape == "Circle")
    Oval <- (shape == "Oval")
    Strip <- (shape == "Strip")

    for (i in 1:n_cells){
      x <- background_sample[i, "Cell.X.Position"]
      y <- background_sample[i, "Cell.Y.Position"]
      pheno <- background_sample[i, "Phenotype"]

      A <- (x - a)^2
      B <- (y - b)^2
      AB <- (x-a)*(y-b)

      if (shape != "Irregular"){
        D <- Circle*(A + B) + Oval*(A + AB + B) + Strip*(A - 1.96*AB + B)

        if (D < R){ # in the region of cluster
          # generate random number to decide the phenotype
          random <- stats::runif(1)

          n_infiltration_types <- length(infiltration_types)

          # default phenotype is cell type of interest of this cluster
          pheno <- cell_type
          # if the random number falls in the range of an infiltration proportion,
          # pheno will be the corresponding infiltraiton type
          n <- 1 # start from the first proportion
          current_p <- 0
          while (n <= n_infiltration_types){
            current_p <- current_p + infiltration_proportions[n]
            if (random <= current_p) pheno <- infiltration_types[n]; break
            n <- n+1
          }

        }
      }

      else { # use heart shape to represent irregular immune cluster
        d <- sqrt(A+B)
        theta <- atan( (y-b)/(x-a) )
        # adjust the calculated angle
        if (AB*(y-b) < 0) theta <- theta + pi  # II or III Quadrant
        else if (AB < 0) theta <- theta + 2*pi  # IV Quadrant

        if (d < 0.55*r+0.45*r*cos(theta) && theta < r_theta + 1 && theta > r_theta){
          random <- stats::runif(1)

          n_infiltration_types <- length(infiltration_types)

          # default phenotype is cell type of interest of this cluster
          pheno <- cell_type
          # if the random number falls in the range of an infiltration proportion,
          # pheno will be the corresponding infiltraiton type
          n <- 1 # start from the first proportion
          current_p <- 0
          while (n <= n_infiltration_types){
            current_p <- current_p + infiltration_proportions[n]
            if (random <= current_p) pheno <- infiltration_types[n]; break
            n <- n+1 }
        }
      }
      background_sample[i, "Phenotype"] <- pheno
    }
  }

  return(background_sample)
}
