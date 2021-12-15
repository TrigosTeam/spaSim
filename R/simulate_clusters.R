#' simulate_clusters
#'
#' @param background_sample Data.Frame The image that the stripes are simulated on
#' @param n_clusters Number of clusters
#' @param bg_type (OPTIONAL) String The name of the background cell type. By
#' default is "Others".
#' @param win (OPTIONAL) owin object output from spatstat::owin function. By
#' default is the window of the background image
#' @param properties_of_clusters List of properties of the clusters
#'
#' @return A data.frame of the simulated image
#' @export

simulate_clusters <- function(background_sample,
                              n_clusters = 2,
                              bg_type = "Others",
                              win = NULL,
                              properties_of_clusters = list(
                                C1 = list(
                                  name_of_cluster_cell = "Tumour",
                                  size = 900,
                                  shape = "Oval",
                                  centre_loc = data.frame("x" = 500, "y" = 500),
                                  infiltration_types = c("Immune1", "Others"),
                                  infiltration_proportions = c(0.1, 0.05)),
                                C2 = list(
                                  name_of_cluster_cell = "Immune1",
                                  size = 900,
                                  shape = "Irregular",
                                  centre_loc = data.frame("x" = 1500, "y" = 500),
                                  infiltration_types = c("Immune2", "Others"),
                                  infiltration_proportions = c(0.1, 0.05))
                              )
)
{

  # Get the window
  # if window is specified, use the specified window
  # otherwise, use the window of the background sample
  if (is.null(win)) {
    X <- max(background_sample$Cell.X.Position)
    Y <- max(background_sample$Cell.Y.Position)
    win <- spatstat::owin(c(0, X), c(0,Y))
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
      seed_point <- spatstat::runifpoint(1, win=win)}
    else seed_point <- centre_loc
    a <- seed_point$x
    b <- seed_point$y

    r <- size
    R <- r^2
    shape <- shape
    r_theta <- stats::runif(1, min = -2 , max = 1) # for the irregular shape

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
