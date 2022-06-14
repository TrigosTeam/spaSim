#' Simulate clusters
#' @description Based on an existing background image, simulate clusters of
#'   cells where the same type of cells aggregate. The default values for the
#'   arguments give an example of cluster simulation which enable an automatic
#'   simulation of clusters without the specification of any argument.
#'
#' @param bg_sample (OPTIONAL) A data.frame or SingleCellExperiment class object
#'   with locations of points representing background cells. Further cell types
#'   will be simulated based on this background sample. The data.frame or the
#'   metadata of the SCE object should have colnames including
#'   "Cell.X.Positions" and "Cell.Y.Positions". By default use the internal
#'   \code{\link{bg1}} background image.
#' @param n_clusters Numeric. Number of clusters. This must match the
#'   `length(cluster_properties)`.
#' @param bg_type (OPTIONAL) String. The name of the background cell type if the
#'   background sample does not have a "Phenotype" column. By default is
#'   "Others".
#' @param cluster_properties List of properties of the clusters. See examples
#'   for the format of this arg.
#' @param plot_image Boolean. Whether the simulated image is plotted.
#' @param plot_categories String Vector specifying the order of the cell
#'   categories to be plotted. Default is NULL - the cell categories under the
#'   "Phenotype" column would be used for plotting.
#' @param plot_colours String Vector specifying the order of the colours that
#'   correspond to the `plot_categories` arg. Default is NULL - the predefined
#'   colour vector would be used for plotting.
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
#' cluster_image <- simulate_clusters(bg_sample = bg1,
#' n_clusters=2, cluster_properties=list(C1=list(name_of_cluster_cell="Tumour",
#' size=300, shape="Oval", centre_loc=data.frame("x"=500, "y"=500),
#' infiltration_types=c("Immune1", "Others"), infiltration_proportions=c(0.1, 0.05)),
#' C2=list(name_of_cluster_cell="Immune1", size=500, shape="Irregular",
#' centre_loc=data.frame("x"=1500,"y"=500), infiltration_types=c("Immune2", "Others"),
#' infiltration_proportions=c(0.1, 0.05))))
simulate_clusters <- function(bg_sample = bg1,
                              n_clusters = 2,
                              bg_type = "Others",
                              cluster_properties = list(
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
                              ),
                              plot_image = TRUE,
                              plot_categories = NULL,
                              plot_colours = NULL
){
    ## CHECK
    if (!is.data.frame(bg_sample) & !methods::is(bg_sample, "SingleCellExperiment")) {
        stop("`bg_sample` should be either a data.frame or a SingleCellExperiment object!")
    }
    if (!is.list(cluster_properties)){
        stop("`cluster_properties` should be a list of lists where each list contains the properties of a cluster!")
    }
    if (length(cluster_properties) != n_clusters){
        stop("`n_clusters` should be the same as the length of `cluster_properties`!")
    }
    for (i in seq_len(length(cluster_properties))){
        if (!setequal(names(cluster_properties[[i]]),
                      c("name_of_cluster_cell", "size", "shape", "centre_loc",
                        "infiltration_types", "infiltration_proportions"))) {
            stop("`cluster_properties` is a list of lists. Each list under `cluster_properties` should contain fields:
`name_of_cluster_cell`, `size`, `shape`, `centre_loc`, `infiltration_types`, `infiltration_proportions`.")
        }
        if (length(cluster_properties[[i]]$infiltration_types) != length(cluster_properties[[i]]$infiltration_proportions)){
            stop(paste("The", i, "th list of `cluster_properties` has different length of `infiltration_types` and `infiltration_proportions`.", sep = ""))
        }
    }

    if (!is.null(plot_colours) & !is.null(plot_categories)){
        if (length(plot_categories) != length(plot_colours)){
            stop("`plot_categories` and `plot_colours` should be of the same length!")}}

    if (methods::is(bg_sample, "SingleCellExperiment")) {
        bg_sample <- data.frame(SummarizedExperiment::colData(bg_sample))}

    # Get the window, use the window of the background sample
    X <- max(bg_sample$Cell.X.Position)
    Y <- max(bg_sample$Cell.Y.Position)
    win <- spatstat.geom::owin(c(0, X), c(0,Y))

    # Default phenotype is specified by bg_type
    # (when background sample does not have Phenotype)
    if (is.null(bg_sample$Phenotype)){
        bg_sample[, "Phenotype"] <- bg_type
    }

    n_cells <- dim(bg_sample)[1]

    for (k in seq_len(n_clusters)) { # for each cluster
        # get the arguments

        cell_type <- cluster_properties[[k]]$name_of_cluster_cell
        size <- cluster_properties[[k]]$size
        shape <- cluster_properties[[k]]$shape
        centre_loc <- cluster_properties[[k]]$centre_loc
        infiltration_types <- cluster_properties[[k]]$infiltration_types
        infiltration_proportions <- cluster_properties[[k]]$infiltration_proportions

        # generate a location as the centre of the cluster
        if (is.null(centre_loc)){
            seed_point <- spatstat.random::runifpoint(1, win=win)}
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

        for (i in seq_len(n_cells)){
            x <- bg_sample[i, "Cell.X.Position"]
            y <- bg_sample[i, "Cell.Y.Position"]
            pheno <- bg_sample[i, "Phenotype"]

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
                        if (random <= current_p) {
                            pheno <- infiltration_types[n]
                            break
                        }
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
                        if (random <= current_p) {
                            pheno <- infiltration_types[n]
                            break
                        }
                        n <- n+1 }
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
