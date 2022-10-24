#' Tissue Image Simulator (TIS)
#'
#' @description Tissue Image Simulator (TIS) integrates the basic simulation
#'   functions in spaSim, including simulating (mixed) background image,
#'   clusters, immune rings, double immune rings and stripes. The patterns are
#'   simulated on separate layers sequentially (e.g. immune rings are simulated
#'   after/on top of background cells). And each layer is also plot
#'   sequentially.
#'
#'   Pattern properties (e.g. `properties_of_clusters`) contain the properties
#'   of a pattern in the format of list where each element is one pattern. These
#'   properties need to be manually defined. Details about the format of the
#'   properties see the examples in \link{simulate_clusters}
#'   \link{simulate_immune_rings} \link{simulate_double_rings}
#'   \link{simulate_stripes}
#'
#' @param bg_sample (OPTIONAL) A data frame or `SpatialExperiment` class object
#'   with locations of points representing background cells. Further cell types
#'   will be simulated based on this background sample. The data.frame or the
#'   `spatialCoords` of the SPE object should have colnames including
#'   "Cell.X.Positions" and "Cell.Y.Positions". By default use the internal
#'   \code{\link{bg1}} background image.
#' @param n_cells (OPTIONAL) Number of background cells to simulate. Only when
#'   `bg_sample` is NULL.
#' @param width (OPTIONAL) Number The width of the image.
#' @param height (OPTIONAL) Number The height of the image.
#' @param min_d (OPTIONAL) Number The minimum distance between two cells.
#' @param oversampling_rate (OPTIONAL) Numeric. The multiplier for oversampling.
#'   Without oversampling, the simulation deletes cells that are within `min_d`
#'   from each other, resulting in a less total number of cells than `n_cells`.
#'   Default is 1.2 (this should be set based on `n_cells` and `min_d`; should
#'   always be larger than 1).
#' @param names_of_bg_cells (OPTIONAL) Vector The cell types of the background
#'   cells. If NULL, the background cells are of one type.
#' @param proportions_of_bg_cells (OPTIONAL) Vector The corresponding proportion
#'   of each cell type in the background cells.
#' @param n_clusters (OPTIONAL) Number of cell clusters. If NULL, no clusters to
#'   simulate.
#' @param properties_of_clusters (OPTIONAL) List of parameters to define the
#'   clusters.
#' @param n_immune_rings (OPTIONAL) Number of immune rings. If NULL, no immune
#'   rings to simulate.
#' @param properties_of_immune_rings (OPTIONAL) List of parameters to define the
#'   immune rings.
#' @param n_double_rings (OPTIONAL) Number of double immune rings. If NULL, no
#'   double rings to simulate.
#' @param properties_of_double_rings (OPTIONAL) List of parameters to define the
#'   double immune rings.
#' @param n_stripe_type (OPTIONAL) Number of stripe (vessel) types. If NULL, no
#'   stripes to simulate.
#' @param properties_of_stripes (OPTIONAL) List of parameters to define the
#'   stripes.
#' @param image_name (OPTIONAL) String to name the output tissue image.
#' @param plot_image Boolean. Whether the simulated image is plotted.
#' @param plot_categories String Vector specifying the order of the cell
#'   categories to be plotted. Default is NULL - the cell categories under the
#'   "Cell.Type" column would be used for plotting.
#' @param plot_colours String Vector specifying the order of the colours that
#'   correspond to the `plot_categories` arg. Default is NULL - the predefined
#'   colour vector would be used for plotting.
#'
#' @return An spe object of the simulated image
#' @export
#' @examples
#' set.seed(610)
#' double_ring_image <- TIS(bg_sample=bg1, n_clusters=1,
#' properties_of_clusters=list(C1=list( name_of_cluster_cell="Tumour",
#' size=300, shape="Oval", centre_loc=data.frame("x"=500, "y"=500),
#' infiltration_types=c("Immune1", "Others"), infiltration_proportions=c(0.1, 0.05))),
#' plot_image=TRUE)

TIS <- function(bg_sample = NULL,
                n_cells = NULL,
                width = NULL,
                height = NULL,
                min_d = NULL,
                oversampling_rate = 1.2,
                names_of_bg_cells = NULL,
                proportions_of_bg_cells = NULL,
                n_clusters = NULL,
                properties_of_clusters = NULL,
                n_immune_rings = NULL,
                properties_of_immune_rings = NULL,
                n_double_rings = NULL,
                properties_of_double_rings = NULL,
                n_stripe_type = NULL,
                properties_of_stripes = NULL,
                image_name = NULL,
                plot_image = FALSE,
                plot_categories = NULL,
                plot_colours = NULL) {
    ## CHECK
    if (is.null(bg_sample)) {
        if(!is.numeric(n_cells) | !is.numeric(width) | !is.numeric(height) |
           !is.numeric(min_d) | !is.numeric(oversampling_rate)){
            stop("One or more of `n_cells`, `width`, `height`, `min_d`, `oversampling_rate` is not numeric!")}
        bg_sample <- simulate_background_cells(n_cells, width, height, min_d,
                                               oversampling_rate = oversampling_rate)
        X <- width
        Y <- height
    }
    if (!is.null(bg_sample)){
        if (!is.data.frame(bg_sample) & !methods::is(bg_sample, "SpatialExperiment")) {
            stop("`bg_sample` should be either a data.frame or a SpatialExperiment object!")
        }
    }
    if (!is.null(plot_colours) & !is.null(plot_categories)){
        if (length(plot_categories) != length(plot_colours)){
            stop("`plot_categories` and `plot_colours` should be of the same length!")}}

    if (is.null(plot_categories)) plot_categories <- unique(bg_sample$Cell.Type)

    if (methods::is(bg_sample,"SpatialExperiment")) {
        bg_sample <- get_colData(bg_sample)}


    image <- bg_sample
    X <- max(bg_sample$Cell.X.Position)
    Y <- max(bg_sample$Cell.Y.Position)

    # get background information
    bg_size <- c(X, Y)

    # simulate bg with mixing types of cells ####
    if (!is.null(names_of_bg_cells)) {
        # check
        if (length(names_of_bg_cells) != length(proportions_of_bg_cells)){
            stop("The length of `names_of_bg_cells` does not match the length of `proportions_of_bg_cells`!")
        }
        image <- simulate_mixing(
            bg_sample = image,
            idents = names_of_bg_cells,
            props = proportions_of_bg_cells,
            plot_image = FALSE
            # plot_image = plot_image,
            # plot_colours = plot_colours
        )
    }
    # simulate clusters ####
    if (!is.null(n_clusters)) {
        #  check
        if (!is.list(properties_of_clusters)){
            stop("`properties_of_clusters` should be a list of lists where each list contains the properties of a cluster!")
        }
        if (length(properties_of_clusters) != n_clusters){
            stop("`n_clusters` should be the same as the length of `properties_of_clusters`!")
        }
        for (i in seq_len(length(properties_of_clusters))){
            if (!setequal(names(properties_of_clusters[[i]]),
                          c("name_of_cluster_cell", "size", "shape", "centre_loc",
                            "infiltration_types", "infiltration_proportions"))) {
                stop("`properties_of_clusters` is a list of lists. Each list under `properties_of_clusters` should contain fields:
`name_of_cluster_cell`, `size`, `shape`, `centre_loc`, `infiltration_types`, `infiltration_proportions`.")}
            if (length(properties_of_clusters[[i]]$infiltration_types) != length(properties_of_clusters[[i]]$infiltration_proportions)){
                stop("The ", i, "th list of `properties_of_clusters` has different length of `infiltration_types` and `infiltration_proportions`.")
            }
        }
        image <- simulate_clusters(
            bg_sample = image,
            n_clusters = n_clusters,
            bg_type = "Others",
            cluster_properties = properties_of_clusters,
            plot_image = FALSE
            # plot_image = plot_image,
            # plot_categories = plot_categories,
            # plot_colours = plot_colours
        )
    }
    # simulate_immune_rings ####
    if (!is.null(n_immune_rings)) {
        #check
        if (!is.list(properties_of_immune_rings)){
            stop("`properties_of_immune_rings` should be a list of lists where each list contains the properties of an immune ring!")
        }
        for (i in seq_len(length(properties_of_immune_rings))){
            if (!setequal(names(properties_of_immune_rings[[i]]),
                          c("name_of_cluster_cell", "size", "shape", "centre_loc",
                            "infiltration_types", "infiltration_proportions",
                            "name_of_ring_cell", "immune_ring_width",
                            "immune_ring_infiltration_types", "immune_ring_infiltration_proportions"))) {
                stop("`properties_of_immune_rings` is a list of lists. Each list under `properties_of_immune_rings` should contain fields:
`name_of_cluster_cell`, `size`, `shape`, `centre_loc`, `infiltration_types`, `infiltration_proportions`,
`name_of_ring_cell`, `immune_ring_width`, `immune_ring_infiltration_types`, `immune_ring_infiltration_proportions`.")
            }
            if (length(properties_of_immune_rings[[i]]$infiltration_types) !=
                length(properties_of_immune_rings[[i]]$infiltration_proportions)){
                stop("The ", i, "th list of `properties_of_immune_rings` has different length of `infiltration_types` and `infiltration_proportions`.")
            }
            if (length(properties_of_immune_rings[[i]]$immune_ring_infiltration_types) !=
                length(properties_of_immune_rings[[i]]$immune_ring_infiltration_proportions)){
                stop("The ", i, "th list of `properties_of_immune_rings` has different length of `immune_ring_infiltration_types` and `immune_ring_infiltration_proportions`.")
            }
        }
        image <- simulate_immune_rings(
            bg_sample = image,
            n_ir = n_immune_rings,
            bg_type = "Others",
            ir_properties = properties_of_immune_rings,
            plot_image = FALSE
            # plot_image = plot_image,
            # plot_categories = plot_categories,
            # plot_colours = plot_colours
        )
    }
    # simulate_double_rings ####
    if (!is.null(n_double_rings)) {
        #check
        if (!is.list(properties_of_double_rings)){
            stop("`properties_of_double_rings` should be a list of lists where each list contains the properties of a double ring!")
        }
        for (i in seq_len(length(properties_of_double_rings))){
            if (!setequal(names(properties_of_double_rings[[i]]),
                          c("name_of_cluster_cell", "size", "shape", "centre_loc",
                            "infiltration_types", "infiltration_proportions",
                            "name_of_ring_cell", "immune_ring_width",
                            "immune_ring_infiltration_types", "immune_ring_infiltration_proportions",
                            "name_of_double_ring_cell", "double_ring_width",
                            "double_ring_infiltration_types", "double_ring_infiltration_proportions"))) {
                stop("`properties_of_double_rings` is a list of lists. Each list under `properties_of_double_rings` should contain fields:
`name_of_cluster_cell`, `size`, `shape`, `centre_loc`, `infiltration_types`, `infiltration_proportions`,
`name_of_ring_cell`, `immune_ring_width`, `immune_ring_infiltration_types`, `immune_ring_infiltration_proportions`,
`name_of_double_ring_cell`, `double_ring_width`, `double_ring_infiltration_types`, `double_ring_infiltration_proportions`.")
            }
            if (length(properties_of_double_rings[[i]]$infiltration_types) !=
                length(properties_of_double_rings[[i]]$infiltration_proportions)){
                stop("The ", i, "th list of `properties_of_double_rings` has different length of `infiltration_types` and `infiltration_proportions`.")
            }
            if (length(properties_of_double_rings[[i]]$immune_ring_infiltration_types) !=
                length(properties_of_double_rings[[i]]$immune_ring_infiltration_proportions)){
                stop("The ", i, "th list of `properties_of_double_rings` has different length of `immune_ring_infiltration_types` and `immune_ring_infiltration_proportions`.")
            }
            if (length(properties_of_double_rings[[i]]$double_ring_infiltration_types) !=
                length(properties_of_double_rings[[i]]$double_ring_infiltration_proportions)){
                stop("The ", i, "th list of `properties_of_double_rings` has different length of `double_ring_infiltration_types` and `double_ring_infiltration_proportions`.")
            }
        }
        image <- simulate_double_rings(
            bg_sample = image,
            n_dr = n_double_rings,
            bg_type = "Others",
            dr_properties = properties_of_double_rings,
            plot_image = FALSE
            # plot_image = plot_image,
            # plot_categories = plot_categories,
            # plot_colours = plot_colours
        )
    }
    # simulate_stripes ####
    if (!is.null(n_stripe_type)) {
        #check
        if (!is.list(properties_of_stripes)){
            stop("`properties_of_stripes` should be a list of lists where each list contains the properties of a stripe type!")
        }
        if (length(properties_of_stripes) != n_stripe_type){
            stop("`n_stripe_type` should be the same as the length of `properties_of_stripes`!")
        }
        for (i in seq_len(length(properties_of_stripes))){
            if (!setequal(names(properties_of_stripes[[i]]),
                          c("number_of_stripes", "name_of_stripe_cell", "width_of_stripe",
                            "infiltration_types", "infiltration_proportions"))) {
                stop("`properties_of_stripes` is be a list of lists. Each list under `properties_of_stripes` should contain fields:
`number_of_stripes`, `name_of_stripe_cell`, `width_of_stripe`, `infiltration_types`, `infiltration_proportions`.")
            }
            if (length(properties_of_stripes[[i]]$infiltration_types) !=
                length(properties_of_stripes[[i]]$infiltration_proportions)){
                stop("The ", i, "th list of `properties_of_stripes` has different length of `infiltration_types` and `infiltration_proportions`.")
            }
        }
        image <- simulate_stripes(
            bg_sample = image,
            n_stripe_type = n_stripe_type,
            stripe_properties = properties_of_stripes,
            plot_image = FALSE
            # plot_image = plot_image,
            # plot_categories = plot_categories,
            # plot_colours = plot_colours
        )
    }
    #####
    # format spe object
    spe <- format_spe(image)
    attr(spe, "name") <- image_name

    if(plot_image){
        if(is.null(plot_categories)) plot_categories <- unique(spe$Cell.Type)
        if (is.null(plot_colours)){
            plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown", "purple", "lightblue",
                              "lightgreen", "yellow", "black", "pink")
        }
        phenos <- plot_categories
        plot_cells(image, phenos, plot_colours[seq_len(length(phenos))], "Cell.Type")
    }

    return(spe)
}

