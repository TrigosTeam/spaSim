#' simulate_stripes
#'
#' @description Based on an existing background image, simulate stripes of cells
#'   representing vessels. The cell types and widths of the stripes can be
#'   specified. The locations of the stripes are randomly simulated. Please
#'   refer to the examples to check what properties of the stripes can be
#'   specified. The default values for the arguments give an example of vessel
#'   simulation which enable an automatic simulation of vessels without the
#'   specification of any argument.
#'
#' @param bg_sample (OPTIONAL) A data frame or `SpatialExperiment` class object
#'   with locations of points representing background cells. Further cell types
#'   will be simulated based on this background sample. The data.frame or the
#'   `spatialCoords()` of the SPE object should have colnames including
#'   "Cell.X.Positions" and "Cell.Y.Positions". By default use the internal
#'   \code{\link{bg1}} background image.
#' @param n_stripe_type Number of stripe types. Should be the same as
#'   `length(stripe_properties`.
#' @param stripe_properties List of the properties of the stripes. See examples
#'   for the format of the properties. Please refer to the examples for the
#'   structure of `stripe_properties`.
#' @param plot_image Boolean. Whether the simulated image is plotted. Default is
#'   TRUE.
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
#'   \code{\link{simulate_clusters}} for cluster simulation, and
#'   \code{\link{simulate_immune_rings}}/\code{\link{simulate_double_rings}} for
#'   immune ring simulation
#'
#' @return A data.frame of the simulated image
#' @export
#' @examples
#' stripe_properties = list(
#' S1 = list(
#'   number_of_stripes = 1,
#'   name_of_stripe_cell = "Others",
#'   width_of_stripe = 80,
#'   infiltration_types = c("Immune"),
#'   infiltration_proportions = c(0.08)
#' ), S2 = list(
#'   number_of_stripes = 1,
#'   name_of_stripe_cell = "Others",
#'   width_of_stripe = 80,
#'   infiltration_types = c("Immune"),
#'   infiltration_proportions = c(0.08)))
#' set.seed(610)
#' stripe_image <- simulate_stripes(bg_sample = bg1, n_stripe_type=2,
#' stripe_properties = stripe_properties, plot_image = TRUE)

simulate_stripes <- function(bg_sample = bg1,
                             n_stripe_type = 2,
                             stripe_properties = list(
                                 S1 = list(
                                     number_of_stripes = 1,
                                     name_of_stripe_cell = "Others",
                                     width_of_stripe = 80,
                                     infiltration_types = c("Immune"),
                                     infiltration_proportions = c(0.08)
                                 ),
                                 S2 = list(
                                     number_of_stripes = 1,
                                     name_of_stripe_cell = "Others",
                                     width_of_stripe = 80,
                                     infiltration_types = c("Immune"),
                                     infiltration_proportions = c(0.08)
                                 )
                             ),
                             plot_image = TRUE,
                             plot_categories = NULL,
                             plot_colours = NULL
){
    ## CHECK
    if (!is.data.frame(bg_sample) & !methods::is(bg_sample, "SpatialExperiment")) {
        stop("`bg_sample` should be either a data.frame or a SpatialExperiment object!")
    }
    if (!is.list(stripe_properties)){
        stop("`stripe_properties` should be a list of lists where each list contains the properties of a stripe type!")
    }
    if (length(stripe_properties) != n_stripe_type){
        stop("`n_stripe_type` should be the same as the length of `stripe_properties`!")
    }
    for (i in seq_len(length(stripe_properties))){
        if (!setequal(names(stripe_properties[[i]]),
                    c("number_of_stripes", "name_of_stripe_cell", "width_of_stripe",
                      "infiltration_types", "infiltration_proportions"))) {
            stop("`stripe_properties` is be a list of lists. Each list under `stripe_properties` should contain fields:
`number_of_stripes`, `name_of_stripe_cell`, `width_of_stripe`, `infiltration_types`, `infiltration_proportions`.")
        }
        if (length(stripe_properties[[i]]$infiltration_types) !=
            length(stripe_properties[[i]]$infiltration_proportions)){
            stop(paste("The", i, "th list of `stripe_properties` has different length of `infiltration_types` and `infiltration_proportions`.", sep = ""))
        }
    }

    if (!is.null(plot_colours) & !is.null(plot_categories)){
        if (length(plot_categories) != length(plot_colours)){
            stop("`plot_categories` and `plot_colours` should be of the same length!")}}

    if (methods::is(bg_sample, "SpatialExperiment")) {
        bg_sample <- get_colData(bg_sample)}
    # get the window
    X <- max(bg_sample$Cell.X.Position)
    Y <- max(bg_sample$Cell.Y.Position)
    win <- spatstat.geom::owin(c(0, X), c(0,Y))
    # default phenotype is "Others"
    if (is.null(bg_sample$Phenotype)){
        bg_sample[, "Phenotype"] <- "Others"
    }
    for (k in seq_len(n_stripe_type)){
        n_stripes <- stripe_properties[[k]]$number_of_stripes
        stripe_cell_type <- stripe_properties[[k]]$name_of_stripe_cell
        stripe_width <- stripe_properties[[k]]$width_of_stripe
        infiltration_types <- stripe_properties[[k]]$infiltration_types
        infiltration_proportions <-
            stripe_properties[[k]]$infiltration_proportions

        # generate intercepts
        random_nums <- stats::runif(n_stripes, min = -max(X,Y), max = max(X,Y))
        for (i in seq_len(dim(bg_sample)[1])){
            x <- bg_sample[i, "Cell.X.Position"]
            y <- bg_sample[i, "Cell.Y.Position"]
            pheno <- bg_sample[i, "Phenotype"]

            p <- tryCatch(which(random_nums ==
                                    max(random_nums[which(random_nums<y-x)])),
                          error=function(e) e, warning=function(w) w)
            if (methods::is(p,"warning") == FALSE) {
                b <- random_nums[p]
                if ( y < x + b + stripe_width ){
                    random <- stats::runif(1)
                    n_infiltration_types <- length(infiltration_types)
                    pheno <- stripe_cell_type

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
                    bg_sample[i, "Phenotype"] <- pheno}}}}
    if (plot_image){
        if(is.null(plot_categories)) plot_categories <-
                unique(bg_sample$Phenotype)
        if (is.null(plot_colours)){
            plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown",
                              "purple", "lightblue", "lightgreen", "yellow",
                              "black", "pink")}
        phenos <- plot_categories
        plot_cells(bg_sample, phenos, plot_colours[seq_len(length(phenos))],
                   "Phenotype")}
    return(bg_sample)
}
