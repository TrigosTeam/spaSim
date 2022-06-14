#' Simulate mixed background image
#' @description Based on an existing background image, simulate mixed cell types
#'   with specified cell types and proportions. The default values for the
#'   arguments give an example of mixed cell type simulation which enable an
#'   automatic simulation of mixed cell types without the specification of any
#'   argument.
#' @param bg_sample (OPTIONAL) A data frame or `SpatialExperiment` class object
#'   with locations of points representing background cells. Further cell types
#'   will be simulated based on this background sample. The data.frame or the
#'   `spatialCoords()` of the SPE object should have colnames including
#'   "Cell.X.Positions" and "Cell.Y.Positions". By default use the internal
#'   \code{\link{bg1}} background image.
#' @param idents String Vector of the mixed cell types.
#' @param props Numeric Vector of the proportions of the mixed cell types.
#' @param plot_image Boolean. Whether the simulated image is plotted.
#' @param plot_colours String Vector specifying the order of the colours that
#'   correspond to the `idents` arg. Default is NULL - the predefined
#'   colour vector would be used for plotting..
#'
#' @family simulate pattern functions
#' @seealso   \code{\link{simulate_background_cells}} for all cell simulation,
#'   \code{\link{simulate_clusters}} for cluster simulation,
#'   \code{\link{simulate_immune_rings}}/\code{\link{simulate_double_rings}} for
#'   immune ring simulation, and \code{\link{simulate_stripes}} for vessel
#'   simulation.
#'
#' @return A data.frame of the simulated image
#' @export
#'
#' @examples
#' set.seed(610)
#' mix_background <- simulate_mixing(bg_sample=bg1,
#' idents=c("Tumour","Immune", "Others"), props=c(0.2, 0.4,  0.4),
#' plot_image=TRUE)
simulate_mixing <- function(bg_sample = bg1,
                            idents = c("Tumour", "Immune", "Others"),
                            props = c(0.2, 0.4, 0.4),
                            plot_image = TRUE,
                            plot_colours = NULL) {

    ## CHECK
    if (length(idents) != length(props)){
        stop("`idents` and `props` should be of the same length!")}
    if (!is.character(idents)){
        stop("`idents` should be character or a character vector!")}
    if (!is.numeric(props)){
        stop("`props` should be numeric or a numeric vector!")}
    if (!is.data.frame(bg_sample) & !methods::is(bg_sample, "SpatialExperiment")) {
        stop("`bg_sample` should be either a data.frame or a SpatialExperiment object!")}
    if (!is.null(plot_colours)){
        if (length(idents) != length(plot_colours)){
            stop("`idents` and `plot_colours` should be of the same length!")}}
    if (methods::is(bg_sample,"SpatialExperiment")) {
        bg_sample <- get_colData(bg_sample)}
    # default cell type is "Others"
    if (is.null(bg_sample$Cell.Type)){
        bg_sample[, "Cell.Type"] <- "Others"
    }

    n_types <- length(idents)
    for (i in seq_len(dim(bg_sample)[1])){
        x <- bg_sample[i, "Cell.X.Position"]
        y <- bg_sample[i, "Cell.Y.Position"]

        random <- stats::runif(1)

        # if the random number falls in the range of an infiltration proportion,
        # pheno will be the corresponding infiltraiton type
        n <- 1 # start from the first proportion
        current_p <- 0
        while (n <= n_types){
            current_p <- current_p + props[n]
            if (random <= current_p) {
                pheno <- idents[n]
                break
            }
            n <- n+1
        }
        bg_sample[i, "Cell.Type"] <- pheno
    }
    if (plot_image){
        if (is.null(plot_colours)){
            plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown", "purple", "lightblue",
                              "lightgreen", "yellow", "black", "pink")
        }
        plot_cells(bg_sample, idents, plot_colours[seq_len(length(idents))], "Cell.Type")
    }

    return(bg_sample)
}
