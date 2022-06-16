#' Simulate multiple background images (mixed cell types)
#'
#' @description Generate a set of background images with different proportions
#'   of mixed cell types all at once. The default values for the arguments give
#'   an example of multiple image simulation which enable an automatic multiple
#'   image simulation without the specification of any argument.
#'
#' @param bg_sample A data frame or `SpatialExperiment` class object with
#'   locations of points representing background cells. Further cell types will
#'   be simulated based on this background sample. The data.frame or the
#'   `spatialCoords()` of the SPE object should have colnames including
#'   "Cell.X.Positions" and "Cell.Y.Positions". By default use the internal
#'   \code{\link{bg1}} background image.
#' @param idents String Vector. Names of the cell types to generate.
#' @param props List. Each element is a vector of proportions of the
#'   corresponding cell type. The length of the vector is how many images to
#'   generate. All vectors should be of the same length, also equal to the
#'   number of images.
#' @param plot_image Boolean. Whether plot the simulated images or not. Default
#'   is TRUE.
#' @param plot_colours String Vector. If plot_image is TRUE, this param is the
#'   corresponding colours for the `idents` arg. Default is NULL, the predefined
#'   colour vector would be used for plotting.
#'
#' @family simulate multiple images functions
#' @seealso \code{\link{multiple_images_with_clusters}} for simulating multiple
#'   images with clusters, and \code{\link{multiple_images_with_immune_rings}}
#'   for simulating multiple images with immune rings.
#'
#' @return A list of SPE objects
#' @export
#' @examples
#' idents <- c("Tumour","Immune","Others")
#' prop1 <- rep(0.1,9)
#' prop2 <- seq(0, 0.4, 0.05)
#' prop3 <- seq(0.9,0.5,-0.05)
#' set.seed(610)
#' bg_image_list <- multiple_background_images(bg_sample = bg1,
#' idents = idents, props = list(prop1, prop2, prop3), plot_image = FALSE)

multiple_background_images <- function(bg_sample = bg1,
                                       idents = c("Tumour", "Immune","Others"),
                                       props = list(rep(0.1, 9),
                                                    seq(0, 0.4, 0.05),
                                                    seq(0.9, 0.5, -0.05)),
                                       plot_image = TRUE,
                                       plot_colours = NULL){
    ## CHECK
    if (!is.data.frame(bg_sample) & !methods::is(bg_sample, "SpatialExperiment")) {
        stop("`bg_sample` should be either a data.frame or a SpatialExperiment object!")
    }
    if (length(idents) != length(props)){
        stop("`idents` and `props` should be of the same length!")
    }
    if (!is.character(idents)){
        stop("`idents` should be character or a character vector!")
    }
    if (!is.list(props)){
        stop("`props` should be a list of numeric vectors. See `?multiple_background_images` for the details of this parameter.")
    }
    if (!is.null(plot_colours)){
        if (length(idents) != length(plot_colours)){
            stop("`idents` and `plot_colours` should be of the same length!")}}
    if (methods::is(bg_sample, "SpatialExperiment")) {
        bg_sample <- get_colData(bg_sample)}
    # default cell type is "Others" (only when there is no "Cell.Type" column in bg_sample)
    if (is.null(bg_sample$Cell.Type)) {bg_sample[, "Cell.Type"] <- "Others"}

    # define the plotting properties
    if (plot_image){
        if (is.null(plot_colours)){
            plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown",
                              "purple", "lightblue", "lightgreen", "yellow",
                              "black", "pink")}}
    n_types <- length(idents)
    # count the image number
    p_idx <- 0
    list.images <- list()

    # loop through the proportions of cell types
    for (prop in props[[1]]){
        p_idx <- p_idx + 1 # this is the p_idx(th) image (also the p_idx proportion)
        # get the vector of proportions for the current image
        props_temp <- c(prop)
        for (k in 2:length(idents)){
            props_temp <- c(props_temp, props[[k]][p_idx])}
        # assign cell type to each cell in the current image
        for (i in seq_len(dim(bg_sample)[1])){
            r <- stats::runif(1)
            # if the random number falls in the range of a proportion,
            # pheno will be the corresponding infiltraiton type
            n <- 1 # start from the first proportion
            current_p <- 0
            while (n <= n_types){
                current_p <- current_p + props_temp[n]
                if (r <= current_p) {
                    pheno <- idents[n]
                    break}
                n <- n+1}
            bg_sample[i, "Cell.Type"] <- pheno
        }

        spe <- format_spe(bg_sample)
        if (plot_image){
            plot_cells(bg_sample, idents, plot_colours[seq_len(length(idents))], "Cell.Type")}
        list.images[[p_idx]] <- spe
    }
    return(list.images)
}
