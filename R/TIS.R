#' Tissue Image Simulator (TIS)
#'
#' @description Tissue Image Simulator (TIS) integrates the basic simulation
#'   functions in spaSim, including simulating (mixed) background image,
#'   clusters, immune rings, double immune rings and stripes. The patterns are
#'   simulated on separate layers sequentially (e.g. immune rings are simulated
#'   after/on top of background cells). And each layer is also plot sequentially.
#'
#'   Pattern properties (e.g. `properties_of_clusters`) contain the properties
#'   of a pattern in the format of list where each element is one pattern. These
#'   properties need to be manually defined. Details about the format of the properties
#'   see the examples in \link{simulate_clusters} \link{simulate_immune_rings}
#'   \link{simulate_double_rings} \link{simulate_stripes}
#'
#' @param bg_sample (OPTIONAL) A data.frame or SingleCellExperiment class object
#'   with locations of points representing background cells. Further cell types
#'   will be simulated based on this background sample. The data.frame or the
#'   metadata of the SCE object should have colnames including
#'   "Cell.X.Positions" and "Cell.Y.Positions". By default use the internal
#'   \code{\link{bg1}} background image.
#' @param n_cells (OPTIONAL) Number of background cells to simulate. Only when
#'   `bg_sample` is NULL.
#' @param width (OPTIONAL) Number The width of the image.
#' @param height (OPTIONAL) Number The height of the image.
#' @param min_d (OPTIONAL) Number The minimum distance between two cells.
#' @param oversampling_rate (OPTIONAL) Numeric. The multiplier for oversampling.
#'   Without oversampling, the simulation deletes cells that are within `min_d`
#'   from each other, resulting in a less number of cells than specified.
#'   Default is 1.2.
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
#'   categories to be plotted.
#' @param plot_colours String Vector specifying the order of the colours that
#'   correspond to the `plot_categories` arg.
#'
#' @return An sce object of the simulated image
#' @export
#' @examples
#' set.seed(610)
#' double_ring_image <- TIS(bg_sample=bg1, n_clusters = 1,
#' properties_of_clusters = list(C1 = list( name_of_cluster_cell = "Tumour",
#' size = 300, shape = "Oval", centre_loc = data.frame("x" = 500, "y" = 500),
#' infiltration_types = c("Immune1", "Others"), infiltration_proportions = c(0.1, 0.05))),
#' plot_image = TRUE)

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
                plot_colours = NULL)
{
  if (is.null(bg_sample)){
    bg_sample <- simulate_background_cells(n_cells, width, height, min_d, 
                                           oversampling_rate = oversampling_rate)
    X <- width
    Y <- height
  }
  if (is.null(plot_colours)){
    plot_colours <- plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown", "purple", "lightblue",
                                      "lightgreen", "yellow", "black", "pink")
  }
  if(is.null(plot_categories)) plot_categories <- unique(image$Phenotype)

  # CHECK is the background sample a dataframe?
  if (!is.data.frame(bg_sample)) {
    bg_sample <- data.frame(SummarizedExperiment::colData(bg_sample))}

  image <- bg_sample
  X <- max(bg_sample$Cell.X.Position)
  Y <- max(bg_sample$Cell.Y.Position)

  # get background information
  bg_size = c(X, Y)

  # simulate bg with mixing types of cells
  if (!is.null(names_of_bg_cells)){
    image <- simulate_mixing(bg_sample = image,
                             idents = names_of_bg_cells,
                             props = proportions_of_bg_cells,
                             plot_image = plot_image,
                             plot_colours = plot_colours)
  }
  # simulate clusters
  if (!is.null(n_clusters)){
    image <- simulate_clusters(bg_sample = image,
                               n_clusters = n_clusters,
                               bg_type = "Others",
                               win = NULL,
                               cluster_properties = properties_of_clusters,
                               plot_image = plot_image,
                               plot_categories = plot_categories,
                               plot_colours = plot_colours)
  }
  # simulate_immune_rings
  if (!is.null(n_immune_rings)){
    image <- simulate_immune_rings(bg_sample = image,
                                   n_ir = n_immune_rings,
                                   bg_type = "Others",
                                   win = NULL,
                                   ir_properties= properties_of_immune_rings,
                                   plot_image = plot_image,
                                   plot_categories = plot_categories,
                                   plot_colours = plot_colours)
  }
  # simulate_double_rings
  if (!is.null(n_double_rings)){
    image <- simulate_double_rings(bg_sample = image,
                                   n_dr = n_double_rings,
                                   bg_type = "Others",
                                   win = NULL,
                                   dr_properties = properties_of_double_rings,
                                   plot_image = plot_image,
                                   plot_categories = plot_categories,
                                   plot_colours = plot_colours)
  }
  # simulate_stripes
  if (!is.null(n_stripe_type)){
    image <- simulate_stripes(bg_sample = image,
                              n_stripe_type = n_stripe_type,
                              win = NULL,
                              stripe_properties = properties_of_stripes,
                              plot_image = plot_image,
                              plot_categories = plot_categories,
                              plot_colours = plot_colours)}

  # format sce object
  sce <- format_sce(image)
  attr(sce, "name") <- image_name

  return(sce)
}
