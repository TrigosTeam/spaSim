#' TIS
#'
#' @description Tissue Image Simulator (TIS) integrates the basic simulation functions
#' in spaSim, including simulating (mixed) background image, clusters, immune rings,
#' double immune rings and stripes. The patterns are simulated on separate layers
#' sequentially (e.g. immune rings are simulated after/on top of background cells).
#'
#' Pattern properties (e.g. `properties_of_clusters`) contain the properties of
#' a pattern in the format of list where each element is one pattern.
#' These properties can be manually defined. Or users can use predefined properties
#' from the package (`C_shape1`, `C_shape2`, `C_Shape3` for clusters; `R_shape1`,
#' `R_shape2`, `R_shape3` for immune rings; `D_shape1` for double rings).
#' Details about the format of the properties see \link{simulate_clusters}
#'
#' @param background_sample (OPTIONAL) Data.frame or SingleCellExperiment object
#' with locations of points representing background cells. If NULL, background cells
#' will be simulated in this function.
#' @param n_cells (OPTIONAL) Number of background cells to simulate. If NULL,
#' already have a background image.
#' @param width (OPTIONAL) Number The width of the image.
#' @param height (OPTIONAL) Number The height of the image.
#' @param min_d (OPTIONAL) Number The minimum distance between two cells.
#' @param names_of_bg_cells (OPTIONAL) Vector The cell types of the background
#' cells. If NULL, the background cells are of one type.
#' @param proportions_of_bg_cells (OPTIONAL) Vector The corresponding proportion
#' of each cell type in the background cells.
#' @param n_clusters (OPTIONAL) Number of cell clusters. If NULL, no clusters to
#' simulate.
#' @param properties_of_clusters (OPTIONAL) List of parameters to define the
#' clusters.
#' @param n_immune_rings (OPTIONAL) Number of immune rings. If NULL, no immune
#' rings to simulate.
#' @param properties_of_immune_rings (OPTIONAL) List of parameters to define the
#' immune rings.
#' @param n_double_rings (OPTIONAL) Number of double immune rings. If NULL, no double
#' rings to simulate.
#' @param properties_of_double_rings (OPTIONAL) List of parameters to define the
#' double immune rings.
#' @param n_stripe_type (OPTIONAL) Number of stripe (vessel) types. If NULL, no
#' stripes to simulate.
#' @param properties_of_stripes (OPTIONAL) List of parameters to define the stripes.
#' @param image_name (OPTIONAL) String to name the output tissue image.
#'
#' @return An sce object of the simulated image
#' @export
#' @examples
#' set.seed(610)
#' double_ring_image <- TIS(background_sample=bg1, n_clusters = 1,
#' properties_of_clusters = list(C1 = list( name_of_cluster_cell = "Tumour",
#' size = 300, shape = "Oval", centre_loc = data.frame("x" = 500, "y" = 500),
#' infiltration_types = c("Immune1", "Others"), infiltration_proportions = c(0.1, 0.05))))
#'
#' # library(SPIAT)
#' # plot_cell_categories(double_ring_image, categories_of_interest =
#' # c("Tumour","Immune1", "Immune2"), colour_vector = c("red", "darkgreen",
#' # "darkblue"), feature_colname = "Phenotype")

TIS <- function(background_sample = NULL,
                n_cells = NULL,
                width = NULL,
                height = NULL,
                min_d = NULL,

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

                image_name = NULL)
{
  if (is.null(background_sample)){
    background_sample <- simulate_background_cells(n_cells, width, height, min_d)
    X <- width
    Y <- height
  }

  # CHECK is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SingleCellExperiment::colData(background_sample))}

  image <- background_sample
  X <- max(background_sample$Cell.X.Position)
  Y <- max(background_sample$Cell.Y.Position)

  # get background information
  bg_size = c(X, Y)

  # simulate bg with mixing types of cells
  if (!is.null(names_of_bg_cells)){
    image <- simulate_mixing(background_sample = image,
                             names_of_mixing = names_of_bg_cells,
                             mixing_degree = proportions_of_bg_cells)
  }
  # simulate clusters
  if (!is.null(n_clusters)){
    image <- simulate_clusters(background_sample = image,
                               n_clusters = n_clusters,
                               bg_type = "Others",
                               win = NULL,
                               properties_of_clusters = properties_of_clusters)
  }
  # simulate_immune_rings
  if (!is.null(n_immune_rings)){
    image <- simulate_immune_rings(background_sample = image,
                                   n_immune_rings = n_immune_rings,
                                   bg_type = "Others",
                                   win = NULL,
                                   properties_of_immune_rings = properties_of_immune_rings)
  }
  # simulate_double_rings
  if (!is.null(n_double_rings)){
    image <- simulate_double_rings(background_sample = image,
                                   n_double_rings = n_double_rings,
                                   bg_type = "Others",
                                   win = NULL,
                                   properties_of_double_rings = properties_of_double_rings)
  }
  # simulate_stripes
  if (!is.null(n_stripe_type)){
    image <- simulate_stripes(background_sample = image,
                              n_stripe_type = n_stripe_type,
                              win = NULL,
                              properties_of_stripes = properties_of_stripes)}

  # format sce object
  sce <- format_sce(image)
  attr(sce, "name") <- image_name
  return(sce)
}
