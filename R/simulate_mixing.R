#' Simulate mixed background image
#' @description Based on an existing background image, simulate mixed cell types
#'   with specified cell types and proportions.
#' @param background_sample (OPTIONAL) Data.Frame. The image that the stripes
#'   are simulated on. By default use the `bg1` background image from the
#'   package.
#' @param names_of_mixing String Vector of the mixed cell types.
#' @param mixing_degree Numeric Vector of the proportions of the mixed cell
#'   types.
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
#' mix_background <- simulate_mixing(background_sample = bg1,
#' names_of_mixing = c("Tumour","Immune", "Others"), mixing_degree = c(0.2, 0.4,  0.4))
#'
#' # library(SPIAT)
#' # plot_cell_categories(mix_background, categories_of_interest = c("Tumour","Immune"),
#' # colour_vector = c("red", "darkgreen"), feature_colname = "Phenotype")

simulate_mixing <- function(background_sample = bg1,
                            names_of_mixing = c("Tumour", "Immune", "Others"),
                            mixing_degree = c(0.2, 0.4, 0.4)) {

  # CHECK is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SummarizedExperiment::colData(background_sample))}

  # default phenotype is "Others"
  if (is.null(background_sample$Phenotype)){
    background_sample[, "Phenotype"] <- "Others"
  }

  n_types <- length(names_of_mixing)
  for (i in seq_len(dim(background_sample)[1])){
    x <- background_sample[i, "Cell.X.Position"]
    y <- background_sample[i, "Cell.Y.Position"]

    random <- stats::runif(1)

    # if the random number falls in the range of an infiltration proportion,
    # pheno will be the corresponding infiltraiton type
    n <- 1 # start from the first proportion
    current_p <- 0
    while (n <= n_types){
      current_p <- current_p + mixing_degree[n]
      if (random <= current_p) {
        pheno <- names_of_mixing[n]
        break
      }
      n <- n+1
    }
    background_sample[i, "Phenotype"] <- pheno
  }

  return(background_sample)
}
