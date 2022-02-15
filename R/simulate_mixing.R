#' simulate_mixing
#'
#' @param background_sample Data.Frame The image that the stripes are simulated on
#' @param names_of_mixing Vector of the mixed cell types
#' @param mixing_degree Vector of the proportions of the mixed cell types
#'
#' @return A data.frame of the simulated image
#' @export
#'
#' @examples
#' set.seed(610)
#' mix_background <- simulate_mixing(background_sample = bg1,
#' names_of_mixing = c("Tumour","Immune", "Others"), mixing_degree = c(0.2, 0.4,  0.4))
#'
#' # NOT RUN
#' # library(SPIAT)
#' # plot_cell_categories(mix_background, categories_of_interest = c("Tumour","Immune"),
#' # colour_vector = c("red", "darkgreen"), feature_colname = "Phenotype")

simulate_mixing <- function(background_sample,
                            names_of_mixing = c("Tumour", "Immune", "Others"),
                            mixing_degree = c(0.2, 0.4, 0.4)) {

  # CHECK is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SingleCellExperiment::colData(background_sample))}

  # default phenotype is "Others"
  if (is.null(background_sample$Phenotype)){
    background_sample[, "Phenotype"] <- "Others"
  }

  n_types <- length(names_of_mixing)
  for (i in 1:dim(background_sample)[1]){
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
