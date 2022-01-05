#' simulate_multiple_background_images
#'
#' @param background_sample Data.frame or SingleCellExperiment object with locations
#' of points representing background cells. Further cell types will be simulated
#' based on this background sample.
#' @param names_of_cell_types Vector of strings Names of the cell types to generate
#' @param proportions_of_cell_types List of vectors. Each vector is a sequence of
#' proportions of the corresponding cell type. The length of the vector is how many
#' images to generate. All vectors should be of the same length, also equal to
#' the number of images
#'
#' @return
#' @export
#'
#' @examples
simulate_multiple_background_images <- function(background_sample,
                                                names_of_cell_types = c("Tumour",
                                                                        "Immune",
                                                                        "Others"),
                                                proportions_of_cell_types = list(
                                                  rep(0.1, 9),
                                                  seq(0, 0.4, 0.05),
                                                  seq(0.9,0.5,-0.05))){
  # CHECK is the background sample a data frame?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SingleCellExperiment::colData(background_sample))}

  # default phenotype is "Others"
  if (is.null(background_sample$Phenotype)){
    background_sample[, "Phenotype"] <- "Others"
  }

  n_types <- length(names_of_cell_types)

  # count the image number
  p_idx <- 0
  list.images <- list()

  # loop through the proportions of cell types
  for (prop in proportions_of_cell_types[[1]]){
    p_idx <- p_idx + 1 # this is the p_idx(th) image (also the p_idx proportion)

    # get the vector of proportions for the current image
    proportions_of_cell_types_temp <- c(prop)
    for (k in 2:length(names_of_cell_types)){
      proportions_of_cell_types_temp <- c(proportions_of_cell_types_temp,
                                          proportions_of_cell_types[[k]][p_idx])
    }

    # assign cell type to each cell in the current image
    print(p_idx)
    for (i in 1:dim(background_sample)[1]){
      r <- stats::runif(1)
      # if the random number falls in the range of a proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_types){
        current_p <- current_p + proportions_of_cell_types_temp[n]
        if (r <= current_p) {
          pheno <- names_of_cell_types[n]
          break
        }
        n <- n+1
      }

      background_sample[i, "Phenotype"] <- pheno
    }
    sce <- format_sce(background_sample)
    plot_cell_categories(sce, categories_of_interest = c("Tumour", "Immune","Others"), colour_vector = c("red" ,"blue", "gray"), feature_colname = "Phenotype")
    list.images[[p_idx]] <- sce
  }
  return(list.images)
}
