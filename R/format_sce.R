#' format_sce
#'
#' @description Format a data frame object into a SingleCellExperiment class
#'   object where the count assay is empty.
#'
#' @param data A data frame where each row contains information about a cell.
#'   The columns of the data frame will become the colData of the
#'   SingleCellExperiment object.
#'
#' @import SingleCellExperiment
#' @return An SingleCellExperiment object

format_sce <- function(data) {

  #CHECK
  if (dim(data)[1]==0){
    print(1)
    stop("No data in the dataframe")
  }

  data[,"pseudo"] <- 0
  assay_data <- data[,"pseudo"]
  data[,"pseudo"] <- NULL
  assay_rownames <- "pseudo"
  assay_colnames <- rownames(data)

  #transpose the matrix so every column is a cell and every row is a marker
  assay_data_matrix <- as.matrix(assay_data)
  colnames(assay_data_matrix) <- NULL
  rownames(assay_data_matrix) <- NULL
  assay_data_matrix_t <- t(assay_data_matrix)

  sce <- SingleCellExperiment(assays = list(counts = assay_data_matrix_t))

  rownames(sce) <- assay_rownames
  colnames(sce) <- assay_colnames

  #Assign the columns
  for (name in colnames(data)){
    sce[[name]] <- data[,name]
  }
  return(sce)
}
