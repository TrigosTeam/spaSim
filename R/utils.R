# convert spe object to a data frame with only colData
get_colData <- function(spe_object){
    formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
    formatted_data <- cbind(formatted_data,
                            data.frame(spatialCoords(spe_object)))
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")

    # delete column `sample_id`
    formatted_data$sample_id <- NULL

    return(formatted_data)
}


#' format_spe
#' @description Format a data frame object into a `SpatialExperiment` class
#'   object where the count assay is empty.
#' @param df A data frame where each row contains information about a cell.
#'   The columns of the data frame will become the colData of the
#'   `SpatialExperiment` object.
#' @import SpatialExperiment
#' @return An SPE object

format_spe <- function(df) {
    #CHECK
    if (dim(df)[1]==0){
        stop("No data in the data frame!")
    }

    assay_data <- rep(0, dim(df)[1])
    assay_rownames <- "pseudo"
    assay_colnames <- rownames(df)

    #transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)

    spe <- SpatialExperiment::SpatialExperiment(
        assays = assay_data_matrix_t, colData = df,
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))

    rownames(spe) <- assay_rownames
    colnames(spe) <- assay_colnames

    return(spe)
}

