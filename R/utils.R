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
