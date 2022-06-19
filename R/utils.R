# convert spe object to a data frame with only colData
get_colData <- function(spe_object){
    formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
    formatted_data <- cbind(formatted_data,
                            data.frame(spatialCoords(spe_object)))

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

#' plot_cells
#' @description Produces a scatter plot of the cells in the tissue. Cells are
#'   coloured categorically by `Cell.Type` column. Cell categories not specified
#'   will be coloured "lightgrey" and labled "Unspecified".
#' @param spe_object SpatialExperiment object or a data.frame that has cell
#'   locations and cell type info.
#' @param categories_of_interest String Vector of cell categories to be
#'   coloured.
#' @param colour_vector String Vector specifying the colours of each cell
#'   type.
#' @param feature_colname String specifying the column the cell categories
#'   belong to.
#' @import dplyr
#' @import ggplot2
#' @return A plot is returned

plot_cells <- function(spe_object, categories_of_interest = NULL,
                       colour_vector = NULL, feature_colname = "Cell.Type") {

    # if plotting the structure, users do not have to enter the params
    # we have stored the categories and colours for them
    if (feature_colname == "Structure" & is.null(categories_of_interest)) {
        categories_of_interest <- c("Border",
                                    "Inside",
                                    "Infiltrated.immune",
                                    "Outside",
                                    "Stromal.immune",
                                    "Internal.margin",
                                    "Internal.margin.immune",
                                    "External.margin",
                                    "External.margin.immune")
        colour_vector <- c("black", "pink", "purple", "yellow", "orange", "lightgreen", "darkgreen", "lightblue", "blue")
    }

    # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
    Cell.X.Position <- Cell.Y.Position <- Category <- NULL

    if (methods::is(spe_object,'SpatialExperiment')) {
        formatted_data <- get_colData(spe_object)
    }
    else formatted_data <- spe_object

    #CHECK
    if (length(categories_of_interest) != length(colour_vector)) {
        stop("The colour vector is not the same length as the cell types of interest")
    }

    # if some categories are not in the data, delete them from the categories_of_interest vector
    # delete the corresponding colours as well
    # return a message informing the deleted category
    for (category in categories_of_interest) {
        if (!(category %in% unique(formatted_data[[feature_colname]]))) {
            cat_idx <- match(category, categories_of_interest)
            categories_of_interest <- categories_of_interest[-cat_idx]
            colour_vector <- colour_vector[-cat_idx]
            methods::show(paste(category, "cells were not found and not plotted."))
        }
    }

    #set all categories of those that aren't in categories_of_interest to be "Unspecified"
    if (any(!formatted_data[[feature_colname]] %in% categories_of_interest)) {
        formatted_data[!formatted_data$Cell.Type %in% categories_of_interest,][[feature_colname]] <- "Unspecified"
    }

    #Assign the colour to corresponding cell types in df
    formatted_data$color <- ""
    for (category in categories_of_interest) {
        idx <- which(categories_of_interest == category)
        formatted_data[formatted_data[[feature_colname]] == category, ]$color <- colour_vector[idx]
    }
    if (any(formatted_data[[feature_colname]] == "Unspecified")) {
        formatted_data[formatted_data[[feature_colname]] == "Unspecified", ]$color <- "lightgrey"
        all_categories <- c(categories_of_interest, "Unspecified")
        all_colours <- c(colour_vector, "lightgrey")
    } else {
        all_categories <- categories_of_interest
        all_colours <- colour_vector
    }

    p <- ggplot(formatted_data, aes_string(x = "Cell.X.Position",
                                           y = "Cell.Y.Position", colour = feature_colname)) +
        geom_point(aes_string(colour = feature_colname), size = 1)
    p <- ggplot(formatted_data, aes_string(x = "Cell.X.Position", y = "Cell.Y.Position",
                                           colour = feature_colname))
    if (any(formatted_data[[feature_colname]] == "Unspecified")) {
        p <- p + geom_point(data=subset(formatted_data, get(feature_colname) =='Unspecified'),
                            aes_string(colour = feature_colname), size = 1) +
            geom_point(data=subset(formatted_data, get(feature_colname) !='Unspecified'),
                       aes_string(colour = feature_colname), size = 1)
    }else{
        p <- p + geom_point(aes_string(colour = feature_colname), size = 1)}

    p <- p +
        guides(alpha = "none") +
        labs(colour = feature_colname) +
        scale_color_manual(breaks = all_categories, values=all_colours) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    methods::show(p)
}

