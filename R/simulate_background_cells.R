#' Simulate background cells
#'
#' @description Simulate cell locations. The 2D locations of the cells are
#'   simulated and plotted in a rectangular window. Users can specify the window
#'   size, cell number and the minimum distance between two cells. All cells
#'   have the same cell type, specified by the "Cell.Type" param.
#'
#' @details There are two options for the background cell distribution model. 1)
#'   Hardcore model for tumour tissues. This model uses `rHardcore`
#'   [`spatstat.random`]. Since `rHardcore` deletes cells that are within the
#'   distance of `min_d` to another cell, resulting in fewer cell specified by
#'   users, we here introduce parameter `oversampling_rate` to generate more
#'   cells than specified. 2) Normal tissues use an evenly-spaced model where
#'   the cells are distributed approximately according to the vertices of a
#'   hexagon. The function accomplishes this by generating cells on a hexagonal
#'   grid and individually applying a bounded uniform random jitter. In our
#'   algorithm, `jitter` is the parameter to define the uniform distribution of
#'   the jitter of the cells from the hexagon vertices. If there is a reference
#'   image, `jitter` can be estimated by comparing the average minimum distance
#'   between cells of the simulated image and the reference image. If without a
#'   reference image, We suggest 0.3 as the default value of `jitter` as this
#'   gives a sensible outcome.
#'
#' @param n_cells Numeric. Number of cells to simulate in the background.
#' @param width,height Numeric. The width and height of the image.
#' @param method String. The distribution model for the background cells.
#'   Options are "Hardcore" for tumour tissues and "Even" for normal tissues.
#'   Default is "Hardcore".
#' @param min_d (OPTIONAL) Numeric. Use when `method` is "Hardcore". The minimum
#'   distance between two cells.
#' @param oversampling_rate (OPTIONAL) Numeric. Use when `method` is "Hardcore".
#'   The multiplier for oversampling. Without oversampling, the simulation
#'   deletes cells that are within `min_d` from each other, resulting in a less
#'   total number of cells than `n_cells`. Default is 1.2 (this should be set
#'   based on `n_cells` and `min_d`; should always be larger than 1).
#' @param jitter (OPTIONAL) Numeric. Use when `method` is "Even". The uniform
#'   distribution parameter to generate the jitter distance for each cell from
#'   the vertices of the hexagon.
#' @param Cell.Type (OPTIONAL) String. The name of the background cell type.
#'   Default is "Others" since there shouldn't be any identity of the background
#'   cells.
#' @param plot_image (OPTIONAL) Boolean. Default is TRUE.
#'
#' @family simulate pattern functions
#' @seealso \code{\link{simulate_mixing}} for mixed background simulation,
#'   \code{\link{simulate_clusters}} for cluster simulation,
#'   \code{\link{simulate_immune_rings}}/\code{\link{simulate_double_rings}} for
#'   immune ring simulation, and \code{\link{simulate_stripes}} for vessel
#'   simulation.
#'
#' @return A data.frame of the simulated background image
#' @export
#'
#' @examples
#' set.seed(610) # set seed for this background image simulation for reproducibility
#' background_image <- simulate_background_cells(n_cells = 5000, width = 2000,
#'                                               height = 2000, method = "Hardcore",
#'                                               min_d = 10,
#'                                               oversampling_rate = 1.5,
#'                                               Cell.Type = "Others",
#'                                               plot_image = TRUE)

simulate_background_cells <- function(n_cells, width, height, method = "Hardcore",
                                      min_d, oversampling_rate = 1.2, jitter = 0.3,
                                      Cell.Type = "Others", plot_image = TRUE){
    ## CHECK
    if(!is.numeric(n_cells) | !is.numeric(width) | !is.numeric(height)){
        stop("One or more of `n_cells`, `width`, `height` is not numeric!")
    }
    if (!is.character(Cell.Type)){
        stop("`Cell.Type` should be of character type!")
    }

    if (method == "Even"){
        if(!is.numeric(jitter)){
            stop("`jitter` should be numeric!")
        }
        result <- quadraticRoots( 2*height*width/sqrt(3),
                                  (height+ sqrt(3)*width)/sqrt(3),
                                  0.5 - n_cells)
        ratio <- 1/result

        x_cells <- ceiling(width/ratio + 0.5)
        y_cells <- ceiling(n_cells/x_cells)

        n_sim <- x_cells*y_cells

        # adding randomness to the locations of the cells
        jitter_x <- stats::runif(n_sim, -jitter, jitter)
        jitter_y <- stats::runif(n_sim, -jitter, jitter)

        # assume the cells are located on hexagons but have little jittering
        x <- rep(seq(1, x_cells), y_cells) + jitter_x + c(rep(0,y_cells), rep(0.5, y_cells))
        y <- (rep(seq(1, y_cells), each = x_cells))*sqrt(3)/2 + jitter_y

        # format the image into an spe
        df <- data.frame(matrix(nrow = n_sim, ncol = 4))
        colnames(df)<-c("Cell.ID", "Cell.X.Position", "Cell.Y.Position", "Cell.Type")
        df$Cell.ID <- 1:n_sim
        df$Cell.X.Position <- x * ratio
        df$Cell.Y.Position <- y * ratio
        df$Cell.Type <- "Others"

        # delete the oversimulated cells
        n_delete <- n_sim - n_cells
        while(n_delete > 0 ){
            df_temp <- df[, c("Cell.X.Position", "Cell.Y.Position")]
            nn <- RANN::nn2(df_temp,df_temp, k = 2, radius = 20)
            id_delete <- which(nn[["nn.dists"]][,2] == min(nn[["nn.dists"]][,2]))[1]

            df <- df[-id_delete, ]
            n_delete <- n_delete - 1
        }
        if (plot_image) plot_cells(df, Cell.Type, "lightgray", "Cell.Type")
        sample <- df
    }

    else if (method == "Hardcore") {
        if(!is.numeric(min_d) | !is.numeric(oversampling_rate)){
            stop("One or more of `min_d`, `oversampling_rate` is not numeric!")
        }
        # need to oversample first
        n_cells_inflated <- n_cells*oversampling_rate

        # calculate the window and intensity
        win <- spatstat.geom::owin(xrange=c(0,width), yrange=c(0,height))
        beta <- n_cells_inflated/(width*height)

        # Hardcore process
        sample <- spatstat.random::rHardcore(beta = beta,R = min_d, W=win)

        # extract point data
        sample <- data.frame(Cell.X.Position = sample$x, Cell.Y.Position = sample$y)

        # if the sampled data is more than the expected number
        if (dim(sample)[1] > n_cells){
            sample <- sample[sample(nrow(sample), n_cells),]
        }

        rownames(sample) <- paste("Cell_",rownames(sample),sep = "")
        sample$Cell.Type <- Cell.Type

        # plot
        if (plot_image) plot_cells(sample, Cell.Type, "lightgray", "Cell.Type")
    }
    return(sample)
}
