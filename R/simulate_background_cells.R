#' Simulate background cells
#'
#' @description Simulate cell locations. The 2D locations of the cells are
#'   simulated and plotted in a rectangular window. Users can specify the window
#'   size, cell number and the minimum distance between two cells. All cells
#'   have the same "Phenotype", specified by the "Phenotype" param. This
#'   function uses `rHardcore` [`spatstat.random`].
#'
#' @param n_cells Numeric. Number of cells to simulate in the background.
#' @param width,height Numeric. The width and height of the image.
#' @param min_d Numeric. The minimum distance between two cells.
#' @param oversampling_rate (OPTIONAL) Numeric. The multiplier for oversampling.
#'   Without oversampling, the simulation deletes cells that are within `min_d`
#'   from each other, resulting in a less number of cells than specified.
#'   Default is 1.2.
#' @param Phenotype (OPTIONAL) String. The name of the background cell type.
#'   Default is "Others".
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
#'                                               height = 2000, min_d = 10,
#'                                               oversampling_rate = 1.5, Phenotype = "Others")

simulate_background_cells <- function(n_cells, width, height, min_d,
                                      oversampling_rate = 1.2,
                                      Phenotype = "Others"){

    # need to oversample first
    n_cells_inflated <- n_cells*oversampling_rate

    # calculate the window and intensity
    win <- spatstat.geom::owin(xrange=c(0,width), yrange=c(0,height))
    beta <- n_cells_inflated/(width*height)

    # Hardcore process
    sample <- spatstat.random::rHardcore(beta = beta,R = min_d, W=win)

    # extract point data
    Cell.X.Position <- sample$x
    Cell.Y.Position <- sample$y
    sample <- data.frame(Cell.X.Position = sample$x, Cell.Y.Position = sample$y)

    # if the sampled data is more than the expected number
    if (dim(sample)[1] > n_cells){
        sample <- sample[sample(nrow(sample), n_cells),]
    }

    rownames(sample) <- paste("Cell_",rownames(sample),sep = "")

    sample$Phenotype <- Phenotype

    # plot
    plot_cells(sample, Phenotype, "lightgray", "Phenotype")

    return(sample)
}
