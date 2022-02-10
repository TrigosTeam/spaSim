#' simulate_background_cells
#'
#' @description Simulate cells without specifying their cell types. The 2D locations
#' of the cells are simulated and plotted in a rectangular window. Users can specify
#' the window size, cell number and the minimum distance between two cells.
#'
#' @param n_cells Numeric Number of cells in the background
#' @param width Numeric The width of the image
#' @param height Numeric The height of the image
#' @param min_d Numeric The minimum distance between two cells
#' @param oversample Numeric The multiplier for oversampling. Without oversampling,
#' the simulation deletes cells that are within min_d from each other, resulting
#' in a less number of cells than specified.
#' @param Phenotype String The name of the background cell type. Default is "Others"
#'
#' @importFrom ggplot2 ggplot aes geom_point
#' @return A data.frame of the simulated image
#' @export
#'
simulate_background_cells <- function(n_cells, width, height, min_d, oversample = 1.2,
                                      Phenotype = "Others"){

  # need to oversample first
  n_cells_inflated <- n_cells*oversample

  # calculate the window and intensity
  win <- spatstat.geom::owin(xrange=c(0,width), yrange=c(0,height))
  beta <- n_cells_inflated/(width*height)

  # Hardcore process
  sample <- spatstat.core::rHardcore(beta = beta,R = min_d, W=win)

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
  g <- ggplot(sample, aes(Cell.X.Position, Cell.Y.Position)) +
    geom_point()
  plot(g)

  return(sample)
}
