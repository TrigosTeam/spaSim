#' simulate_background_cells
#'
#' @param n_cells Number of cells in the background
#' @param width Number The width of the image
#' @param height Number The height of the image
#' @param min_d Number The minimum distance between two cells
#' @param oversample Number The multiplier for oversampling
#'
#' @importFrom ggplots ggplot aes geom_point
#' @return
#' @export
#'
simulate_background_cells <- function(n_cells, width, height, min_d, oversample = 1.2){

  # need to oversample first
  n_cells_inflated <- n_cells*oversample

  # ccalculate the window and intensity
  win <- spatstat::owin(xrange=c(0,width), yrange=c(0,height))
  beta <- n_cells_inflated/(width*height)

  # Hardcore process
  sample <- spatstat::rHardcore(beta = beta,R = min_d, W=win)

  # extract point data
  Cell.X.Position <- sample$x
  Cell.Y.Position <- sample$y
  sample <- data.frame(Cell.X.Position = sample$x, Cell.Y.Position = sample$y)

  # if the sampled data is more than the expected number
  if (dim(sample)[1] > n_cells){
    sample <- sample[sample(nrow(sample), n_cells),]
  }

  rownames(sample) <- paste("Cell_",rownames(sample),sep = "")

  # plot
  g <- ggplot(sample, aes(Cell.X.Position, Cell.Y.Position)) +
    geom_point()
  plot(g)

  return(sample)
}
