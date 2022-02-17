#' Simulate multiple images with clusters
#'
#' @description Generate a set of images with different cluster properties.
#' @param background_sample Data.frame or SingleCellExperiment object with
#'   locations of points representing background cells. Further cell types will
#'   be simulated based on this background sample.
#' @param cluster_shape List. Properties of clusters to simulate. Either choose
#'   one of the pre-designed shapes (1 or 2 for tumour cluster or 3 for immune
#'   cluster).
#' @param infiltration Numeric Vector. The degree of infiltration. If
#'   numeric, all simulated images have the same infiltration degree. If vector,
#'   images with a range of different infiltrations will be simulated.
#' @param cluster_size Numeric Vector. The size of the cluster. If numeric,
#'   all simulated images have the same cluster size. If vector, images with a
#'   range of different cluster sizes will be simulated.
#' @param cluster_loc_x Numeric or Vector. The X location of the cluster center.
#'   If numeric, all simulated images have the same center X location. If
#'   vector, images with a range of different center locations will be
#'   simulated.
#' @param cluster_loc_y Numeric or Vector of the same length of `cluster_loc_x`.
#'   The Y location of the cluster center.
#' @param plot.image Boolean Whether plot the simulated images or not.Default is
#'   TRUE.
#'
#' @family simulate multiple images functions
#' @seealso \code{\link{multiple_background_images}} for simulating multiple
#'   mixed background images, and
#'   \code{\link{multiple_images_with_immune_rings}} for simulating multiple
#'   images with immune rings.
#'
#' @return A list of sce objects
#' @export
#' @examples
#' set.seed(610)
#' cluster_image_list <- multiple_images_with_clusters(background_sample = bg1,
#' cluster_shape=2, infiltration = c(0.1, 0.3), cluster_size = c(300,600), cluster_loc_x = 0,
#' cluster_loc_y = 0, plot.image = TRUE)

multiple_images_with_clusters <- function(background_sample = bg1,
                                         cluster_shape = 2,
                                         infiltration = 0.1,
                                         cluster_size = seq(200,1000,100),
                                         cluster_loc_x = 0,
                                         cluster_loc_y = 0,
                                         plot.image = TRUE){
  ## CHECK
  # is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SummarizedExperiment::colData(background_sample))}

   # count the image number
  i <- 0
  list.images <- list()

  # choose a cluster shape (predefined in the package)
  if (cluster_shape == 1){
    n_clusters <- 3
    infiltration_types <- c("Immune", "Others")
    properties_of_clusters_temp <- C_shape1
    # if the cluster size is too small, adjust the locations of some of the sub shapes
    size_threshold <- 240
  }
  else if(cluster_shape == 2){
    n_clusters <- 2
    infiltration_types <- c("Immune", "Others")
    properties_of_clusters_temp <- C_shape2
    # if the cluster size is too small, adjust the locations of some of the sub shapes
    size_threshold <- 100
  }
  else if (cluster_shape == 3){
    n_clusters <- 1
    infiltration_types <- c("Immune", "Others")
    properties_of_clusters_temp <- C_shape3
    size_threshold <- 0
  }

  # loop through the cluster size and infiltration proportion ranges
  for (size in cluster_size){
    # if the cluster size is too small, adjust the locations of some of the sub shapes
    if (size < size_threshold) {
    }
    for (infil in infiltration){
      y_idx <- 0
      for (loc_x in cluster_loc_x) {
        i <- i + 1 # image count

        # find the corresponding y loc
        y_idx <- y_idx + 1
        loc_y <- cluster_loc_y[y_idx]

        # change the properties of the cluster based on the current loop
        properties_of_clusters <- properties_of_clusters_temp
        for (k in 1:length(properties_of_clusters)){
          # change the sizes of the default shape
          properties_of_clusters[[k]]$size <- properties_of_clusters[[k]]$size + size
          # change the default infiltration proportions
          properties_of_clusters[[k]]$infiltration_proportions[1] <- infil
          # change the default centre locations
          properties_of_clusters[[k]]$centre_loc[1] <- properties_of_clusters[[k]]$centre_loc[1] + loc_x
          properties_of_clusters[[k]]$centre_loc[2] <- properties_of_clusters[[k]]$centre_loc[2] + loc_y
        }

        # simulate the image
        image <- TIS(background_sample = background_sample,
                     n_clusters = n_clusters,
                     properties_of_clusters = properties_of_clusters)
        if (plot.image){
          plot_cells(image, c("Tumour","Immune"),c("red","darkgreen"), "Phenotype")
        }
        list.images[[i]] <- image
      }
    }
  }
  return(list.images)
}
