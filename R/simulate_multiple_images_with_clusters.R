#' simulate_multiple_images_with_clusters
#'
#' @param background_sample
#' @param cluster_shape
#' @param infiltration
#' @param cluster_size
#' @param cluster_loc_x
#' @param cluster_loc_y
#' @param plot.image
#'
#' @return
#' @export

simulate_multiple_images_with_clusters <- function(background_sample = bg1,
                                                   cluster_shape = 2,
                                                   infiltration = 0.1,
                                                   cluster_size = seq(200,1000,100),
                                                   cluster_loc_x = 0,
                                                   cluster_loc_y = 0,
                                                   plot.image = TRUE){
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
                     properties_of_clusters = properties_of_clusters,
        )
        plot_cell_categories(image, c("Tumour","Immune"),c("red","darkgreen"), "Phenotype")
        list.images[[i]] <- image
      }
    }
  }
  return(list.images)
}
