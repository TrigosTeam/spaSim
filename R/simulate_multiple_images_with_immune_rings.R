#' simulate_multiple_images_with_immune_rings
#'
#' @param background_sample
#' @param cluster_size
#' @param ring_shape
#' @param infiltration
#' @param ring_width
#' @param cluster_loc_x
#' @param cluster_loc_y
#'
#' @return
#' @export
#'
#' @examples
simulate_multiple_images_with_immune_rings <- function(background_sample = bg1,
                                                       cluster_size = 200,
                                                       ring_shape = 1,
                                                       infiltration = 0,
                                                       ring_width = seq(50,100,10),
                                                       cluster_loc_x = 0,
                                                       cluster_loc_y = 0){
  ## CHECK
  # is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SingleCellExperiment::colData(background_sample))}

  # count the image number
  i <- 0
  list.images <- list()

  # choose a cluster shape (predefined in the package)
  if (ring_shape == 1){
    n_immune_rings <- 2
    properties_of_immune_rings_temp <- R_shape1
    # if the cluster size is too small, adjust the locations of some of the sub shapes
    size_threshold <- 100
  }
  else if(ring_shape == 2){
    n_immune_rings <- 2
    properties_of_immune_rings_temp <- R_shape2
    # if the cluster size is too small, adjust the locations of some of the sub shapes
    size_threshold <- 100
  }


  # loop through the properties
  for (size in cluster_size){
    # if the cluster size is too small, adjust the locations of some of the sub shapes
    if (size < size_threshold) {
    }
    for (infil in infiltration){
      for (width in ring_width){
        y_idx <- 0
        for (loc_x in cluster_loc_x) {
          i <- i + 1 # image count

          # find the corresponding y loc
          y_idx <- y_idx + 1
          loc_y <- cluster_loc_y[y_idx]

          # change the properties of the cluster based on the current loop
          properties_of_immune_rings <- properties_of_immune_rings_temp
          for (k in 1:length(properties_of_immune_rings)){
            # change the sizes of the default shape
            properties_of_immune_rings[[k]]$size <- properties_of_immune_rings[[k]]$size + size
            # change the default infiltration proportions
            properties_of_immune_rings[[k]]$infiltration_proportions[1] <- infil
            # change the default centre locations
            properties_of_immune_rings[[k]]$centre_loc[1] <- properties_of_immune_rings[[k]]$centre_loc[1] + loc_x
            properties_of_immune_rings[[k]]$centre_loc[2] <- properties_of_immune_rings[[k]]$centre_loc[2] + loc_y
            # change the immune ring width
            properties_of_immune_rings[[k]]$immune_ring_width <- width
          }

          # simulate the image
          image <- TIS(background_sample = background_sample,
                       n_immune_rings = n_immune_rings,
                       properties_of_immune_rings = properties_of_immune_rings,
          )
          plot_cell_categories(image, c("Tumour","Immune"),c("red","darkgreen"), "Phenotype")
          list.images[[i]] <- image
        }
      }
    }
  }
  return(list.images)
}
