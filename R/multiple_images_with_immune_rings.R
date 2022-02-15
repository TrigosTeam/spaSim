#' multiple_images_with_immune_rings
#'
#' @param background_sample Data.frame or SingleCellExperiment object with locations
#' of points representing background cells. Further cell types will be simulated
#' based on this background sample.
#' @param cluster_size Numeric or Vector. The size of the cluster. If numeric,
#' all simulated images have the same cluster size. If vector, images with
#' a range of different cluster sizes will be simulated.
#' @param ring_shape List. Properties of immune rings to simulate. Either choose
#' one of the pre-designed shapes (1, 2, or 3), or manually
#' define the properties -- format:
#' R_shape3 = list(I1 = list(name_of_cluster_cell = "Tumour",size = 100,
#' shape = "Circle",centre_loc = data.frame("x" = 1000, "y" = 800),
#' infiltration_types = c("Immune", "Others"), infiltration_proportions = c(0, 0.1),
#' name_of_ring_cell = "Immune", immune_ring_width = 50, immune_ring_infiltration_types
#' = c("Tumour", "Others"), immune_ring_infiltration_proportions = c(0.3, 0.15)),
#' I2 = list(name_of_cluster_cell = "Tumour",size = 50,shape = "Oval",centre_loc =
#' data.frame("x" = 1100, "y" = 1050), infiltration_types = c("Immune", "Others"),
#' infiltration_proportions = c(0, 0.1), name_of_ring_cell = "Immune",
#' immune_ring_width = 50, immune_ring_infiltration_types = c("Tumour","Others"),
#' immune_ring_infiltration_proportions = c(0.3, 0.15)))
#' @param infiltration Numeric or Vector. The degree of infiltration.
#' @param ring_width Numeric or Vector. The width of the immune ring.
#' @param cluster_loc_x Numeric or Vector. The X location of the cluster center.
#' If numeric, all simulated images have the same center X location. If vector,
#' images with a range of different center locations will be simulated.
#' @param cluster_loc_y Numeric or Vector of the same length of `cluster_loc_x`.
#' The Y location of the cluster center.
#' @param ring_infiltration Numeric or Vector. The degree of tumour infiltration
#' in the region of immune rings. (Only applied to cell type "Tumour")
#' @param plot.image Boolean Whether plot the simulated images or not.Default is TRUE.
#'
#' @return A list of sce objects
#' @export
#' @examples
#' set.seed(610)
#' ring_image_list <- multiple_images_with_immune_rings(background_sample = bg1,
#' ring_shape = 1, infiltration = 0, ring_width = seq(50,100,10),
#' cluster_size = 300, cluster_loc_x = 0, cluster_loc_y = 0,
#' ring_infiltration = seq(0, 0.2,0.05), plot.image = TRUE)

multiple_images_with_immune_rings <- function(background_sample = bg1,
                                             cluster_size = 200,
                                             ring_shape = 1,
                                             infiltration = 0,
                                             ring_width = seq(50,100,10),
                                             cluster_loc_x = 0,
                                             cluster_loc_y = 0,
                                             ring_infiltration = seq(0, 0.2,0.05),
                                             plot.image = TRUE){
  ## CHECK
  # is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SingleCellExperiment::colData(background_sample))}

  # count the image number
  i <- 0
  list.images <- list()

  # choose a cluster shape (predefined in the package)
  if (ring_shape == 1){
    n_immune_rings <- 3
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
  else if(ring_shape == 3){
    n_immune_rings <- 2
    properties_of_immune_rings_temp <- R_shape3
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
          # find the corresponding y loc
          y_idx <- y_idx + 1
          loc_y <- cluster_loc_y[y_idx]
          for (ring_infil in ring_infiltration){
            i <- i + 1 # image count

            # change the properties of the cluster based on the current loop
            properties_of_immune_rings <- properties_of_immune_rings_temp
            for (k in 1:length(properties_of_immune_rings)){
              # change the sizes of the default shape
              properties_of_immune_rings[[k]]$size <-
                properties_of_immune_rings[[k]]$size + size
              # change the default infiltration proportions
              properties_of_immune_rings[[k]]$infiltration_proportions[1] <- infil
              # change the default centre locations
              properties_of_immune_rings[[k]]$centre_loc[1] <-
                properties_of_immune_rings[[k]]$centre_loc[1] + loc_x
              properties_of_immune_rings[[k]]$centre_loc[2] <-
                properties_of_immune_rings[[k]]$centre_loc[2] + loc_y
              # change the immune ring width
              properties_of_immune_rings[[k]]$immune_ring_width <- width
              # change the immune ring infiltration proportion
              properties_of_immune_rings[[k]]$immune_ring_infiltration_proportions[1] <-
                ring_infil
            }

            # simulate the image
            image <- TIS(background_sample = background_sample,
                         n_immune_rings = n_immune_rings,
                         properties_of_immune_rings = properties_of_immune_rings)
            if (plot.image){
              plot_cells(image, c("Tumour","Immune"),c("red","darkgreen"), "Phenotype")
            }
            list.images[[i]] <- image
          }
        }
      }
    }
  }
  return(list.images)
}
