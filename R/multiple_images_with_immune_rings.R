#' Simulate multiple images with immune rings
#'
#' @description Generate a set of images with different immune ring properties.
#' @param bg_sample A data.frame or SingleCellExperiment class object
#'   with locations of points representing background cells. Further cell types
#'   will be simulated based on this background sample. The data.frame or the
#'   metadata of the SCE object should have colnames including
#'   "Cell.X.Positions" and "Cell.Y.Positions". By default use the internal
#'   \code{\link{bg1}} background image.
#' @param cluster_size Numeric Vector. The size of the cluster. If numeric, all
#'   simulated images have the same cluster size. If vector, images with a range
#'   of different cluster sizes will be simulated. The size should not exceed
#'   the limit of the image sides.
#' @param ring_shape Number. Choose from one of the following pre-designed
#'   shapes (1,2 or 3). The pre-designed shape contains information of the cell
#'   names of the cluster, the infiltration cell types, the proportions of
#'   infiltration, the cluster size, the ring width, the proportions of
#'   infiltrated cells into immune rings and the cluster centre locations. In
#'   order to simulate a set of images, use the arguments below to specify the
#'   ranges of the properties. The predefined cell types can not be changed,
#'   while users can change them manually after the simulation.
#' @param prop_infiltration Numeric Vector. The degree of infiltration. If
#'   numeric, all simulated images have the same infiltration degree. If vector,
#'   images with a range of different infiltration proportions will be
#'   simulated.
#' @param ring_width Numeric Vector. The width of the immune ring. If numeric,
#'   all simulated images have the same ring width. If vector, images with a
#'   range of different ring widths will be simulated.
#' @param cluster_loc_x Numeric or Vector. The X location of the cluster center
#'   offset. If numeric, all simulated images have the same center X location.
#'   If vector, images with a range of different center locations will be
#'   simulated.
#' @param cluster_loc_y Numeric or Vector of the same length of `cluster_loc_x`.
#'   The Y location of the cluster center offset.
#' @param prop_ring_infiltration Numeric or Vector. The degree of tumour
#'   infiltration in the region of immune rings.
#' @param plot_image Boolean Whether plot the simulated images or not.Default is
#'   TRUE.
#' @param plot_categories String Vector specifying the order of the cell
#'   cateories to be plotted.
#' @param plot_colours String Vector specifying the order of the colours that
#'   correspond to the `plot_categories` arg.
#'
#' @family simulate multiple images functions
#' @seealso \code{\link{multiple_background_images}} for simulating multiple
#'   mixed background images, and \code{\link{multiple_images_with_clusters}}
#'   for simulating multiple images with clusters.
#'
#' @return A list of sce objects
#' @export
#'
#' @examples
#' set.seed(610)
#' ring_image_list <- multiple_images_with_immune_rings(bg_sample = bg1,
#' ring_shape = 1, prop_infiltration = 0, ring_width = seq(50,100,10),
#' cluster_size = 300, cluster_loc_x = 0, cluster_loc_y = 0,
#' prop_ring_infiltration = seq(0, 0.2,0.05), plot_image = TRUE)

multiple_images_with_immune_rings <- function(bg_sample = bg1,
                                             cluster_size = 200,
                                             ring_shape = 1,
                                             prop_infiltration = 0,
                                             ring_width = seq(50,100,10),
                                             cluster_loc_x = 0,
                                             cluster_loc_y = 0,
                                             prop_ring_infiltration = seq(0, 0.2,0.05),
                                             plot_image = TRUE,
                                             plot_categories = NULL,
                                             plot_colours = NULL){
  
  ## CHECK
  # is the background sample a dataframe?
  if (!is.data.frame(bg_sample)) {
    bg_sample <- data.frame(SummarizedExperiment::colData(bg_sample))}

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
    for (infil in prop_infiltration){
      for (width in ring_width){
        y_idx <- 0
        for (loc_x in cluster_loc_x) {
          # find the corresponding y loc
          y_idx <- y_idx + 1
          loc_y <- cluster_loc_y[y_idx]
          for (ring_infil in prop_ring_infiltration){
            i <- i + 1 # image count

            # change the properties of the cluster based on the current loop
            properties_of_immune_rings <- properties_of_immune_rings_temp
            for (k in seq_len(length(properties_of_immune_rings))){
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
            image <- TIS(bg_sample = bg_sample,
                         n_immune_rings = n_immune_rings,
                         properties_of_immune_rings = properties_of_immune_rings,
                         plot_categories = plot_categories)
            
            if (plot_image){
              if(is.null(plot_categories)) plot_categories <- unique(image$Phenotype)
              if (is.null(plot_colours)){
                plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown", "purple", "lightblue",
                                  "lightgreen", "yellow", "black", "pink")}
              phenos <- plot_categories
              plot_cells(image, phenos, plot_colours[seq_len(length(phenos))], "Phenotype")
            }
            list.images[[i]] <- image
          }
        }
      }
    }
  }
  return(list.images)
}
