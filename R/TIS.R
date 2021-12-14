##### function #####
TIS <- function(background_sample = NULL,
                n_cells = NULL,
                width = NULL,
                height = NULL,
                min_d = NULL,

                names_of_bg_cells = NULL,
                proportions_of_bg_cells = NULL,

                n_clusters = NULL,
                properties_of_clusters = NULL,

                n_immune_rings = NULL,
                properties_of_immune_rings = NULL,

                n_stripe_type = NULL,
                properties_of_stripes = NULL,

                image_name = NULL)
{
  if (is.null(background_sample)){
    background_sample <- simulate_background_cells(n_cells, width, height, min_d)
    X <- width
    Y <- height
  }

  image <- background_sample
  X <- max(background_sample$Cell.X.Position)
  Y <- max(background_sample$Cell.Y.Position)

  # get background information
  bg_size = c(X, Y)

  # simulate bg with mixing types of cells
  if (!is.null(names_of_bg_cells)){
    image <- simulate_mixing(background_sample = image,
                             names_of_mixing = names_of_bg_cells,
                             mixing_degree = proportions_of_bg_cells,
                             shape = "Rectangle",
                             size = bg_size,
                             centre_loc = data.frame("x" = X/2, "y" = Y/2),
                             win = owin(c(0,X),c(0,Y)))
  }
  # simulate clusters
  if (!is.null(n_clusters)){
    image <- simulate_clusters(background_sample = image,
                               n_clusters = n_clusters,
                               bg_type = "Others",
                               win = NULL,
                               properties_of_clusters = properties_of_clusters)
  }
  # simulate_immune_rings
  if (!is.null(n_immune_rings)){
    image <- simulate_immune_rings(background_sample = image,
                                   n_immune_rings = n_immune_rings,
                                   bg_type = "Others",
                                   win = NULL,
                                   properties_of_immune_rings = properties_of_immune_rings)
  }
  # simulate_stripes
  if (!is.null(n_stripe_type)){
    image <- simulate_stripes(background_sample = image,
                              n_stripe_type = n_stripe_type,
                              win = NULL,
                              properties_of_stripes = properties_of_stripes)}

  # format sce object
  sce <- format_colData_to_sce(image)
  attr(sce, "name") <- image_name
  return(sce)
}
