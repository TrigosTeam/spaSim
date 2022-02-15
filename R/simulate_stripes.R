#' simulate_stripes
#'
#' @description Based on an existing background image, simulate stripes of cells
#'   representing vessels.
#'
#' @param background_sample (OPTIONAL) Data.Frame. The image that the stripes
#'   are simulated on. By default use the internal `bg1` background image.
#' @param n_stripe_type Number of stripe types.
#' @param win (OPTIONAL) `owin` object from spatstat.geom owin method. By
#'   default it is the window of the background image.
#' @param properties_of_stripes List of the properties of the stripes. See
#'   examples for the format of the properties.
#'
#' @family simulate pattern functions
#' @seealso   \code{\link{simulate_background_cells}} for all cell simulation,
#'   \code{\link{simulate_mixing}} for mixed background simulation,
#'   \code{\link{simulate_clusters}} for cluster simulation, and
#'   \code{\link{simulate_immune_rings}}/\code{\link{simulate_double_rings}} for
#'   immune ring simulation
#'
#' @return A data.frame of the simulated image
#' @export
#' @examples
#' properties_of_stripes = list(
#' S1 = list(
#'   number_of_stripes = 1,
#'   name_of_stripe_cell = "Others",
#'   width_of_stripe = 80,
#'   infiltration_types = c("Immune"),
#'   infiltration_proportions = c(0.08)
#' ), S2 = list(
#'   number_of_stripes = 1,
#'   name_of_stripe_cell = "Others",
#'   width_of_stripe = 80,
#'   infiltration_types = c("Immune"),
#'   infiltration_proportions = c(0.08)))
#' set.seed(610)
#' stripe_image <- simulate_stripes(background_sample = bg1, n_stripe_type=2,
#' win = NULL, properties_of_stripes = properties_of_stripes)
#'
#' # library(SPIAT)
#' # plot_cell_categories(stripe_image, categories_of_interest = c("Others","Immune"),
#' # colour_vector = c("gray", "darkgreen"), feature_colname = "Phenotype")

simulate_stripes <- function(background_sample = bg1,
                             n_stripe_type = 2,
                             win = NULL,
                             properties_of_stripes = list(
                               S1 = list(
                                 number_of_stripes = 1,
                                 name_of_stripe_cell = "Others",
                                 width_of_stripe = 80,
                                 infiltration_types = c("Immune"),
                                 infiltration_proportions = c(0.08)
                               ),
                               S2 = list(
                                 number_of_stripes = 1,
                                 name_of_stripe_cell = "Others",
                                 width_of_stripe = 80,
                                 infiltration_types = c("Immune"),
                                 infiltration_proportions = c(0.08)
                               )
                             )
){
  # CHECK is the background sample a dataframe?
  if (!is.data.frame(background_sample)) {
    background_sample <- data.frame(SingleCellExperiment::colData(background_sample))}
  # get the window
  if (is.null(win)) {
    X <- max(background_sample$Cell.X.Position)
    Y <- max(background_sample$Cell.Y.Position)
    win <- spatstat.geom::owin(c(0, X), c(0,Y))
  }


  # default phenotype is "Others"
  if (is.null(background_sample$Phenotype)){
    background_sample[, "Phenotype"] <- "Others"
  }

  for (k in 1:n_stripe_type){
    n_stripes = properties_of_stripes[[k]]$number_of_stripes
    stripe_cell_type = properties_of_stripes[[k]]$name_of_stripe_cell
    stripe_width = properties_of_stripes[[k]]$width_of_stripe
    infiltration_types = properties_of_stripes[[k]]$infiltration_types
    infiltration_proportions = properties_of_stripes[[k]]$infiltration_proportions

    # generate intercepts
    random_nums <- stats::runif(n_stripes, min = -max(X,Y), max = max(X,Y))

    for (i in 1:dim(background_sample)[1]){
      x <- background_sample[i, "Cell.X.Position"]
      y <- background_sample[i, "Cell.Y.Position"]
      pheno <- background_sample[i, "Phenotype"]

      p <- tryCatch(which(random_nums == max(random_nums[which(random_nums<y-x)])),
                    error=function(e) e, warning=function(w) w)

      if (methods::is(p,"warning") == FALSE) {
        b <- random_nums[p]
        if ( y < x + b + stripe_width ){
          random <- stats::runif(1)
          n_infiltration_types <- length(infiltration_types)
          pheno <- stripe_cell_type

          n <- 1 # start from the first proportion
          current_p <- 0
          while (n <= n_infiltration_types){
            current_p <- current_p + infiltration_proportions[n]
            if (random <= current_p) {
              pheno <- infiltration_types[n]
              break
            }
            n <- n+1
          }
          background_sample[i, "Phenotype"] <- pheno
        }
      }
    }
  }

  return(background_sample)
}
