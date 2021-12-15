#' simulate_stripes
#'
#' @param background_sample Data.Frame The image that the stripes are simulated on
#' @param n_stripe_type Number of the types of the stripes
#' @param win (OPTIONAL) owin object from spatstat owin method. By default it is
#' the window of the background image.
#' @param properties_of_stripes List of the properties of the stripes
#'
#' @importFrom spatstat owin
#' @return
#' @export
simulate_stripes <- function(background_sample = sample2,
                             n_stripe_type,
                             win = NULL,
                             properties_of_stripes = list(
                               S1 = list(
                                 number_of_stripes = 1,
                                 name_of_stripe_cell = "Others",
                                 width_of_stripe = 80,
                                 infiltration_types = c("Immune1"),
                                 infiltration_proportions = c(0.08)
                               ),
                               S2 = list(
                                 number_of_stripes = 1,
                                 name_of_stripe_cell = "Others",
                                 width_of_stripe = 80,
                                 infiltration_types = c("Immune1"),
                                 infiltration_proportions = c(0.08)
                               )
                             )
){
  # get the window
  if (is.null(win)) {
    X <- max(background_sample$Cell.X.Position)
    Y <- max(background_sample$Cell.Y.Position)
    win <- owin(c(0, X), c(0,Y))
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
    random_nums <- runif(n_stripes, min = -max(X,Y), max = max(X,Y))

    for (i in 1:dim(background_sample)[1]){
      x <- background_sample[i, "Cell.X.Position"]
      y <- background_sample[i, "Cell.Y.Position"]
      pheno <- background_sample[i, "Phenotype"]

      p <- tryCatch(which(random_nums == max(random_nums[which(random_nums<y-x)])),
                    error=function(e) e, warning=function(w) w)

      if (is(p,"warning") == FALSE) {
        b <- random_nums[p]
        if ( y < x + b + stripe_width ){
          random <- runif(1)
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
