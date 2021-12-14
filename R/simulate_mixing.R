simulate_mixing <- function(background_sample = sample2,
                            names_of_mixing = c("Tumour", "Immune", "Others"),
                            shape = "Rectangle",
                            size = c(2500,2000),
                            mixing_degree = c(0.2, 0.4, 0.4),
                            centre_loc = data.frame("x" = c(1250), "y" = c(1000)),
                            win = owin(c(0,2500),c(0,2000))) {

  # default phenotype is "Others"
  if (is.null(background_sample$Phenotype)){
    background_sample[, "Phenotype"] <- "Others"
  }

  # generate a location as the centre of the cluster
  if (is.null(centre_loc)){
    seed_point <- runifpoint(1, win=win)}
  else seed_point <- centre_loc
  a <- seed_point$x
  b <- seed_point$y

  n_types <- length(names_of_mixing)

  # size of the rectangle
  w <- 0.5*size[1]
  h <- 0.5*size[2]

  for (i in 1:dim(background_sample)[1]){
    x <- background_sample[i, "Cell.X.Position"]
    y <- background_sample[i, "Cell.Y.Position"]

    if (x > a - w && x < a + w && y > b - h && y < b + h){
      random <- runif(1)

      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_types){
        current_p <- current_p + mixing_degree[n]
        if (random <= current_p) {
          pheno <- names_of_mixing[n]
          break
        }
        n <- n+1
      }

      background_sample[i, "Phenotype"] <- pheno
    }
  }
  return(background_sample)
}
