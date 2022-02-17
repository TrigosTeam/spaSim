test_that("simulate_background_cells works", {
  set.seed(610)
  bg <- simulate_background_cells(n_cells = 5000,
                                  width = 2000,
                                  height = 2000,
                                  min_d = 10,
                                  oversample = 1.5,
                                  Phenotype="Others")

  expect_identical(bg, bg1)
})

test_that("simulate_multiple_background_images works",{
  set.seed(610)
  imageL <- multiple_background_images(background_sample = bg1,
                                       names_of_cell_types = c("Tumour","Immune", "Others"),
                                       proportions_of_cell_types = list(
                                         rep(0.1, 3), seq(0.3, 0.4, 0.05),
                                         seq(0.6, 0.5, -0.05)),
                                       plot.image = F)
  sce <- imageL[[1]]
  # test if return a list of 3 objects
  expect_length(imageL, 3)
  # test if each object is an sce
  expect_equal(class(imageL[[1]])[[1]], "SummarizedExperiment")
  # test if there are "Tumour" and "Immune", "Others" cells under the "Phenotype" column
  expect_setequal(colnames(colData(sce)),
                  c("Cell.X.Position", "Cell.Y.Position", "Phenotype", "pseudo"))
  expect_setequal(unique(sce$Phenotype), c("Tumour", "Immune", "Others"))
})

test_that("TIS works for generating background image", {
  set.seed(610)
  image <- TIS(background_sample = NULL,
              n_cells = 5000,
              width = 2000,
              height = 2000,
              min_d = 10,
              names_of_bg_cells = c("Tumour","Immune", "Others"),
              proportions_of_bg_cells = c(0.1, 0.3, 0.6))
  # test the class of the result
  expect_equal(class(image)[[1]], "SummarizedExperiment")
  # test if there are more than 4000 cells in the simulated image
  data <- data.frame(colData(image))
  expect_gt(dim(data)[1], 4000)
  # test if there are more than 200 Tumour cells
  expect_gt(dim(data[data$Phenotype == "Tumour",])[1], 200)
})



