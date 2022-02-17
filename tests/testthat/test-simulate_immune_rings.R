test_that("simulate_immune_rings works", {
  bg_immune_ring <- simulate_immune_rings()

  #test if returns a data.frame
  expect_equal(class(bg_immune_ring), "data.frame")
  # test if "Tumour", "Immune1", "Immune2", "Others" exist under "Phenotype" column
  expect_setequal(unique(bg_immune_ring$Phenotype),
                  c("Tumour", "Immune1", "Immune2", "Others"))
})

test_that("simulate_double_rings works", {
  double_ring <- simulate_double_rings()

  #test if returns a data.frame
  expect_equal(class(double_ring), "data.frame")
  # test if "Tumour", "Immune1", "Immune2", "Others" exist under "Phenotype" column
  expect_setequal(unique(double_ring$Phenotype),
                  c("Tumour", "Immune1", "Immune2", "Others"))
})

test_that("multiple_images_with_immune_rings works", {
  imageL <- multiple_images_with_immune_rings(background_sample = bg1,
                                              cluster_size = 200,
                                              ring_shape = 1,
                                              infiltration = 0,
                                              ring_width = seq(50,100,50),
                                              cluster_loc_x = 0,
                                              cluster_loc_y = 0,
                                              ring_infiltration = seq(0, 0.2,0.2),
                                              plot.image = F)
  sce <- imageL[[1]]

  # test if return a list of 2 objects
  expect_length(imageL, 4)
  # test if each object is an sce
  expect_equal(class(sce)[[1]], "SummarizedExperiment")
  # test if there are "Tumour" and "Immune", "Others" cells under the "Phenotype" column
  expect_setequal(colnames(colData(sce)),
                  c("Cell.X.Position", "Cell.Y.Position", "Phenotype","lab", "pseudo"))
  expect_setequal(unique(sce$Phenotype), c("Tumour", "Immune", "Others"))

})

test_that("TIS works for simulating immune rings", {
  image <- TIS(background_sample = bg1,
               n_immune_rings = 2,
               properties_of_immune_rings = R_shape2)

  # test the class of the result
  expect_equal(class(image)[[1]], "SummarizedExperiment")

  # test if there are "Tumour" and "Immune", "Others" cells under the "Phenotype" column
  data <- data.frame(colData(image))
  expect_setequal(colnames(data),
                  c("Cell.X.Position", "Cell.Y.Position", "Phenotype","lab", "pseudo"))
  expect_setequal(unique(image$Phenotype), c("Tumour", "Immune", "Others"))
})


test_that("TIS works for simulating double rings", {
  image <- TIS(background_sample = bg1,
               n_double_rings = 2,
               properties_of_double_rings = D_shape1)

  # test the class of the result
  expect_equal(class(image)[[1]], "SummarizedExperiment")

  # test if there are "Tumour" and "Immune", "Others" cells under the "Phenotype" column
  data <- data.frame(colData(image))
  expect_setequal(colnames(data),
                  c("Cell.X.Position", "Cell.Y.Position", "Phenotype","lab", "pseudo"))
  expect_setequal(unique(image$Phenotype), c("Tumour", "Immune1", "Immune2", "Others"))
})
