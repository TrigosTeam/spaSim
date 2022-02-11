test_that("simulate_clusters works", {
  bg_cluster <- simulate_clusters()

  # test if returns a data.frame
  expect_equal(class(bg_cluster), "data.frame")
  # test if "Tumour", "Immune1", "Immune2", "Others" exist under "Phenotype" column
  expect_setequal(unique(bg_cluster$Phenotype), c("Tumour", "Immune1", "Immune2", "Others"))
})

test_that("multiple_images_with_clusters works", {
  imageL <- multiple_images_with_clusters(background_sample = bg1,
                                          cluster_shape = 2,
                                          infiltration = 0.1,
                                          cluster_size = seq(400,500,100),
                                          cluster_loc_x = 0,
                                          cluster_loc_y = 0,
                                          plot.image = TRUE)
  sce <- imageL[[1]]

  # test if return a list of 2 objects
  expect_length(imageL, 2)
  # test if each object is an sce
  expect_equal(class(sce)[[1]], "SingleCellExperiment")
  # test if there are "Tumour" and "Immune", "Others" cells under the "Phenotype" column
  expect_setequal(colnames(colData(sce)),
                  c("Cell.X.Position", "Cell.Y.Position", "Phenotype", "pseudo"))
  expect_setequal(unique(sce$Phenotype), c("Tumour", "Immune", "Others"))

})

test_that("TIS works for simulating clusters", {
  image <- TIS(background_sample = bg1,
               n_clusters = 3,
               properties_of_clusters = C_shape1,
               image_name = "cluster_image")

  # test the class of the result
  expect_equal(class(image)[[1]], "SingleCellExperiment")

  # test if there are "Tumour" and "Immune", "Others" cells under the "Phenotype" column
  data <- data.frame(colData(image))
  expect_setequal(colnames(data),
                  c("Cell.X.Position", "Cell.Y.Position", "Phenotype", "pseudo"))
  expect_setequal(unique(image$Phenotype), c("Tumour", "Immune", "Others"))

  # test if the "name" attribute of the image is "cluster_image"
  expect_identical(attr(image, "name"), "cluster_image")
})
