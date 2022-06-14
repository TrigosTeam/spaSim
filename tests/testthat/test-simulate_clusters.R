test_that("simulate_clusters works", {
    bg_cluster <- simulate_clusters()

    # test if returns a data.frame
    expect_equal(class(bg_cluster), "data.frame")
    # test if "Tumour", "Immune1", "Immune2", "Others" exist under "Cell.Type" column
    expect_setequal(unique(bg_cluster$Cell.Type), c("Tumour", "Immune1", "Immune2", "Others"))
})

test_that("multiple_images_with_clusters works", {
    imageL <- multiple_images_with_clusters(bg_sample = bg1,
                                            cluster_shape = 2,
                                            prop_infiltration = 0.1,
                                            cluster_size = seq(400,500,100),
                                            cluster_loc_x = 0,
                                            cluster_loc_y = 0,
                                            plot_image = FALSE)
    spe <- imageL[[1]]

    # test if return a list of 2 objects
    expect_length(imageL, 2)
    # test if each object is an spe
    expect_equal(class(spe)[[1]], "SpatialExperiment")
    # test if there are "Tumour" and "Immune", "Others" cells under the "Cell.Type" column
    expect_setequal(colnames(SummarizedExperiment::colData(spe)), c("Cell.Type", "sample_id"))
    expect_setequal(colnames(SpatialExperiment::spatialCoords(spe)), c("Cell.X.Position", "Cell.Y.Position"))
    expect_setequal(unique(spe$Cell.Type), c("Tumour", "Immune", "Others"))

})

test_that("TIS works for simulating clusters", {
    image <- TIS(bg_sample = bg1,
                 n_clusters = 3,
                 properties_of_clusters = C_shape1,
                 image_name = "cluster_image")

    # test the class of the result
    expect_equal(class(image)[[1]], "SpatialExperiment")

    # test if there are "Tumour" and "Immune", "Others" cells under the "Cell.Type" column
    expect_setequal(colnames(SummarizedExperiment::colData(image)), c("Cell.Type", "sample_id"))
    expect_setequal(colnames(SpatialExperiment::spatialCoords(image)), c("Cell.X.Position", "Cell.Y.Position"))
    expect_setequal(unique(image$Cell.Type), c("Tumour", "Immune", "Others"))

    # test if the "name" attribute of the image is "cluster_image"
    expect_identical(attr(image, "name"), "cluster_image")
})
