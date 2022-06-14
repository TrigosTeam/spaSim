test_that("simulate_immune_rings works", {
    bg_immune_ring <- simulate_immune_rings()

    #test if returns a data.frame
    expect_equal(class(bg_immune_ring), "data.frame")
    # test if "Tumour", "Immune1", "Immune2", "Others" exist under "Cell.Type" column
    expect_setequal(unique(bg_immune_ring$Cell.Type),
                    c("Tumour", "Immune1", "Immune2", "Others"))
})

test_that("simulate_double_rings works", {
    double_ring <- simulate_double_rings()

    #test if returns a data.frame
    expect_equal(class(double_ring), "data.frame")
    # test if "Tumour", "Immune1", "Immune2", "Others" exist under "Cell.Type" column
    expect_setequal(unique(double_ring$Cell.Type),
                    c("Tumour", "Immune1", "Immune2", "Others"))
})

test_that("multiple_images_with_immune_rings works", {
    imageL <- multiple_images_with_immune_rings(bg_sample = bg1,
                                                cluster_size = 200,
                                                ring_shape = 1,
                                                prop_infiltration = 0,
                                                ring_width = seq(50,100,50),
                                                cluster_loc_x = 0,
                                                cluster_loc_y = 0,
                                                prop_ring_infiltration = seq(0, 0.2,0.2),
                                                plot_image = FALSE)
    spe <- imageL[[1]]

    # test if return a list of 4 objects
    expect_length(imageL, 4)
    # test if each object is an spe
    expect_equal(class(spe)[[1]], "SpatialExperiment")
    # test if there are "Tumour" and "Immune", "Others" cells under the "Cell.Type" column
    expect_setequal(colnames(SummarizedExperiment::colData(spe)), c("Cell.Type", "sample_id"))
    expect_setequal(colnames(SpatialExperiment::spatialCoords(spe)), c("Cell.X.Position", "Cell.Y.Position"))
    expect_setequal(unique(spe$Cell.Type), c("Tumour", "Immune", "Others"))

})

test_that("TIS works for simulating immune rings", {
    ir_properties <-list(
        I1 = list(
            name_of_cluster_cell = "Tumour",
            size = 600,
            shape = "Circle",
            centre_loc = data.frame("x" = 930, "y" = 1000),
            infiltration_types = c("Immune1", "Immune2", "Others"),
            infiltration_proportions = c(0.15, 0.05, 0.05),
            name_of_ring_cell = "Immune1",
            immune_ring_width = 150,
            immune_ring_infiltration_types = c("Others"),
            immune_ring_infiltration_proportions = c(0.15)
        ),
        I2 = list(
            name_of_cluster_cell = "Tumour",
            size = 500,
            shape = "Oval",
            centre_loc = data.frame("x" = 1330, "y" = 1100),
            infiltration_types = c("Immune1", "Immune2", "Others"),
            infiltration_proportions = c(0.15, 0.05, 0.05),
            name_of_ring_cell = "Immune1",
            immune_ring_width = 150,
            immune_ring_infiltration_types = c("Others"),
            immune_ring_infiltration_proportions = c(0.15)
        )
    )
    image <- TIS(bg_sample = bg1,
                 n_immune_rings = 2,
                 properties_of_immune_rings = ir_properties)

    # test the class of the result
    expect_equal(class(image)[[1]], "SpatialExperiment")
    # test if there are "Tumour" and "Immune", "Others" cells under the "Cell.Type" column
    expect_setequal(colnames(SummarizedExperiment::colData(image)), c("Cell.Type", "sample_id"))
    expect_setequal(colnames(SpatialExperiment::spatialCoords(image)), c("Cell.X.Position", "Cell.Y.Position"))
    expect_setequal(unique(image$Cell.Type), c("Tumour", "Immune1","Immune2", "Others"))
})


test_that("TIS works for simulating double rings", {
    image <- TIS(bg_sample = bg1,
                 n_double_rings = 2,
                 properties_of_double_rings = D_shape1)

    # test the class of the result
    expect_equal(class(image)[[1]], "SpatialExperiment")

    # test if there are "Tumour" and "Immune", "Others" cells under the "Cell.Type" column
    expect_setequal(colnames(SummarizedExperiment::colData(image)), c("Cell.Type", "sample_id"))
    expect_setequal(colnames(SpatialExperiment::spatialCoords(image)), c("Cell.X.Position", "Cell.Y.Position"))
    expect_setequal(unique(image$Cell.Type), c("Tumour", "Immune1", "Immune2", "Others"))
})
