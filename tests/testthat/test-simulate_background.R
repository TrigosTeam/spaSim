test_that("simulate_background_cells works", {
    set.seed(610)
    bg <- simulate_background_cells(n_cells = 5000,
                                    width = 2000,
                                    height = 2000,
                                    min_d = 10,
                                    oversampling_rate = 1.5,
                                    Cell.Type="Others")

    expect_identical(bg, bg1)
})

test_that("simulate_multiple_background_images works",{
    set.seed(610)
    imageL <- multiple_background_images(bg_sample = bg1,
                                         idents = c("Tumour","Immune", "Others"),
                                         props = list(
                                             rep(0.1, 3), seq(0.3, 0.4, 0.05),
                                             seq(0.6, 0.5, -0.05)),
                                         plot_image = FALSE)
    spe <- imageL[[1]]
    # test if return a list of 3 objects
    expect_length(imageL, 3)
    # test if each object is an spe
    expect_equal(class(imageL[[1]])[[1]], "SpatialExperiment")
    # test if there are expected information in the colData and spatialCoords
    expect_setequal(colnames(SummarizedExperiment::colData(spe)), c("Cell.Type", "sample_id"))
    expect_setequal(colnames(SpatialExperiment::spatialCoords(spe)), c("Cell.X.Position", "Cell.Y.Position"))
    # test if there are "Tumour" and "Immune", "Others" cells under the "Cell.Type" column
    expect_setequal(unique(spe$Cell.Type), c("Tumour", "Immune", "Others"))
})

test_that("TIS works for generating background image", {
    set.seed(610)
    image <- TIS(bg_sample = NULL,
                 n_cells = 5000,
                 width = 2000,
                 height = 2000,
                 min_d = 10,
                 names_of_bg_cells = c("Tumour","Immune", "Others"),
                 proportions_of_bg_cells = c(0.1, 0.3, 0.6),
                 plot_image = FALSE)
    # test the class of the result
    expect_equal(class(image)[[1]], "SpatialExperiment")
    # test if there are more than 4000 cells in the simulated image
    data <- data.frame(SummarizedExperiment::colData(image))
    expect_gt(dim(data)[1], 4000)
    # test if there are more than 200 Tumour cells
    expect_gt(dim(data[data$Cell.Type == "Tumour",])[1], 200)
})



