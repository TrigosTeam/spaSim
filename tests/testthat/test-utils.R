test_that("format_spe works", {
    spe <-  format_spe(bg1)

    #test if returns an spe
    expect_equal(class(spe)[[1]], "SpatialExperiment")
    # test if spatial coordinates are there
    expect_equal(colnames(SpatialExperiment::spatialCoords(spe)),
                 c("Cell.X.Position", "Cell.Y.Position"))
    # test if all cells are there
    expect_setequal(unique(SummarizedExperiment::colData(spe)$Cell.Type),  "Others")
})


test_that("get_colData works", {
    spe <-  format_spe(bg1)
    df <- get_colData(spe)

    #test if returns a data.frame
    expect_equal(class(df), "data.frame")
    # test if "Tumour", "Immune1", "Immune2", "Others" exist under "Cell.Type" column
    expect_setequal(unique(df$Cell.Type),"Others")

    # test if df is the same as bg1
    df <- tibble::column_to_rownames(df, "Cell.ID")
    df <- df[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Type" )]
    expect_equal(df, bg1)
})
