test_that("marker indices remain correct after subsetting", {
    data(cell_x_adt)
    data(cell_x_feature)

    tenx <- grepl("10X_pbmc_5k", cell_x_feature$batch)
    cell_x_adt_subs <- cell_x_adt[tenx, , drop = FALSE]
    cell_x_feature_subs <- cell_x_feature[tenx, , drop = FALSE]

    system.time(res <- ADTnorm(
        cell_x_adt = cell_x_adt_subs,
        cell_x_feature = cell_x_feature_subs,
        trimodal_marker = "CD45RA",
        marker_to_process = c("CD3", "CD14", "CD45RA"),
        save_intermediate_fig = FALSE
        ))



})
