test_that("ADTnorm works", {
    data(cell_x_adt)
    data(cell_x_feature)
    save_outpath <- tempdir()
    study_name <- "ADTnorm_demoRun"

    suppressWarnings({
      res <- ADTnorm(
        cell_x_adt = cell_x_adt,
        cell_x_feature = cell_x_feature,
        save_outpath = save_outpath,
        study_name = study_name,
        marker_to_process = c("CD3", "CD4", "CD8")
      )
    })

    expect_type(res, "double")
    expect_equal(nrow(res), 422682)
    expect_equal(ncol(res), 3)
})


test_that("incorrect inputs cause errors in ADTnorm", {
    data(cell_x_adt)
    data(cell_x_feature)

    adt_norm_args <- list(save_outpath = tempdir(),
                          study_name = "ADTnorm_demoRun",
                          cell_x_adt = cell_x_adt,
                          cell_x_feature = cell_x_feature,
                          save_outpath = save_outpath,
                          marker_to_process = c("CD3", "CD4", "CD8"))

    # Bimodal marker isn't a column in cell_x_adt
    expect_error(do.call(ADTnorm, utils::modifyList(adt_norm_args,
                                    list(bimodal_marker = "HELLO"))))
    # Trimodal marker isn't a column in cell_x_adt
    expect_error(do.call(ADTnorm, utils::modifyList(adt_norm_args,
                                           list(trimodal_marker = "HELLO"))))

    # invalid peak type
    expect_error(do.call(ADTnorm, utils::modifyList(adt_norm_args,
                                          list(peak_type = "Matterhorn"))))

    # invalid landmark
    expect_error(do.call(ADTnorm, utils::modifyList(adt_norm_args,
                                    list(landmark_align_type = "negpeak"))))

    # marker_to_process isn't in cell_x_adt
    expect_error(do.call(ADTnorm, utils::modifyList(adt_norm_args,
                              list(marker_to_process = c("CD4", "CD500")))))

    adt_norm_args$save_outpath <- NULL

    ## save_outpath is null and save_intermediate_rds is TRUE
    expect_error(do.call(ADTnorm,
                         utils::modifyList(adt_norm_args,
                                        list(save_intermediate_rds = TRUE))))

    ## save_outpath is null and save_intermediate_rds is TRUE
    expect_error(do.call(ADTnorm,
                         utils::modifyList(adt_norm_args,
                                           list(save_intermediate_fig = TRUE))))


    # Remove cell_x_adt and cell_x_feature
    adt_norm_args$save_outpath <- tempdir()
    adt_norm_args$cell_x_adt <- NULL
    adt_norm_args$cell_x_feature <- NULL

    c_x_adt <- cell_x_adt[1:1000, ]

    # cell_x_feature is missing
    expect_error(do.call(ADTnorm,
                         modifyList(adt_norm_args, list(cell_x_adt = c_x_adt))))

    # cell_x_adt and cell_x_feature have unequal dimensions
    expect_error(do.call(ADTnorm,
                         modifyList(adt_norm_args,
                                    list(cell_x_adt = c_x_adt,
                                         cell_x_feature = cell_x_feature))))
})
