# .adjust_peak_indices ----
test_that(".adjust_peak_indices works as expected", {
    peak_locs <- list(A = c(0.5, 1),
                      B = c(0.55, 1.05),
                      C = c(0.1, 0.49, 1.2),
                      D = c(0.4))
    expected_res <- matrix(c(0.5, NA, 1,
                             0.55, NA, 1.05,
                             0.1, 0.49, 1.2,
                             0.4, NA, NA), byrow = TRUE, nrow = 4,
                           dimnames = list(c("A", "B", "C","D"), NULL))

    res <- .adjust_peak_indices(peak_locs)
    expect_equal(res, expected_res)

    # If we set a positive peak, it should be adjusted to the third position
    res <- .adjust_peak_indices(peak_locs, pos_samples = "D")

    expected_res_adj <- matrix(c(0.5, NA, 1,
                                 0.55, NA, 1.05,
                                 0.1, 0.49, 1.2,
                                 NA, NA, 0.4), byrow = TRUE, nrow = 4,
                               dimnames = list(c("A", "B", "C","D"), NULL))
    expect_equal(res, expected_res_adj)


    # If we set a positive peak but there is more than one peak,
    # nothing should happen
    res <- .adjust_peak_indices(peak_locs, pos_samples = "C")
    expect_equal(res, expected_res)
})

# .all_negative_peaks ----
test_that(".all_negative_peaks works as expected", {
    landmarks <- matrix(c(0.5, NA, 1,
                          0.55, NA, 1.05,
                          0.1, 0.49, 1.2,
                          NA, NA, NA), byrow = TRUE, nrow = 4,
                        dimnames = list(c("A", "B", "C","D"), NULL))

    exp_landmarks <- matrix(c(0.05, 0, 0.06, NA), nrow = 4,
                            dimnames = list(c("A", "B", "C","D"), NULL))
    expect_equal(.all_negative_peaks(landmark), exp_landmarks)
})

# test that peak alignment works if there is one extra positive peak
