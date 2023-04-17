process_log_file <- function(log_fname){
    f <- readr::read_delim(log_fname, show_col_types=FALSE)
    colnames(f) <- gsub(" ", "", colnames(f))
    f <- f %>%
      dplyr::mutate(diff_mode_midpoint = abs(peak_mode - peak_midpoint),
                    peak_width = peak_right - peak_left,
                    diff_on_range = diff_mode_midpoint/peak_width) %>%
      dplyr::group_by(sample, peak_n) %>%
      dplyr::mutate(final = n() == seq_along(peak_n)) %>%
      dplyr::ungroup()

  #    tidyr::pivot_wider(id_cols=c("sample", "bandwidth", "noise"),
  #                       names_from="peak_n",
  #                       values_from = c("peak_midpoint", "peak_mode",
  #                                       "diff_mode_midpoint"),
  #                       values_fill=NA)
}


# load_all()
# log_fname <- "~/Desktop/peak_midpoint_log.txt"
# peak_mode_res = get_peak_midpoint(cell_x_adt, cell_x_feature, log_fname, "CD19", 5)

#cell_x_adt_norm <- ADTnorm(
#  cell_x_adt = cell_x_adt,
#  cell_x_feature = cell_x_feature,
#  save_outpath = save_outpath,
#  study_name = run_name,
#  trimodal_marker = c("CD4", "CD45RA"),
#  positive_peak = list(ADT = "CD3", sample = "buus_2021_T"),
#  save_intermediate_fig = TRUE,
#  logdir = log_outpath
#)


#  marker_to_process = c("CD3", "CD4", "CD8", "CD45RA"),


# log_files <- list.files(log_outpath, full.names = TRUE)

# If positive peak, don't run through peak optimisation pipeline
