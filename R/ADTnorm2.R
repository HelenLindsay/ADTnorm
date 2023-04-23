ADTnorm2 <- function(cell_x_adt, cell_x_feature, save_outpath=NULL,
                    study_name="ADTnorm", marker_to_process=NULL,
                    positive_peak=NULL, trimodal_marker=NULL,
                    bw_smallest_bi=1.1, bw_smallest_tri=0.8,
                    peak_type="midpoint",
                    multi_sample_per_batch=FALSE, shoulder_valley=FALSE,
                    shoulder_valley_slope=-0.5, valley_density_adjust=3,
                    landmark_align_type="negPeak_valley_posPeak",
                    midpoint_type="valley", neg_candidate_thres=asinh(8/5 + 1),
                    lower_peak_thres=0.001, brewer_palettes="Set1",
                    save_intermediate_rds=FALSE, save_intermediate_fig=TRUE,
                    detect_outlier_valley=FALSE, target_landmark_location=NULL,
                    clean_adt_name=FALSE, verbose=FALSE, log_dir=NULL,
                    proximity=TRUE)
{

    if (isTRUE(clean_adt_name)) { # standardise names ----
        new_nms <- clean_adt_name(colnames(cell_x_adt))
        marker_to_process <- new_nms[match(marker_to_process,
                                           colnames(cell_x_adt))]
        colnames(cell_x_adt) = new_nms
    }

    cell_x_feature <- .format_cell_x_feature(cell_x_feature)

    # Check validity of inputs ----
    .checkInputsADTnorm(cell_x_adt, cell_x_feature, trimodal_marker,
                        save_intermediate_rds, save_intermediate_fig,
                        save_outpath, peak_type, landmark_align_type,
                        marker_to_process, multi_sample_per_batch,
                        positive_peak)

    # Set target locations ----
    if (! is.null(target_landmark_location)){
        target_landmark_location <- .setTargetLocation(target_landmark_location)
    }

    # make arcsinh_transform optional ----
    #if(input_raw_counts){
    #    ## Check that cell_x_adt is integers only
    #    matrix = as.matrix(cell_x_adt)
    #    as_int = as.matrix(cell_x_adt)
    #    mode(as_int) <- "integer"
    #    if(!all(as_int[!is.na(as_int)] == matrix[!is.na(matrix)]) | !all(matrix[!is.na(matrix)]>=0)){
    #        stop("When using input_raw_counts, please only input positive",
    #             "integers for expression values.")
    #    }
    #
    #    cell_x_adt = arcsinh_transform(cell_x_adt = cell_x_adt) ## Arcsinh transformation
    #} ----

    cell_x_adt = arcsinh_transform(cell_x_adt = cell_x_adt)

    all_marker_name = colnames(cell_x_adt)

    # Select peak function
    get_peak <- ifelse(peak_type == "mode", get_peak_mode, get_peak_midpoint)

    bw_res <- .setBandwidth(all_marker_name, trimodal_marker,
                            bw_smallest_bi, bw_smallest_tri)

    cell_x_adt_norm <- .selectMarkers(cell_x_adt, marker_to_process)

    #######
    # TO DO: n_expected_peak
    # Do the subsetting in this function
    #######

    if (! is.null(log_dir)) dir.create(log_dir, recursive=TRUE)

    for (adt in all_marker_name) {
       print(adt)

       ######################
       if (! is.null(log_dir)){
          log_file <- file.path(log_dir, adt)
       } else {
          log_file <- NULL
       }
       #####################

       bwFac_smallest <- bw_res$bwFac_smallest[[adt]]
       marker_type <- bw_res$marker_type[[adt]]

       peak_mode_res = get_peak(cell_x_adt, cell_x_feature, log_file,
                                adt, marker_type,
                                bwFac_smallest=bwFac_smallest,
                                positive_peak=positive_peak,
                                neg_candidate_thres=neg_candidate_thres,
                                lower_peak_thres=lower_peak_thres,
                                proximity=proximity)

       # get valley ----
       peak_valley_list <- get_valley_location(
           cell_x_adt, cell_x_feature, adt, peak_mode_res,
           shoulder_valley, positive_peak, multi_sample_per_batch,
           adjust=valley_density_adjust, min_fc=20,
           shoulder_valley_slope=shoulder_valley_slope,
           neg_candidate_thres=neg_candidate_thres) # proximity=TRUE

       valley_location_res <- peak_valley_list$valley_location_list
       peak_mode_res <- peak_valley_list$peak_landmark_list

    print("stopping")
    return(peak_mode_res)

    if (detect_outlier_valley) {
      valley_location_res <- detect_impute_outlier_valley(
          valley_location_res,
          adt_marker_select, cell_x_adt, cell_x_feature,
          scale = 3, method = "MAD", nearest_neighbor_n = 3,
          nearest_neighbor_threshold = 0.75)
    }

    # Plot and save ----

    n_samples = length(levels(cell_x_feature$sample))
    adt_study = sprintf("%s_%s", adt_marker_select_name, study_name)
    fig_dir <- file.path(save_outpath, "/figures")
    rds_dir = file.path(save_outpath, "/RDS")

    # fig_dir needed for this function?
    density_plot <- plot_adt_density_with_peak_valley_each(
        cell_x_adt[, adt_marker_select], cell_x_feature,
        peak_landmark_list = peak_mode_res,
        valley_landmark_list = valley_location_res,
        brewer_palettes = brewer_palettes,
        parameter_list = list(bw = 0.1,
                              method_label = "Arcsinh Transformation"))

    peak_valley <- list(peak_landmark_list = peak_mode_res,
                        valley_landmark_list = valley_location_res)

    if (isTRUE(save_intermediate_rds)) {
        if ( ! dir.exists(rds_dir)) {
            dir.create(rds_dir, recursive = TRUE)
        }
        saveRDS(peak_valley,
                file = sprintf("%s/peak_valley_raw_%s.rds", rds_dir, adt_study),
                compress = FALSE)
        saveRDS(density_plot,
                file = sprintf("%s/density_raw_%s.rds", rds_dir, adt_study),
                compress = FALSE)
    }

    if (isTRUE(save_intermediate_fig)) {
        if (! dir.exists(fig_dir) ) {
            dir.create(fig_dir, recursive = TRUE)
        }
        fig_height <- ceiling(n_samples * 0.4)
        grDevices::pdf(sprintf("%s/ArcsinhTransform_%s.pdf", fig_dir, adt_study),
                       width = 11, height = fig_height)
        print(density_plot)
        grDevices::dev.off()
    }


    # Landmark location ----
    landmark_matrix <- landmark_fill_na(peak_landmark_list=peak_mode_res,
                                    valley_landmark_list=valley_location_res,
                                    landmark_align_type=landmark_align_type,
                                    midpoint_type=midpoint_type,
                                    neg_candidate_thres=neg_candidate_thres)

    # what is the difference between target_landmark_location and peak_mode_res?
    target_landmark <- .setTargetLandmark(target_landmark_location,
                                          landmark_align_type)


    # Peak alignment resolution ----
    peak_alignment_res = peak_alignment(cell_x_adt[, adt_marker_select],
                                        cell_x_feature, landmark_matrix,
                                        target_landmark = target_landmark)

    cell_x_adt_norm[, adt_marker_select] = peak_alignment_res[[1]]

    if (ncol(peak_alignment_res[[2]]) == 2) {
      peak_mode_norm_res = t(t(peak_alignment_res[[2]][, 1]))
      valley_location_norm_res = t(t(peak_alignment_res[[2]][,2]))
    }
    else if (ncol(peak_alignment_res[[2]]) == 3) {
      peak_mode_norm_res = peak_alignment_res[[2]][, c(1,3)]
      valley_location_norm_res = t(t(peak_alignment_res[[2]][,2]))
    }
    else if (ncol(peak_alignment_res[[2]]) == 5) {
      peak_mode_norm_res = peak_alignment_res[[2]][, c(1, 3, 5)]
      valley_location_norm_res = peak_alignment_res[[2]][, c(2, 4)]
    }


    # Normalised density plot ----
    density_norm_plot <- plot_adt_density_with_peak_valley_each(
        cell_x_adt_norm[, adt_marker_select], cell_x_feature,
        peak_landmark_list=peak_mode_norm_res,
        valley_landmark_list=valley_location_norm_res,
        brewer_palettes=brewer_palettes,
        parameter_list=list(bw = 0.2, method_label = "ADTnorm"))

    if (isTRUE(save_intermediate_rds)) {
        saveRDS(density_norm_plot,
                file = sprintf("%s/density_ADTnorm_%s.rds", rds_dir, adt_study),
                compress = FALSE)
    }
    if (isTRUE(save_intermediate_fig)) {
        grDevices::pdf(sprintf("%s/ADTnorm_%s.pdf", fig_dir, adt_study),
                       width = 11, height = ceiling(n_samples * 0.4))
        print(density_norm_plot)
        grDevices::dev.off()
    }
  }


  colnames(cell_x_adt_norm) = all_marker_name[adt_marker_index_list]
  return(cell_x_adt_norm)
}


# .checkInputsADTnorm ----
.checkInputsADTnorm <- function(cell_x_adt, cell_x_feature, trimodal_marker,
                                save_intermediate_rds, save_intermediate_fig,
                                save_outpath, peak_type, landmark_align_type,
                                marker_to_process, multi_sample_per_batch,
                                positive_peak){

    if (nrow(cell_x_adt) < ncol(cell_x_adt)) {
        warning("Please check if the ADT raw count matrix has cell as row and ",
                "adt marker as column.")
    }

    # Check that the cell_x_adt and cell_x_feature have the same number of cells
    if (! nrow(cell_x_adt) == nrow(cell_x_feature)){
        stop("rows (cells) in cell_x_feature should match rows in cell_x_adt")
    }

    if (isTRUE(multi_sample_per_batch) & !
        ("batch" %in% colnames(cell_x_feature))){
        stop("cell_x_feature must contain a column named batch if",
             "multi_sample_per_batch is TRUE")
    }

    if (sum( !(trimodal_marker %in% colnames(cell_x_adt))) > 0) {
        stop("Trimodal marker names must match column names of the input ",
             "ADT count matrix.")
    }

    # If results should be saved, save_outpath must not be null ----
    if (is.null(save_outpath)){
        msg <- "Please provide the save_outpath to save the intermediate "
        if (isTRUE(save_intermediate_rds)) { stop(msg, "results as rds.") }
        if (isTRUE(save_intermediate_fig)){ stop(msg, "figures as pdf.") }
    }

    if (! (peak_type %in% c("mode", "midpoint"))) {
       stop("Please specify the peak type to be either 'mode' or 'midpoint'.")
    }

    allowed_landmarks = c("negPeak", "negPeak_valley",
                          "negPeak_valley_posPeak", "valley")

    if (! landmark_align_type %in% allowed_landmarks) {
        stop("Please provide one of the landmark_align_type from: ",
             toString(allowed_landmarks))
    }

    # Check that the markers to process are in the cell_x_adt matrix
    if (! is.null(marker_to_process)){
        missing_marker <- setdiff(marker_to_process, colnames(cell_x_adt))
        if (length(missing_marker) > 0){
            stop("The following marker_to_process are not columns of ",
                 "cell_x_adt:", paste(missing_marker, sep = ", "))
        }
    }

    # Check that all sample names in positive_peak are in cell_x_feature
    pos_samples <- unique(unname(unlist(positive_peak)))
    if (! all(pos_samples %in% cell_x_feature$sample)){
        stop("All sample names in 'positive_peak' must be in the column
             'sample' of 'cell_x_feature'")
    }
}


# .format_cell_x_feature ----
.format_cell_x_feature <- function(cell_x_feature){
    # Check that cell_x_feature has the required columns
    if (! "sample" %in% colnames(cell_x_feature)){
        stop("cell_x_feature must contain a column named sample")
    }

    # Make cell_x_feature a factor if not already
    if (! is.factor(cell_x_feature$sample) ){
        cell_x_feature$sample = factor(cell_x_feature$sample,
                                       levels = unique(cell_x_feature$sample))
    }

    # Droplevels in case this is a subset
    cell_x_feature$sample <- droplevels(cell_x_feature$sample)

    return(cell_x_feature)
}


# .setTargetLocation ----
.setTargetLocation <- function(landmark_loc){
    msg <- "Will align negative peak to %s and right-most positive peak to %s"
    if (landmark_loc == "fixed") {
        landmark_loc = c(1, 5)
    } else if (! length(landmark_loc) == 2 | landmark_loc[2] < landmark_loc[1]){
            stop("Please provide two elements vector to ",
                 "target_landmark_location where the first element is smaller!")
    }

    message(sprintf(msg, landmark_loc[1], landmark_loc[2]))
    return(landmark_loc)
}


# .setBandwidth ----
.setBandwidth <- function(marker_names, trimodal_marker, bw_smallest_bi,
                          bw_smallest_tri, special_cases=c("CD4","CD3","CD8"),
                          special_bw=0.8){

    n <- length(marker_names)
    bwFac_smallest <- structure(rep(bw_smallest_bi, n), names = marker_names)
    marker_type <- structure(rep("bimodal", n), names = marker_names)

    bwFac_smallest[marker_names %in% trimodal_marker] <- bw_smallest_tri
    marker_type[marker_names %in% trimodal_marker] <- "trimodal"

    bwFac_smallest[marker_names %in% special_cases] <- special_bw
    marker_type[match("CD4", marker_names)] <- "CD4"

    return(list(bwFac_smallest = bwFac_smallest, marker_type = marker_type))
}


# .selectMarkers ----
# Subset cell_x_adt and report which markers will be used
.selectMarkers <- function(cell_x_adt, marker_to_process){
    msg <- "ADTnorm will process "
    if (is.null(marker_to_process)) {
        message(sprintf("%s all the ADT markers from the ADT matrix: %s\n",
                        msg, toString(colnames(cell_x_adt))))
        cell_x_adt_norm = cell_x_adt
    }

    else {
        # Note: checked that all marker_to_process are in colnames cell_x_adt
        messsage(sprintf("%s the following ADT markers as provided: ",
                         msg, paste(marker_to_process, collapse = ", "), "\n"))
        cell_x_adt_norm = cell_x_adt[, marker_to_process, drop = FALSE]
    }
    return(cell_x_adt_norm)
}


# .setTargetLandmark ----
.setTargetLandmark <- function(peaks, landmark_align_type){
    if (is.null(peaks)) { return(NULL) }

    if (landmark_align_type == "negPeak") { return(peaks[1]) }

    mean_loc <- round(mean(peaks), 1)

    if (landmark_align_type == "valley"){ return(mean_loc) }
    if (landmark_align_type == "negPeak_valley") {
        return(c(peaks[1], mean_loc))
    }
    if (landmark_align_type == "negPeak_valley_posPeak") {
        if (ncol(peak_mode_res) == 1) { return(c(peaks[1], mean_loc)) }
        if (ncol(peak_mode_res) == 2) {
            return(c(peaks[1], mean_loc, peaks[2]))
        }
        return(c(peaks[1], (peaks[1] + mean_loc)/2, peaks[2]))
    }

    return(NULL)
}


# .setPeakAlignRes ----
.setPeakAlignRes <- function(){
    # Check that it cannot be 1

    # Peak alignment resolution ----
    pk_align = peak_alignment(cell_x_adt[, adt_marker_select], cell_x_feature,
                              landmark_matrix, target_landmark=target_landmark)

    cell_x_adt_norm[, adt_marker_select] = pk_align[[1]]

    # Odd number columns are peaks, even numbers are valleys
    col_idx <- seq_len(ncol(pk_align[[2]]))
    peak_mode_norm_res = pk_align[[2]][, col_idx %% 2 == 1, drop = FALSE]
    valley_location_norm_res = pk_align[[2]][, col_idx %% 2 == 0, drop = FALSE]

    return(list(peak_mode_norm_res=peak_mode_norm_res,
                valley_location_norm_res=valley_location_norm_res))
}
