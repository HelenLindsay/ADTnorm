# .adt_range ----
.adt_range <- function(adt_expr, range_extension=0.15){
  range_diff <- diff(range(adt_expr, na.rm = TRUE))
  from = min(adt_expr, na.rm = TRUE) - range_diff * range_extension
  to = max(adt_expr, na.rm = TRUE) + range_diff * range_extension
  return(c(from, to))
}

# .bw_by_zero_prop ----
.bw_by_zero_prop <- function(zero_prop){
  # different bandwidth w.r.t the zero proportion.
  if (zero_prop > 0.5) { return(3.1) }
  if (zero_prop > 0.3) { return(3) }
  return(2)
}

# .random_noise ----
.random_noise <- function(counts){
  return(counts + stats::rnorm(nrow(counts), mean = 0, sd = 0.05))
}

# .adjust_peak_indices ----
# If running via ADTnorm function as intended, validity of pos_samples has
# already been checked.
.adjust_peak_indices <- function(peak_locs, pos_samples=NULL, proximity=FALSE){
    n_peaks <- lengths(peak_locs)
    n_peaks[is.na(peak_locs)] <- 0
    max_peaks <- max(n_peaks)
    peak_align <- unlist(lapply(n_peaks, seq_len)) # idxs
    total_peaks <- rep(n_peaks, n_peaks) # totals corresponding to each idx
    peaks <- unlist(peak_locs)

    ## initiate landmark to record peak mode location
    landmark = matrix(NA, ncol = max_peaks, nrow = length(peak_locs))
    rownames(landmark) = names(peak_locs)

    # If we find >1 peaks, set the last peak to max_peaks
    peak_align[total_peaks > 1 & peak_align == total_peaks] <- max_peaks

    # If positive samples are explicitly specified, set to max idxs
    if (! is.null(pos_samples)){
        is_pos <- names(peak_locs) %in% pos_samples
        peak_align[is_pos & total_peaks == 1] <- max_peaks
    } else if (isTRUE(proximity) & any(total_peaks != max_peaks)){
        # Otherwise, set to closest average
        # TO DO: remove temp
        temp <- .adjust_by_proximity(peaks, total_peaks, max_peaks)
        if (! identical(temp, peak_align)){
            message("peak alignment changed by proximity")
            print(peak_align[peak_align != temp])
            print(temp[peak_align != temp])
        }
        peak_align <- temp
    }

    row_offset <- rep(seq_along(n_peaks), n_peaks)
    landmark[cbind(row_offset, peak_align)] <- peaks

    return(landmark)
}


# Average peaks across samples that have max_peaks, find closest peak
# in samples that have less
.adjust_by_proximity <- function(peaks, total_peaks, max_peaks){
    adjust <- total_peaks != max_peaks
    av_locs <- rowMeans(matrix(peaks[!adjust], nrow = max_peaks))
    print("adjust_by_proximity, av_locs")
    print(av_locs)
    peak_to_av <- abs(outer(peaks, av_locs, FUN = "-"))
    new_idxs <- apply(peak_to_av, 1, which.min)
    # CHECK IF THERE ARE DUPLICATES
    return(new_idxs)
}


# .all_negative_peaks ----
.all_negative_peaks <- function(landmark){
    landmarkAllMedian = stats::median(landmark[!is.na(landmark)])
    has_peak <- rowSums(! is.na(landmark)) > 0
    row_mins <- apply(abs(landmark[has_peak, , drop=FALSE] - landmarkAllMedian),
                      1, min, na.rm = TRUE)
    landmark[] <- NA
    landmark[has_peak] <- row_mins
    return(landmark[, 1, drop=FALSE])
}
