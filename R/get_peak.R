# TO DO: DEFINITION OF positive_peak changed!
# cell_x_feature and cell_x_adt pre-subset for sample

#' Get the peak landmark locations
#'
#' This function detects the peak landmark locations for each sample per ADT
#' marker using either the peak mode or midpoint. Using the peak midpoint can
#' be more stable across samples and less affected by the bandwidth.
#' Using the peak mode can be more accurate in determining the peak location if
#' the band width is generally proper and the local peak density is not too
#' discrete.
#'
#' @param expr Single-column matrix of arcsinh transformed ADT counts per cell
#' @param samples vector of sample names, matching adt
#' @param bwFac_smallest The smallest band width parameter value.
#' Recommend 1.1 for general bi-modal ADT markers except CD3, CD4 and CD8.
#' @param n_expected_peaks How many peaks are expected based on researchers'
#' prior knowledge (e.g. CD4 usually has 3 peaks) or preliminary observation on
#' particular data to be processed.  Should be 2 or 3.
#' @param positive_peak A vector of sample name(s) to specify that when a
#' single peak is detected it should be aligned to positive peaks. For example,
#' for samples that only contain T cells, the only CD3 peak should be aligned
#' to positive peaks of other samples.
#' @param neg_candidate_thres The upper bound for the negative peak. Users can
#' refer to their IgG samples to obtain the minimal upper bound of the IgG
#' sample peak. It can be one of the values of asinh(4/5+1), asinh(6/5+1),
#' or asinh(8/5+1) if the right 95% quantile of IgG samples are large.
#' @param lower_peak_thres The minimal ADT marker density height to call it a
#' real peak. Set it to 0.01 to avoid suspicious positive peaks.
#' Set it to 0.001 or smaller to include some small but often real
#' positive peaks, especially for markers like CD19.
#' @import dplyr
#' @export
#' @examples
#' \dontrun{
#' get_peak(
#'   expr = cell_x_adt["CD3",],
#'   samples = cell_x_feature$sample,
#'   neg_candidate_thres = asinh(6/5+1)
#' )
#' }
get_peak = function(expr, samples, n_expected_peaks = 2,
                    peak_type = c("mode", "midpoint"),
                    bwFac_smallest = 1.1, positive_peak = NULL,
                    neg_candidate_thres = asinh(10/5 + 1),
                    lower_peak_thres = 0.001) {

    if (ncol(adt) > 1){ stop("expr should be a single column matrix") }
    if (! is.factor(sample)){ stop("sample should be a factor") }
    expr <- as.matrix(expr)

    border <- 0.01
    peak_type <- match.arg(peak_type)
    sample_names <- levels(samples)
    expr <- as.matrix(expr)

    peak_loc = list()

    # get ADT value range with a slight extension across all samples
    adt_range <- .adt_range(expr)

    # get peak for each sample of this processing ADT marker
    for (sample_name in sample_names){

        keep_cells <- sample == sample_name & ! is.na(expr)
        if (! any(keep_cells)){
          peak_loc[[sample_name]] = NA
          next
        }

        fcs_count = expr[keep_cells, , drop=FALSE]
        n_unique_vals = nrow(unique(fcs_count))

        if (n_unique_vals == 1){
            ## only one value for this marker
            peak_loc[[sample_name]] = NA
            message(sample_name, "-Single Value!")
            next
        }

        # get proportion of cells near zero to diagnose the neg peak enrichment
        zero_prop = sum(fcs_count < 2) / nrow(fcs_count)

        # if most are around 0 and there are very few unique values:
        # add random small number
        if (zero_prop > 0.95 & n_unique_vals < 50){
            fcs_count = .random_noise(fcs_count)
        }

        fcs = flowCore::flowFrame(fcs_count)

        # Run workflow -------------------------------------------------

        # add ADT here?
        fcs_obj <- do.call(.setup_flowframe,
                           c(list(fcs = fcs, zero_prop = zero_prop,
                                  n_unique_vals = n_unique_vals),
                             workflow$starting_vals))

        peak_info <- .run_workflow(fcs_obj, workflow, adt_range)

        # -------------------------------------------------
        # -------------------------------------------------


            ## User defined the marker that is known to usually have multiple peaks (n = 2)
            if (adt_marker_index %in% bimodal_marker_index) {
              if (zero_prop > 0.3 && length(peak_info$peaks[, "x"]) < 2) {
                    ## less than 2 peaks and zero proportion is larger than 0.3,
                    # use finer bandwidth:fres1 instead of fres2

                    peak_info = flowStats::curvPeaks(x = fres1, dat = adt_expression,
                                                     borderQuant = 0,
                                                     from = adt_range$from,
                                                     to = adt_range$to)

                     if (length(peak_info$peaks[, "x"]) == 2) {
                         ## peak number ==2 output results.
                         y_sum = peak_info$peaks[, 'y'] %>% sum
                         res = peak_info$peaks[, "x"]

                         fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt,
                                                        bwFac = bwFac_smallest))

                         peak_info = flowStats::curvPeaks(x = fres0,
                                                  dat = adt_expression,
                                                  borderQuant = 0,
                                                  from = adt_range$from,
                                                  to = adt_range$to)

                if((length(peak_info$peaks[, "x"]) >= 2) &&
                   (sum(peak_info$peaks[, 'y']) > y_sum) &&
                   (.diff_first_two(peak_info) > 0.3)){
                  ## if using the smallest bw obtain better peak mode, switch from fres1 to fres0 results
                  res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
                }


              } else if (length(peak_info$peaks[, "x"]) > 2) {
                ## using new bandwidth, too many peaks, consider filtering out very low density peaks
                res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
              } else if (length(peak_info$peaks[, "x"]) < 2) {
                ## try with smallest bw to get more peak modes
                fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = bwFac_smallest))
                peak_info = flowStats::curvPeaks(
                  x = fres0,
                  dat = adt_expression,
                  borderQuant = 0,
                  from = adt_range$from,
                  to = adt_range$to
                )
                if (.na_or_close_peaks(peak_info, 0.5)){

                  ## smallest bw may lead to NA midpoint or peaks that are too close due to discrete value
                  peak_info = flowStats::curvPeaks(
                    x = fres1,
                    dat =  adt_expression,
                    borderQuant = 0,
                    from = adt_range$from,
                    to = adt_range$to
                  )
                  res = peak_info$peaks[, "x"]
                }else if (length(peak_info$peaks[, "x"]) >= 2) {
                  ## using new bandwidth, too many peaks, consider filtering out very low density peaks
                  res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
                }else{
                  ## still one peak left
                  res = peak_info$peaks[, "x"]
                }
              }
            } else if (zero_prop <= 0.3 && length(peak_info$peaks[, "x"]) < 2) {
              ## less than 2 peaks and small zero proportion, user finer bandwidth: fres0 instead of fres1
              fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = bwFac_smallest)) ## 1.5
              peak_info = flowStats::curvPeaks(x = fres0,
                                                dat =  adt_expression,
                                                borderQuant = 0,
                                                from = adt_range$from,
                                                to = adt_range$to)

              peak_info$peaks[, "x"] = peak_info$peaks[, "x"]
              if (..na_or_close_peaks(peak_info, 0.5)){

                peak_info = flowStats::curvPeaks(x = fres1,
                                                  dat =  adt_expression,
                                                  borderQuant = 0,
                                                  from = adt_range$from,
                                                  to = adt_range$to)

                res = peak_info$peaks[, "x"]
              }else if (length(peak_info$peaks[, "x"]) <= 2) {
                res = peak_info$peaks[, "x"]
              } else if (length(peak_info$peaks[, "x"]) > 2) {
                res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
              }
            } else {
              ## no other cases?
              res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
            }
          } else {
            ## not in user defined marker list, can have 1 peaks. Filter very low density peaks
            res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
          }
        } ## end of other marker processing


        ## all the multiple peaks are around 0
        if (length(res) > 1 && zero_prop <= 0.3 &&
            (sum(res < neg_candidate_thres) == length(res))) {
          # To do - check what the current bw is

          ## use broader bandwidth to merge multiple peaks around 0. Use fres2 instead fres1
          peak_infoTmp = flowStats::curvPeaks(
            x = fres2,
            dat =  adt_expression,
            borderQuant = border,
            from = adt_range$from,
            to = adt_range$to
          )
          # peak_infoTmp$midpoints = peak_infoTmp$peaks[, "x"]

          # Redo check as bw has changed
          if ((adt %in% bimodal_marker_index) && (length(peak_infoTmp$midpoints) == 2)) {
            ## if user define this marker to have 2 peaks.
            resTmp = peak_infoTmp$midpoints
          } else {
            # If it's supposed to have 2 peaks but doesn't, filter
            resTmp = peak_infoTmp$midpoints[peak_infoTmp$peaks[, "y"] > lower_peak_thres]
          }

          indTmp = which(!is.na(resTmp))
          resTmp = resTmp[indTmp]

          if (length(resTmp) > 1 && (sum(resTmp < 2) < length(resTmp))) {
            res = resTmp
          }
        }

        ## remove small negative peak around 0
        if (length(res) > 1 && zero_prop < 0.3 &&
            (sum(res < neg_candidate_thres) < length(res))) {
          if (peak_info$peaks[1, "x"] < 0.9 && peak_info$peaks[1, "y"] < 1 &&
              peak_info$peaks[2, "x"] > 2 &&
              peak_info$peaks[2, "y"] / peak_info$peaks[1, "y"] > 5) {
            res = res[-1]
          }
        }

        ## all the peaks around 2 and zero proportion very large.
        ## Highly likely to have only one peak.
        if (length(res) > 1 && zero_prop > 0.5 &&
            (sum(res < neg_candidate_thres) == length(res))) {
          max_peak <- which(peak_info$peaks[, "y"] == max(peak_info$peak[, "y"]))
          res = peak_info$peaks[, "x"][max_peak]
        }

        ## record peak mode
        peak_loc[[sample_name]] = res


   ## end of for loop for sample_name in sample_names

  # -----------------------------------
  # TO DO: check that positive peak is defined as this function expects
  landmark <- .adjust_peak_indices(peak_locs, positive_peak)

  ## if all the peaks are within 1 it is highly likely that there
  ## is only one negative peak
  if (max(landmark[!is.na(landmark)]) < neg_candidate_thres) {
    landmark <- .all_negative_peaks(landmark)
  }
  return(landmark)
}


# .run_curvPeaks ----
.run_curvPeaks <- function(fcs, adt, bwFac, border=border, adt_range){
  fres = flowCore::filter(fcs,
                          flowStats::curv1Filter(adt, bwFac = bwFac))
  peak_info = flowStats::curvPeaks(x=fres, dat=fcs@exprs[, adt],
                                   borderQuant=border,
                                   from=adt_range$from,
                                   to=adt_range$to)
  return(peak_info)
}

# .setup_flowframe ----
.setup_flowframe <- function(fcs, zero_prop, n_unique_vals, bwFac, border){
    # bwFac may be either a number, or a function of zero_prop
    if (is.function(bwFac)) { bwFac <- bwFac(zero_prop) }

    obj <- list(fcs = fcs,
                zero_prop = zero_prop,
                n_unique_vals = n_unique_vals,
                bwFac = bwFac,
                border = border)
    return(obj)
}

# .check_npeaks ----
# Return a function for checking the number of peaks with some args prefilled
# (e.g. `<=`, `==`)
.check_npeaks <- function(compare_f, n_expected, msg=NULL){
    check_peaks <- function(peak_info){
        result <- compare_f(length(peak_info$peaks[, "x"]), n_expected)
        if (! is.null(msg) & isFALSE(result)) { message(msg) }
        return(result)
    }
    return(check_peaks)
}

.check_npeaks_test <- function(peak_info, f, n_expected, msg=NULL, ...){
    result <- f(nrow(peak_info$peaks), n_expected)
    dot_res <- sapply(list(...), eval)
    print(c(result, dot_res))
}

.check_npeaks_2 <- function(peak_info, expr, n_expected, msg=NULL, ...){
    eval(expr)
  #dot_res <- lapply(list(...), eval)
}

## TO DO: LOG BANDWIDTH AND PEAKS BEFORE CHANGING
.log_peaks <- function(obj, peak_info){
    # Log bw, n_peaks, peak_loc, peak_sum?
}

# .update_flowframe ----
.update_flowframe <- function(obj, ...){
    args <- list(random = FALSE, bwFac = NULL, border = NULL)
    args <- modifyList(args, ...)

    if (isTRUE(args$random)){ obj <- .random_noise(obj) }
    if (! is.null(args$bwFac)) { obj$bwFac <- args$bwFac }
    #update_args <- match.call()
    #updata_args <- update_args[names(update_args) %in% c(b)]
    return(obj)
}


# .get_peaks ----
.get_peaks <- function(peak_loc, peak_type = "peaks", lower_peak_thres=0){
    if (is.null(peak_loc)){ return(NA) }

    peak_ind = peak_info[[peak_type]][, "y"] > lower_peak_thres
    res = peak_info[[peak_type]][, "x"][peak_ind]
    return(res)
}


# .diff_first_two ----
.diff_first_two <- function(peak_info){
    if (length(peak_info$peaks[, "x"]) < 2) { return(NA) }
    return(peak_info$peaks[, "x"][2] - peak_info$peaks[, "x"][1])
}

# .na_or_close_peaks ----
.na_or_close_peaks <- function(peak_info, peak_diff){
    any(is.na(peak_info$peaks[, "x"])) || .diff_first_two(peak_info) < peak_diff
}



compare_bw <- function(peak_info, obj, bw2, compare_func){
    # Compare bw currently in obj with second bw according to compare func


}


# Workflows --------------------------------------------------------------

.run_workflow <- function(fcs_obj, workflow, adt_range, ...){
    for (i in seq_along(workflow$checks)){
        print(i)
        peak_info <- .run_curvPeaks(fcs_obj$fcs, colnames(adt), fcs_obj$bwFac,
                                    fcs_obj$border, adt_range)
        result <- workflow$checks[[i]](peak_info, ...)
        if (isFALSE(result)){
          fcs_obj <- .update_flowframe(fcs_obj, workflow$check_false[[i]])
        }
        if (isTRUE(result)){ break }
    }
    return(peak_info)
}

# NOTHING DONE WITH THRESHOLD YET
trimodal_workflow <- function(bwFac_smallest, threshold = NULL){
    tri_msg <- "Didn't find 3 peaks.  Trying larger bandwidth\n"
    return(list(starting_vals = list(bwFac = bwFac_smallest, border = 0.01),
                checks = list(.check_npeaks(`==`, 3, msg = tri_msg)),
                check_false = list(list(bwFac = bwFac_smallest + 0.5))))
}

unknown_workflow <- function(bwFac_smallest){
    return(list(setup = list(bwFac = .bw_by_zero_prop),
                checks = c(.check_npeaks(`>=`, 1),
                           .check_npeaks(`>=`, 1),
                           .check_npeaks(`==`, 2)),
                # Args for .update_flowframe
                check_false = c(list(bwFac = bwFac_smallest),
                                list(random = TRUE, bwFac = 3.1)),
                check_true = c(NA, NA)))
}


bimodal_workflow <- function(){
    return(list(setup = list(bwFac = .bw_by_zero_prop),
                # unknown workflow is also run for bimodal...

                checks = c(.check_npeaks(`==`, 2), # return, no threshold
                           .check_npeaks(`>=`, 2), # return, w threshold
                           .check_npeaks(`<`, 2, expr(zero_prop > 0.3))),

                # Args for .update_flowframe
                check_false = c(NA, NA, NA),
                check_true = c(NA, NA, c(bwFac = 2, border = 0))
                ))
}


high_zero_few_peaks_workflow <- function(){
    # two peaks -> try a different bandwidth
    list(setup = list(bwFac = 2, border = 2),
         checks = c(.check_npeaks(`==`, 2), # if two peaks, change bw, update peak_info
                    .compare_peak_sep(),
                    .check_npeaks(`>`, 2), # return w threshold
                    .check_npeaks(`<`, 2), # set bw smaller
                    .na_or_close_peaks(0.5), # go back to fres1 = 3
         check_false = c(NA, NA, NA),
         check_true = c(update_peak_info, return, return_thres,
                        c(bwFac = bwFac_smallest, border = 0),
                        c(bwFac = 3, border = 0))))

}

# } else if (length(peak_info$peaks[, "x"]) > 2) {
#   ## using new bandwidth, too many peaks, consider filtering out very low density peaks
#   res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
# } else if (length(peak_info$peaks[, "x"]) < 2) {
#   ## try with smallest bw to get more peak modes
#   fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = bwFac_smallest))
#   peak_info = flowStats::curvPeaks(
#     x = fres0,
#     dat = adt_expression,
#     borderQuant = 0,
#     from = adt_range$from,
#     to = adt_range$to
#   )



# compare bwFac 2 and bwFac_smallest
# peaks with smallest >= 2 & sum of peak heights increased
# and difference >3 between second and first peak locs with bw_smallest
#  - then use smallest with threshold
# if now have too many peaks, use threshold



all_peaks_small_workflow <- function(){

}




trimodal <- trimodal_workflow(1.1)
cd4 <- trimodal_workflow(1.1, threshold = 0.001)

# Notes ----
# workflow:
# setup
# checks
# check_false
# check_true
# for seq_along checks - if false check_false[1] else check_true[1] & break

#bwFac_smallest = 1.1
#neg_candidate_thres = asinh(10/5 + 1)
#lower_peak_thres = 0.001


# For CD4 in test data - changing to bwFac_smallest + 0.5 doesn't change n_peaks
# Does change peak midpoint location
