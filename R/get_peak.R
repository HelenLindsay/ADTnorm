#' Get the peak landmark locations
#'
#' This function detect the peak landmark locations for each sample per ADT
#' marker using either the peak mode or midpoint. Using the peak midpoint can
#' be more stable across samples and less affected by the bandwidth.
#' Using the peak mode can be more accurate in determining the peak location if
#' the band width is generally proper and the local peak density is not too
#' discrete.
#'
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers
#' (columns) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such
#' as cell type, sample, and batch related information.
#' @param adt Markers to normalize
#'
#' @param bwFac_smallest The smallest band width parameter value.
#' Recommend 1.1 for general bi-modal ADT markers except CD3, CD4 and CD8.
#' @param bimodal_marker Name(s) of the ADT markers that tend to have two
#' peaks based on researchers' prior knowledge or preliminary observation on the
#' particular data to be processed.
#' @param trimodal_marker Name(s) of the ADT markers that tend to have three
#' peaks based on researchers' prior knowledge (e.g. CD4) or
#' preliminary observation on particular data to be processed.
#' @param positive_peak A list variable containing a vector of ADT marker(s)
#' and a corresponding vector of sample name(s) in matching order to specify
#' that the uni-peak detected should be aligned to positive peaks. For example,
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
#'
#' @param cd3_index Index for CD3 marker. If the ADT marker is named "CD3",
#' the program can detect the index automatically. If users named CD3 in a
#' unique way, please provide the index for CD3 markers in the colname of cell_x_adt.
#' @param cd4_index Index for CD4 marker. If the ADT marker is named "CD4",
#' the program can detect the index automatically. If users named CD4 in a
#' unique way, please provide the index for CD4 markers in the colname of cell_x_adt.
#' @param cd8_index Index for CD8 marker. If the ADT marker is named "CD8",
#' the program can detect the index automatically. If users named CD8 in a
#' unique way, please provide the index for CD8 markers in the colname of cell_x_adt.
#'
#' @importFrom flowCore flowFrame
#' @importFrom flowStats curv1Filter
#' @export
#' @examples
#' \dontrun{
#' get_peak_midpoint(
#'   cell_x_adt = cell_x_adt,
#'   cell_x_feature = cell_x_feature,
#'   adt = "CD3",
#'   neg_candidate_thres = asinh(6/5+1)
#' )
#' }
# require(flowStats)
# require(dplyr)
get_peak = function(cell_x_adt, cell_x_feature, adt = NULL,
                    peak_type = c("mode", "midpoint"), bwFac_smallest = 1.1,
                    bimodal_marker = NULL, trimodal_marker = NULL,
                    positive_peak = NULL, neg_candidate_thres = asinh(10/5 + 1),
                    lower_peak_thres = 0.001, cd3 = NULL,
                    cd4 = NULL, cd8 = NULL) {

    ## set parameters ----
    bwFac = 1.2
    border = 0.01
    peakNr = NULL
    densities = NULL
    n = 201
    indices = FALSE

    peak_type <- match.arg(peak_type)
    sample_names = levels(cell_x_feature$sample)

    # get ADT value range with a slight extension
    range_diff <- diff(range(cell_x_adt[, adt], na.rm = TRUE))
    from = min(cell_x_adt[, adt], na.rm = TRUE) - range_diff * 0.15
    to = max(cell_x_adt[, adt], na.rm = TRUE) + range_diff * 0.15

    peak_num = 0
    peak_loc = list()
    peak_region = list()

    # working here ----

    cell_x_adt <- as.matrix(cell_x_adt)

    # get peak for each sample of this processing ADT marker
    for(sample_name in sample_names){

      fcs_count <- cell_x_adt[cell_x_feature$sample == sample_name, adt,
                              drop = FALSE]
      fcs_count <- fcs_count[! is.na(fcs_count), ]

      if (nrow(fcs_count) == 0){
          peak_loc[[sample_name]] = NA
          peak_region[[sample_name]] = NA
          next
      }

      fcs = flowCore::flowFrame(fcs_count)
      n_unique_vals = nrow(unique(fcs_count))

      if (n_unique_vals == 1){
          ## only one value for this marker
          peak_loc[[sample_name]] = NA
          peak_region[[sample_name]] = matrix(NA, ncol = 2, nrow = 1)
          message(sample_name, "-Single Value!")
          next
      }

      # get proportion of cells near zero to diagnose the neg peak enrichment
      zero_prop = sum(fcs_count < 2) / nrow(fcs_count)

      # NOT NECESSARY?
      adt_expression = cell_x_adt[cell_ind, adt] ## adt value for this marker and this sample

      # if most are around 0 and there are very few unique values:
      # add random small number
      if (zero_prop > 0.95){
          if(n_unique_vals < 50){
              fcs_count= fcs_count + stats::rnorm(nrow(fcs_count), mean = 0, sd = 0.05)
          }
      }

      # different bandwidth w.r.t the zero proportion.
      if (zero_prop > 0.5) {
          fres = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = 3.1))
      } else if (zero_prop > 0.3) {
          fres = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = 3))
      } else {
          fres = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = 2))
      }

      ######

      # is fres calculated above used??

      ########

      # Can we use the flowFrame here?

      run_curvPeaks <- function(bwFac, border = border){
          fres = flowCore::filter(fcs,
                                  flowStats::curv1Filter(adt, bwFac = bwFac))
          peak_info = flowStats:::curvPeaks(x = fres, dat = fcs,
                                            borderQuant = border, from = from,
                                            to = to)
          return(peak_info)
      }


      ## processing CD4 ----
      if (adt == cd4) {
          peak_info = run_curvPeaks(bwFac_smallest)

          if (length(peak_info$peaks[, "x"]) != 3){ ## if not obtain 3 peaks, better to use a larger bw
              peak_info = run_curvPeaks(bwFac_smallest + 0.5)
          }

          # -----------------
          # do we only need the peak info??
          # ---------

          peak_ind = peak_info$peaks[, "y"] > lower_peak_thres
          res = peak_info$peaks[, "x"][peak_ind]
          res_region = peak_info$regions[peak_ind, ]

      } else if(adt %in% trimodal_marker){ ## trimodal marker

          peak_info = run_curvPeaks(bwFac_smallest)

          if(length(peak_info$peaks[, "x"]) != 3){ ## if not obtain 3 peaks, better to use a larger bw
            peak_info = run_curvPeaks(bwFac_smallest + 0.5)
          }

          res = peak_info$peaks[, "x"] # should this be midpoint??
          res_region = peak_info$regions
      } else { ## other marker
          peak_info = flowStats:::curvPeaks(x = fres, dat = adt_expression,
                                            borderQuant = border, from = from,
                                            to = to)

          ## if no peak is detected
          if(is.na(peak_info$midpoints[1])){
            ## try using the smallest bw
            fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = bwFac_smallest))
            peak_info = flowStats:::curvPeaks(x = fres0, dat = adt_expression,                                                borderQuant = 0, from = from, to = to)

            ## if still no peak detected
            if(is.na(peak_info$midpoints[1])){
              adt_expression = adt_expression + stats::rnorm(length(adt_expression), mean = 0, sd = 0.05)
              fcs_count = adt_expression %>% t %>% t %>% as.matrix()
              colnames(fcs_count) = adt
              fcs = flowCore::flowFrame(fcs_count)
              ## update fres
              fres1 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = 2))
              fres2 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = 3))
              fres3 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = 3.1))
              fres = fres3
              peak_info = flowStats:::curvPeaks(x = fres, dat = adt_expression,
                                                borderQuant = border,
                                                from = from, to = to)
            }
      }

      # HERE -----

          ## User defined the marker that is known to usually have multiple peaks (n = 2)
          if (adt_marker_index %in% bimodal_marker_index) {
            if (length(peak_info$peaks[, "x"]) == 2) {
              ## 2 peaks, perfect!
              res = peak_info$peaks[, "x"]
              res_region = peak_info$regions
            } else if (length(peak_info$peaks[, "x"]) > 2) {
              ## more than 2 peaks, consider filtering out very low density peaks
              peak_ind = peak_info$peaks[, "y"] > lower_peak_thres
              res = peak_info$peaks[, "x"][peak_ind]
              res_region = peak_info$regions[peak_ind, ]

            } else if (zero_prop > 0.3 && length(peak_info$peaks[, "x"]) < 2) {
              ## less than 2 peaks and zero proportion is larger than 0.3, use finer bandwidth:fres1 instead of fres2
              peak_info = flowStats:::curvPeaks(x = fres1, dat = adt_expression,
                                                borderQuant = 0, from = from,
                                                to = to)

              if (length(peak_info$peaks[, "x"]) == 2) {
                ## peak number ==2 output results.
                y_sum = peak_info$peaks[, 'y'] %>% sum
                res = peak_info$peaks[, "x"]
                res_region = peak_info$regions

                fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = bwFac_smallest))

                peak_info = flowStats:::curvPeaks(x = fres0,
                                                  dat = adt_expression,
                                                  borderQuant = 0,
                                                  from = from,
                                                  to = to)

                if((length(peak_info$peaks[, "x"]) >= 2) && (sum(peak_info$peaks[, 'y']) > y_sum) && (peak_info$peaks[, "x"][2] - peak_info$peaks[, "x"][1] > 0.3)){
                  ## if using the smallest bw obtain better peak mode, switch from fres1 to fres0 results
                  res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
                  res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                }


              } else if (length(peak_info$peaks[, "x"]) > 2) {
                ## using new bandwidth, too many peaks, consider filtering out very low density peaks
                res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
                res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
              } else if (length(peak_info$peaks[, "x"]) < 2) {
                ## try with smallest bw to get more peak modes
                fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = bwFac_smallest))
                peak_info = flowStats:::curvPeaks(
                  x = fres0,
                  dat = adt_expression,
                  borderQuant = 0,
                  from = from,
                  to = to
                )
                if(any(is.na(peak_info$peaks[, "x"])) || (length(peak_info$peaks[, "x"]) >= 2 && peak_info$peaks[, "x"][2] - peak_info$peaks[, "x"][1] < 0.5)){

                  ## smallest bw may lead to NA midpoint or peaks that are too close due to discrete value
                  peak_info = flowStats:::curvPeaks(
                    x = fres1,
                    dat =  adt_expression,
                    borderQuant = 0,
                    from = from,
                    to = to
                  )
                  res = peak_info$peaks[, "x"]
                  res_region = peak_info$regions
                }else if (length(peak_info$peaks[, "x"]) >= 2) {
                  ## using new bandwidth, too many peaks, consider filtering out very low density peaks
                  res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
                  res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
                }else{
                  ## still one peak left
                  res = peak_info$peaks[, "x"]
                  res_region = peak_info$regions
                }
              }
            } else if (zero_prop <= 0.3 && length(peak_info$peaks[, "x"]) < 2) {
              ## less than 2 peaks and small zero proportion, user finer bandwidth: fres0 instead of fres1
              fres0 = flowCore::filter(fcs, flowStats::curv1Filter(adt, bwFac = bwFac_smallest)) ## 1.5
              peak_info = flowStats:::curvPeaks(x = fres0,
                                                dat =  adt_expression,
                                                borderQuant = 0,
                                                from = from,
                                                to = to)

              peak_info$peaks[, "x"] = peak_info$peaks[, "x"]
              if(any(is.na(peak_info$peaks[, "x"])) || (length(peak_info$peaks[, "x"]) >= 2 && peak_info$peaks[, "x"][2] - peak_info$peaks[, "x"][1] < 0.5)){
                peak_info = flowStats:::curvPeaks(x = fres1,
                                                  dat =  adt_expression,
                                                  borderQuant = 0,
                                                  from = from, to = to)

                res = peak_info$peaks[, "x"]
                res_region = peak_info$regions
              }else if (length(peak_info$peaks[, "x"]) <= 2) {
                res = peak_info$peaks[, "x"]
                res_region = peak_info$regions
              } else if (length(peak_info$peaks[, "x"]) > 2) {
                res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
                res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
              }
            } else {
              ## no other cases?
              res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
              res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
            }
          } else {
            ## not in user defined marker list, can have 1 peaks. Filter very low density peaks
            res = peak_info$peaks[, "x"][peak_info$peaks[, "y"] > lower_peak_thres]
            res_region = peak_info$regions[peak_info$peaks[, "y"] > lower_peak_thres, ]
          }
        } ## end of other marker processing


        ## all the multiple peaks are around 0
        if (length(res) > 1 && zero_prop <= 0.3 && (sum(res < neg_candidate_thres) == length(res))) {
          ## use broader bandwidth to merge multiple peaks around 0. Use fres2 instead fres1
          peak_infoTmp = flowStats:::curvPeaks(
            x = fres2,
            dat =  adt_expression,
            borderQuant = border,
            from = from,
            to = to
          )
          # peak_infoTmp$midpoints = peak_infoTmp$peaks[, "x"]

          if ((adt %in% bimodal_marker_index) && (length(peak_infoTmp$midpoints) == 2)) {
            ## if user define this marker to have 2 peaks.
            resTmp = peak_infoTmp$midpoints
            res_regionTmp = peak_infoTmp$regions
          } else {
            resTmp = peak_infoTmp$midpoints[peak_infoTmp$peaks[, "y"] > lower_peak_thres]
            res_regionTmp = peak_infoTmp$regions[peak_infoTmp$peaks[, "y"] > lower_peak_thres, ]
          }

          indTmp = which(!is.na(resTmp))
          resTmp = resTmp[indTmp]
          if (is.null(nrow(res_regionTmp))) {
            res_regionTmp = res_regionTmp %>%
              as.matrix() %>%
              t()
          }
          res_regionTmp = res_regionTmp[indTmp, ]
          if (length(resTmp) > 1 && (sum(resTmp < 2) < length(resTmp))) {
            res = resTmp
            res_region = res_regionTmp
          }
        }

        ## remove small negative peak around 0
        if (length(res) > 1 && zero_prop < 0.3 && (sum(res < neg_candidate_thres) < length(res))) {
          if (peak_info$peaks[1, "x"] < 0.9 && peak_info$peaks[1, "y"] < 1 && peak_info$peaks[2, "x"] > 2 && peak_info$peaks[2, "y"] / peak_info$peaks[1, "y"] > 5) {
            res = res[-1]
            res_region = res_region[-1, ]
          }
        }

        ## all the peaks around 2 and zero proportion very large. Highly likely to have only one peak.
        if (length(res) > 1 && zero_prop > 0.5 && (sum(res < neg_candidate_thres) == length(res))) {
          res = peak_info$peaks[, "x"][which(peak_info$peaks[, "y"] == max(peak_info$peak[, "y"]))]
          res_region = peak_info$regions[which(peak_info$peaks[, "y"] == max(peak_info$peak[, "y"])), ]
        }

        ## record peak mode and peak regions
        peak_num = max(peak_num, length(res))
        peak_loc[[sample_name]] = res
        peak_region[[sample_name]] = matrix(NA, ncol = 2, nrow = length(res))
        peak_region[[sample_name]][1:length(res), ] = res_region

      }

  } ## end of for loop for sample_name in sample_names

  ## initiate landmark to record peak mode location
  landmark = matrix(NA, ncol = peak_num, nrow = length(sample_names))
  landmarkRegion = list()
  for (i in 1:peak_num) {
    landmarkRegion[[i]] = matrix(NA, ncol = 2, nrow = length(sample_names))
    rownames(landmarkRegion[[i]]) = sample_names
  }
  rownames(landmark) = sample_names
  for (i in names(peak_loc)) { ## go through samples
    if (!is.na(peak_loc[[i]][1])) { ## if the peak modes are detected
      peak_locNum = length(peak_loc[[i]])

      if (peak_locNum == 1) {

        ## check if the only peak should be positive peak
        pos_marker_index = which(paste0("tmpName", positive_peak$ADT_index) == adt)
        pos_sample_index = which(positive_peak$sample == names(peak_loc)[i])

        if (length(intersect(pos_marker_index, pos_sample_index)) > 0) {
          landmark[i, min(2, peak_num)] = peak_loc[[i]]
          landmarkRegion[[min(2, peak_num)]][i, ] = peak_region[[i]]
        } else {
          landmark[i, 1] = peak_loc[[i]]
          landmarkRegion[[1]][i, ] = peak_region[[i]]
        }
      } else if (peak_locNum == 2) {
        landmark[i, c(1, max(2, peak_num))] = peak_loc[[i]]

        landmarkRegion[[1]][i, ] = peak_region[[i]][1, ]
        landmarkRegion[[max(2, peak_num)]][i, ] = peak_region[[i]][2, ]
      } else if (peak_locNum == 3) {
        landmark[i, c(1, 2, max(3, peak_num))] = peak_loc[[i]]
        landmarkRegion[[1]][i, ] = peak_region[[i]][1, ]
        landmarkRegion[[2]][i, ] = peak_region[[i]][2, ]
        landmarkRegion[[max(3, peak_num)]][i, ] = peak_region[[i]][3, ]
      } else {
        landmark[i, 1:peak_locNum] = peak_loc[[i]]
        for (k in 1:peak_locNum) {
          landmarkRegion[[k]][i, ] = peak_region[[i]][k, ]
        }
      }
    }
  }

  ## if all the peaks are within 1 - highly likely that there is only one negative peak
  if (max(landmark[!is.na(landmark)]) < neg_candidate_thres) {
    landmark_new = matrix(NA, ncol = 1, nrow = nrow(landmark))
    landmarkAllMedian = stats::median(landmark[!is.na(landmark)])
    for (i in 1:nrow(landmark)) {
      landmark_nonNA = landmark[i, !is.na(landmark[i, ])]
      if (length(landmark_nonNA) > 0) {
        landmark_new[i, ] = landmark[i, which.min(abs(landmark_nonNA - landmarkAllMedian))]
      } else {
        landmark_new[i, ] = NA
      }
    }
    landmark = landmark_new
  }
  rownames(landmark) = sample_names


  return(landmark)
}

