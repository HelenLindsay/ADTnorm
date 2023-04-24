#' Identify the valley outliers and impute by valley by closet neighbor
#' samples on the graph.
#'
#' This function identifies the valleys that tend to be outliers compared to
#' other valley locations and tries to find the closest samples that have similar
#' density distribution to impute the valley. If no neighbor sample is detected,
#' the valley will remain as original.
#' @param valley_location_res Matrix of valley landmark locations with rows
#' being samples and columns being the valleys.
#' @param adt_marker_select The marker whose valley need to be imputed.
#' Find the neighbor samples whose density distribution is close to the target
#' sample of the same ADT marker.
#' @param cell_x_adt Matrix of ADT raw counts in cells (rows) by ADT markers
#' (columns) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such
#' as cell type, sample, and batch related information.
#' @param scale Scale level to defining outlier. Larger scale value corresponds
#' to more severe ourliers.
#' @param method Outlier detection methods, choose from "MAD"
#' (Median Absolute Deviation) or "IQR" (InterQuartile Range).
#' @param nearest_neighbor_n Number of top nearest neighbor samples to detect.
#' @param nearest_neighbor_threshold Threshold to call neighbor samples.
#' @examples
#' \dontrun{
#' detect_impute_outlier_valley(valley_location_res, cell_x_feature)
#' }
# require(EMDomics)
# require(dplyr)
#' @export
detect_impute_outlier_valley <- function(valley_location_res, adt_marker_select,
                                         cell_x_adt, cell_x_feature, scale=3,
                                         method="MAD", nearest_neighbor_n=3,
                                         nearest_neighbor_threshold=0.75){

    method = match.arg(method, choices = c("MAD", "IQR"))

    # Get batch information
    valley_df <- valley_location_res %>%
        data.frame %>%
        dplyr::mutate(sample = rownames(valley_location_res)) %>%
        dplyr::left_join(valley_df, cell_x_feature %>%
                           dplyr::select(sample, batch) %>%
                           unique,
                         by = "sample")

    ## within each batch find the valley outlier and impute by the nearest
    ## neighbor samples' valley
    for (batch_each in unique(cell_x_feature$batch)){

        ## get the sample id within each batch
        sample_select <- which(valley_df$batch == batch_each)

        if (length(sample_select) > 2){ ## more than two samples per batch
            ## for each valley
            for (c in seq_along(valley_location_res)){

                ## choose outlier detection method
                valley_loc <- valley_location_res[sample_select, c]
                if (method == "MAD"){
                    row_index <- sample_select[.outlier_mad(valley_loc, scale)]
                } else if (method == "IQR"){
                    row_index <- sample_select[.outlier_mad(valley_loc, scale)]
                }

                ## for each detected outlier sample, find the nearest neighbors
                if(length(row_index) > 0){
                    message(sprintf("Outlier valley for sample: %s Valley: %s\n",
                                 valley_df$sample[row_index], c))

                  for (target_sample in valley_df$sample[row_index]){
                        ## neighbors at most 3 and earth mover distance <= 0.75 by default
                        target_neighbors <- get_neighbors(target_sample,
                                  adt_marker_select,
                                  cell_x_adt,
                                  cell_x_feature,
                                  nearest_neighbor_n = nearest_neighbor_n,
                                  nearest_neighbor_threshold = nearest_neighbor_threshold)

                        msg <- "Outlier sample %s nearest neighbors: %s\n"
                        message(sprintf(msg, target_sample, target_neighbors))

                        ## if there is qualified neighbors to impute
                        ## otherwise, this is a unique sample marker distribution.
                        ## Leave original valley value.
                        if(length(target_neighbors) > 0){
                            valley_location_res[target_sample, c] <-
                              valley_location_res[target_neighbors, c] %>%
                              stats::median
                        }
                    }
                }
            }
        }

      }
    return(valley_location_res)
}

.outlier_mad <- function(valley_loc, scale){
    abs_dev <- abs(valley_loc - stats::median(valley_loc, na.rm=TRUE))
    outlier <- abs_dev > stats::mad(valley_loc, na.rm=TRUE) * scale
    return(which(outlier))
}

.outlier_iqr <- function(valley_loc, scale){
    outlier <- stats::quantile(valley_loc, 0.75, na.rm=TRUE) +
                  scale * stats::IQR(valley_loc, na.rm=TRUE) < valley_loc
    return(which(outlier))
}
