#' Plot the expression density profile for ONE ADT marker with identifies peak
#' and valley locations
#'
#' This function plots adt expression density profile with identifies peak and
#' valley locations for only one ADT marker. Each track is a sample. Color by batch
#' @param adt_count Matrix of ADT raw counts in cells (rows) by one target ADT
#' marker (column) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such
#' as cell type, sample, and batch related information.
#' @param peak_landmark_list Matrix of peak landmark locations with rows being
#' samples and columns being the peaks.
#' @param valley_landmark_list Matrix of valley landmark locations with rows
#' being samples and columns being the valleys.
#' @param brewer_palettes Set the color scheme of color brewer.
#' @param parameter_list Users can specify: "run_label" to give name for this
#' run; "bw" to adjust the band width of the density plot.
#' @importFrom ggridges geom_density_ridges
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr pivot_longer
#' @export
#' @examples
#' \dontrun{
#' plot_adt_density_with_peak_valley_each(
#'   cell_x_adt,
#'   cell_x_feature,
#'   peak_landmark_list = peak_mode_norm_res,
#'   valley_landmark_list = valley_location_norm_res,
#'   brewer_palettes = "Set1",
#'   parameter_list = list(bw = 0.1, run_label = "ADTnorm")
#' )
#' }
# require(ggplot2)
# require(RColorBrewer)
# require(tidyr)
# require(ggridges)
# require(ggpubr)
plot_adt_density_with_peak_valley_each = function(adt_count, cell_x_feature,
                                                  peak_landmark_list,
                                                  valley_landmark_list,
                                                  brewer_palettes = "Set1",
                                                  parameter_list = NULL) {

    if (is.null(parameter_list)) { parameter_list <- list() }

    parameter_list_name = names(parameter_list)

    run_label = ""
    bw = 1

    if ("run_label" %in% parameter_list_name) {
        run_label = parameter_list[["run_label"]]
    }
    if("bw" %in% parameter_list_name){
        bw = parameter_list$bw
    }

    # If there is no batch, add a dummy variable
    if (! "batch" %in% colnames(cell_x_feature)){ cell_x_feature$batch <- 1 }

    tmpProfile = data.frame(counts = adt_count) %>%
        dplyr::mutate(sample = cell_x_feature$sample,
                      batch = cell_x_feature$batch)

    make_loc_df <- function(loc_df){
        loc_df %>%
            dplyr::mutate(sample = levels(cell_x_feature$sample),
                          peaks = seq_along(levels(cell_x_feature$sample)),
                          peaky = 0.5) %>%
            tidyr::pivot_longer(cols = c(colnames(loc_df)), values_to = "peakx")
    }

    peak_location = make_loc_df(as.data.frame(peak_landmark_list))
    valley_location = make_loc_df(as.data.frame(valley_landmark_list))

    n_batch = length(unique(tmpProfile$batch))
    brewer_pal = RColorBrewer::brewer.pal(8, brewer_palettes)
    fillColor = grDevices::colorRampPalette(brewer_pal)(n_batch)

    resPlot = ggplot(tmpProfile, aes(x = counts, y = sample)) +
        ggridges::geom_density_ridges(aes(fill = factor(batch)),
                                      bandwidth = bw) +
        geom_segment(data = peak_location,
                     aes(x = peakx, xend = peakx, y = peaks,
                         yend = peaky + peaks),
                     linewidth = 1) +
        geom_segment(data = valley_location,
                     aes(x = peakx, xend = peakx, y = peaks,
                         yend = peaky + peaks),
                     linewidth = 1, color = "grey") +
        theme_bw(base_size = 20) +
        xlab(run_label) +
        ylab("") +
        scale_fill_manual(values = fillColor) +
        theme(axis.text.x = element_text(angle = 90)) +
        guides(fill = "none")

    return(resPlot)
}
