# plot_adt_density ----
#' Plot the expression density profile with identified peak and valley locations
#'
#' This function plots adt expression density profile with identifies peak and
#' valley locations. Each track is a sample, colored by batch.
#' @param cell_x_count Matrix of ADT raw counts in cells (rows) by one target ADT
#' marker (column) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such
#' as cell type, sample, and batch related information.
#' @param adt_marker_select The target ADT marker(s) that the density plot is
#' about. Leaving it NULL will generate figures for all the ADT markers available
#' in cell_x_adt. Default: NULL
#' @param peak_landmarks Matrix of peak landmark locations with rows being
#' samples and columns being the peaks.  Default: NULL
#' @param valley_landmarks Matrix of valley landmark locations with rows
#' being samples and columns being the valleys.  Default: NULL
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
#' plot_adt_density(
#'   cell_x_adt,
#'   cell_x_feature,
#'   peak_landmark_list = peak_mode_norm_res,
#'   valley_landmark_list = valley_location_norm_res,
#'   adt_marker_select = c("CD3", "CD4", "CD8", "CD19"),
#'   brewer_palettes = "Set1",
#'   parameter_list = list(bw = 0.1, run_label = "ADTnorm")
#' )
#' }
plot_adt_density = function(adt_count, cell_x_feature,
                            peak_landmark_list = NULL,
                            valley_landmark_list = NULL,
                            brewer_palettes = "Set1", parameter_list = NULL) {

    # Set run_label and bw (bandwidth) if specified ----
    run_label = ""
    bw = 1
    if (is.null(parameter_list)) { parameter_list <- list() }

    if ("run_label" %in% names(parameter_list)) {
      run_label = parameter_list$run_label
    }
    if("bw" %in% parameter_list_name){ bw = parameter_list$bw }

    # If there is no batch, add a dummy variable ----
    if (! "batch" %in% colnames(cell_x_feature)){ cell_x_feature$batch <- 1 }

    tmpProfile = .formatCounts(cell_x_adt, cell_x_feature)
    fillColor = .make_fill_color(brewer_palettes,
                                 length(unique(tmpProfile$batch)))

    p = .adt_density_plot(tmpProfile, bw, fillColor)

    # If there are peaks and / or valleys, add landmarks to the plot
    if (! is.null(peak_landmarks)){ p <- .add_landmarks(p, peak_landmarks) }
    if (! is.null(valley_landmarks)){ p <- .add_landmarks(p, valley_landmarks) }

    return(p)
}

# .formatCounts -----
.formatCounts <- function(cell_x_adt, cell_x_feature){
    tmpProfile = data.frame(counts = cell_x_adt) %>%
      dplyr::mutate(sample = cell_x_feature$sample,
                    batch = cell_x_feature$batch)
    return(tmpProfile)
}


# .adt_density_plot -----
#' Plot ADT density per sample
#'
#' Creates a density ridge plot of ADT per sample, where samples are rows and
#' ADTs are facets.  Also adds theme elements and fill colors.
.adt_density_plot <- function(tmpProfile, bw, fillColor, base_size = 20){
    resPlot = ggplot(tmpProfile, aes(x = counts, y = sample)) +
      ggridges::geom_density_ridges(aes(fill = factor(batch)), bandwidth = bw) +
      theme_bw(base_size = base_size) +
      xlab(run_label) +
      ylab("") +
      scale_fill_manual(values = fillColor) +
      theme(axis.text.x = element_text(angle = 90)) +
      guides(fill = "none")

    return(resPlot)
}

# .add_landmarks -----
#'
#'Add vertical lines indicating peak/valley locations
#'
#'Given a matrix of landmark locations and a ggplot object p, add vertical geom
#'segments corresponding to peak locations
#'
#'@param p A ggplot object
#'@param locs A matrix of landmark locations, e.g. as returned by ADTnorm
.add_landmarks <- function(p, locs, lwd = 1){
    loc_df <- .make_loc_df(as.data.frame(locs))

    p <- p +
      geom_segment(data = peak_location, linewidth = lwd,
                 aes(x = peakx, xend = peakx, y = peaks,
                     yend = peaky + peaks)) +
      geom_segment(data = valley_location, linewidth = lwd,
                   aes(x = peakx, xend = peakx, y = peaks,
                       yend = peaky + peaks), color = "grey")
    return(p)
}

# .make_loc_df ----
#'
#'Format a matrix of locations into x-y coords to add to an ADT density plot
#'
#'@param locs A matrix of landmark locations, e.g. as returned by ADTnorm
.make_loc_df <- function(locs){
    locs %>%
      dplyr::mutate(sample = levels(cell_x_feature$sample),
                    peaks = seq_along(levels(cell_x_feature$sample)),
                    peaky = 0.5) %>%
      tidyr::pivot_longer(cols = c(colnames(loc_df)), values_to = "peakx")
}


# .make_fill_color ----
.make_fill_color <- function(brewer_palettes, n_batch){
    brewer_pal = RColorBrewer::brewer.pal(8, brewer_palettes)
    fillColor = grDevices::colorRampPalette(brewer_pal)(n_batch)
    return(fillColor)
}
