#' Plot the expression density profile for ONE ADT marker
#'
#' This function plots adt expression density profile for only one ADT marker.
#' Each track is a sample. Color by batch
#' @param adt_count Matrix of ADT raw counts in cells (rows) by one target ADT
#' marker (column) format.
#' @param cell_x_feature Matrix of cells (rows) by cell features (columns) such
#' as cell type, sample, and batch related information.
#' @param brewer_palettes Set the color scheme of color brewer.
#' @param parameter_list Users can specify: "run_label" to give name for this
#' run; "bw" to adjust the band width of the density plot.
#' @export
#' @examples
#' \dontrun{
#' plot_adt_density_each(
#'  cell_x_adt,
#'  cell_x_feature,
#'  brewer_palettes = "Set1",
#'  parameter_list = list(bw = 0.1, run_label = "ADTnorm")
#' )
#' }
# require(ggplot2)
# require(RColorBrewer)
# require(tidyr)
# require(ggridges)
# require(ggpubr)
plot_adt_density_each = function(adt_count, cell_x_feature, brewer_palettes,
                                 parameter_list = NULL) {
    if (is.null(parameter_list)) {
        return("parameter_list is NULL!")
    }
    parameter_list_name = names(parameter_list)

    run_label = ""
    bw = 1
    if (!is.null(parameter_list)) {
        if ("run_label" %in% parameter_list_name) {
            run_label = parameter_list[["run_label"]]
        }
        if("bw" %in% parameter_list_name){
            bw = parameter_list$bw
        }
    }

    # If there is no batch, add a dummy variable
    if (! "batch" %in% colnames(cell_x_feature)){ cell_x_feature$batch <- 1 }

    tmpProfile = data.frame(counts = adt_count) %>%
        mutate(
            sample = rep(cell_x_feature$sample, 1),
            batch = rep(cell_x_feature$batch, 1)
        )

    fillColor = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, brewer_palettes))(length(unique(tmpProfile$batch)))
    resPlot = ggplot(tmpProfile, aes(x = counts, y = sample)) +
        ggridges::geom_density_ridges(aes(fill = factor(batch)), bandwidth = bw) +
        theme_bw(base_size = 20) +
        xlab(run_label) +
        ylab("") +
        scale_fill_manual(values = fillColor) +
        theme(axis.text.x = element_text(angle = 90)) +
        guides(fill = "none")

    return(resPlot)
}
