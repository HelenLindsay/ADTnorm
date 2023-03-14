#' Arcsine transformation.
#'
#' This function arcsine transforms the input cell_x_adt matrix with co-factor
#' 5.  The definition of this function is x_new <- asinh(a + b * x) + c)
#' @param cell_x_adt Matrix where rows are cells and columns are ADT markers.
#' @param parameter_list Parameter list for a: positive double that corresponds
#' to a shift about 0; b: positive double that corresponds to a scale factor;
#' c: positive double. By default a = 1, b = 1/5 and c = 0.
#' @export
#' @examples
#' \dontrun{
#' arcsinh_transform(cell_x_adt)
#' }
arcsinh_transform = function(cell_x_adt, parameter_list = NULL){
    ## parameters
    params = list(a = 1, b = 1/5, c = 0, transformationId = "ln-transformation")

    if (!is.null(parameter_list)){
        params <- modifyList(params, parameter_list)
    }

    ## transformation
    asinhTrans = do.call(flowCore::arcsinhTransform, params)

    ## output
    out = asinhTrans(cell_x_adt)
    return(out)
}
