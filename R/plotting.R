#' Make a Waterfall Plot
#'
#' @param my_variants
#' @param clindata
#' @param gene_lab_size
#'
#' @return
#' @export
#'
#' @examples
make_waterfall_plot <- function(my_variants, clindata, gene_lab_size = 10){
  # in long format with colnames = c("sample", "variable", "value")

  clindata <- dplyr::rename(clindata, sample = Sample) %>%
    tidyr::gather("variable", "value", -sample) %>%
    filter(!grepl("Age", variable)) %>%
    # rename(sample = Sample) %>%
    identity()

  clinvars <- unique(clindata$value)

  clinVarCol <- paletteer::paletteer_d("ggsci", "default_igv")[1:length(clinvars)] %>%
    set_names(clinvars)

  clinVarOrder <- clinvars

  input_wat <- my_variants %>%
    ungroup() %>%
    dplyr::select(sample, gene = SYMBOL,
                  variant_class = Consequence,
                  amino.acid.change = HGVSp) %>%
    dplyr::mutate(variant_class = replace_na(variant_class, "unknown")) %>%
    dplyr::mutate(sample = factor(sample, levels = sort(unique(.$sample))))

  mutation_priority <- as.character(unique(input_wat$variant_class))

  grob1 <- genvisr::waterfall(input_wat,
                     fileType = "Custom",
                     variant_class_order = mutation_priority,
                     out = "grob",
                     mainPalette = paletteer::paletteer_d("ggsci", "default_igv")[1:length(unique(input_wat$variant_class))],
                     clinData = clindata,
                     clinLegCol = 5,
                     clinVarCol = clinVarCol,
                     clinVarOrder = clinVarOrder,
                     mainXlabel = TRUE,
                     main_geneLabSize = gene_lab_size,
                     mainLabelSize = 4,
                     plotMutBurden = FALSE,
                     plot_proportions = FALSE
  )

  # return(list(input_wat, clindata))

}
