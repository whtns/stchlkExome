#' Make a waterfall plot
#'
#' @param my_variants
#' @param mainPalette
#' @param subset_genes
#' @param sample_order
#' @param ...
#' @param clindata
#'
#' @return
#' @export
#'
#' @examples
make_waterfall_plot <- function(my_variants, mainPalette = "ggsci::default_igv", clindata, subset_genes = NULL, sample_order = NULL, clinLegCol = 5, mainLabelSize = 12, ...){
  # in long format with colnames = c("sample", "variable", "value")

  if (!is.null(subset_genes)){
    my_variants <- dplyr::filter(my_variants, gene %in% subset_genes)
  }

  clindata <- clindata %>%
    tidyr::gather("variable", "value", -sample) %>%
    dplyr::filter(!grepl("age|author", variable)) %>%
    dplyr::filter(sample %in% unique(my_variants$sample)) %>%
    # rename(sample = Sample) %>%
    identity()

  clinvars <- unique(clindata$value)

  clinVarCol <- paletteer::paletteer_d(mainPalette)[1:length(clinvars)] %>%
    set_names(clinvars)

  clinVarOrder <- clinvars

  input_wat <- my_variants %>%
    ungroup() %>%
    dplyr::select(sample, gene,
                  variant_class = Consequence,
                  amino.acid.change = HGVSp) %>%
    dplyr::mutate(variant_class = replace_na(variant_class, "unknown")) %>%
    dplyr::mutate(sample = factor(sample, levels = sort(unique(.$sample)))) %>%
    identity()

  mutation_priority <- as.character(unique(input_wat$variant_class))
  waterfall_plot <- waterfall(input_wat,
                     fileType = "Custom",
                     variant_class_order = mutation_priority,
                     mainPalette = paletteer::paletteer_d(mainPalette)[1:length(unique(input_wat$variant_class))],
                     clinData = clindata,
                     clinLegCol = clinLegCol,
                     clinVarCol = clinVarCol,
                     clinVarOrder = clinVarOrder,
                     mainXlabel = TRUE,
                     plotMutBurden = FALSE,
                     plot_proportions = FALSE,
                     sampOrder = sample_order,
                     mainLabelSize = mainLabelSize,
                     ...
  )

  return(waterfall_plot)

}
