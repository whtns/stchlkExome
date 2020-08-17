#' prep stachelek clindata
#'
#' @param clindata
#'
#' @return
#' @export
#'
#' @examples
prep_stachelek_clindata <- function(clindata){
  clindata[is.na(clindata)] <- "0"

  # clindata_cl = mutate(clindata, sample = paste0(sample, "-CL")) %>%
  #   identity()
  #
  # clindata_t = mutate(clindata, sample = paste0(sample, "-T")) %>%
  #   identity()
  #
  # clindata = dplyr::bind_rows(clindata_cl, clindata_t)

  clindata <-
    clindata %>%
    dplyr::mutate(sample_type = case_when(series == "CHLA-VC-RB" ~ "CL;T",
                                            series == "CHLA-RB" ~ "CL")) %>%
    dplyr::mutate(sample = gsub(".*RB-|CHLA-", "", sample)) %>%
    dplyr::mutate(sample_type = stringr::str_split(sample_type, ";")) %>%
    tidyr::unnest(sample_type) %>%
    dplyr::mutate(sample = paste0(sample, "-", sample_type)) %>%
    # dplyr::mutate(sample = gsub(".*RB-|CHLA-", "", sample)) %>%
    identity()

}

#' Title
#'
#' @param variants
#' @param clindata
#'
#' @return
#' @export
#'
#' @examples
extract_clindata <- function(variants, clindata) {
  clindata_cols <- c(colnames(clindata), "author")

  clindata = dplyr::ungroup(variants) %>%
    dplyr::select(one_of(clindata_cols)) %>%
    dplyr::mutate(series = dplyr::coalesce(series, author)) %>%
    dplyr::distinct()
}

#' Title
#'
#' @param variant_df
#'
#' @return
#' @export
#'
#' @examples
list_recurrent_genes <- function(variant_df){
  variant_df <-
    prior_author_variants %>%
    group_by(sample, gene) %>%
    dplyr::filter(row_number() == 1) %>%
    group_by(gene) %>%
    summarize(n = dplyr::n()) %>%
    mutate(freq = n / number_of_samples) %>%
    # filter(gene == "BCOR") %>%
    mutate(total = number_of_samples) %>%
    dplyr::arrange(desc(n)) %>%
    identity()

}

#' Run WebGestaltR
#'
#' @param gene_list
#' @param gene_list_path
#' @param enrichDatabase
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
run_webgestaltr <- function(gene_list, gene_list_path, enrichDatabase = "geneontology_Biological_Process", enrichMethod = "ORA", topThr = 25) {
  gene_list_out <- fs::path(path_dir(gene_list_path), paste0(path_file(path_ext_remove(gene_list_path)), "_results"))
  dir_create(gene_list_out)
  write_delim(gene_list, gene_list_path)

  test0 <- WebGestaltR(enrichMethod = enrichMethod,
                       enrichDatabase = enrichDatabase,
                       interestGeneFile = gene_list_path,
                       interestGeneType = "genesymbol",
                       referenceSet = "genome",
                       outputDirectory = gene_list_out,
                       sigMethod = "top",
                       topThr = topThr,
                       fdrThr = 1)
  return(test0)

}

