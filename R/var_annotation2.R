
#' print length of table (to show how many variant present)
#'
#' @param variant_df
#' @param tname
#'
#' @return
#' @export
#'
#' @examples
print_tl <- function(variant_df, tname){
  variant_df <- dplyr::group_by(variant_df, sample, paramRangeID) %>%
    dplyr::filter(row_number() == 1)

  print(glue(tname, " ", dim(variant_df)[[1]], sep = " "))
}


#' split mutect2 variant list into sample type dataframe list
#'
#' @param m2_somatic_vars
#' @param m2_filt_variant_paths
#'
#' @return
#' @export
#'
#' @examples
split_m2_somatic_vars <- function(m2_somatic_vars, m2_filt_variant_paths = NULL){

  sample_types <- c("reynolds_CL", "vc_CL", "vc_T")
  sample_type_regexes <- c("^[[:digit:]]{3}-CL", "^[[:digit:]]{2}-CL", "[[:digit:]]-T" )

  m2_somatic_vars_list <- purrr::map(sample_type_regexes, ~ dplyr::filter(m2_somatic_vars, grepl(.x, sample)))

  # purrr::map2(m2_somatic_vars_list, m2_filt_variant_paths, write_csv, sep = "\t")

  names(m2_somatic_vars_list) <- sample_types

  return(m2_somatic_vars_list)


}

#' Tidy Strelka Indels
#'
#' @param my_vcf
#'
#' @return
#' @export
#'
#' @examples
strelka_tidy_indels <- function(my_vcf){
  # browser()
  my_vcf@rowRanges$paramRangeID <- names(my_vcf)
  csq <- data.frame(mcols(ensemblVEP::parseCSQToGRanges(my_vcf))) %>%
    tibble::rownames_to_column("paramRangeID") %>%
    mutate(paramRangeID = gsub(".[0-9]$", "", paramRangeID))

  evcf <- S4Vectors::expand(my_vcf)

  vcf_data <- list(
    rrs = as_tibble(SummarizedExperiment::rowRanges(evcf)),

    TAR = as_tibble(geno(evcf)$TAR) %>%
      set_names(paste0("TAR.", colnames(.))),

    TIR = as_tibble(geno(evcf)$TIR) %>%
      set_names(paste0("TIR.", colnames(.))),

    gnomad.AF = as_tibble(VariantAnnotation::info(evcf)$gno_af_all) %>%
      dplyr::select("gnomad.AF" = "value")

  )

  vcf_data <- dplyr::bind_cols(vcf_data)

  vcf_data <- full_join(csq, vcf_data, by = c("paramRangeID"))

  vcf_data <- vcf_data %>%
    dplyr::mutate(AD.TUMOR.2 = TIR.TUMOR.1, AD.TUMOR.1 = TAR.TUMOR.1) %>%
    dplyr::mutate(AD.NORMAL.2 = TIR.NORMAL.1, AD.NORMAL.1 = TAR.NORMAL.1) %>%
    dplyr::mutate(AF.TUMOR = AD.TUMOR.2 / (AD.TUMOR.1 + AD.TUMOR.2)) %>%
    dplyr::mutate(AF.NORMAL = AD.NORMAL.2 / (AD.NORMAL.1+  AD.NORMAL.2)) %>%
    identity()

  return(vcf_data)

}

#' Antijoin Mutect Variants with Panel of Normals
#'
#' @param my_variants
#' @param my_pon
#'
#' @return
#' @export
#'
#' @examples
pon_subtract_m2_vars <- function (my_variants, my_pon) {
  # browser()
  # my_pon <- dplyr::filter(my_pon, AF > 0.35)
  my_pon <- dplyr::filter(my_pon, AF > 0.05)

  variants <- my_variants %>%
    anti_join(my_pon, by = c("seqnames", "start", "end", "REF", "ALT", "GT.TUMOR" = "GT")) %>%
    dplyr::filter(!grepl("[a-z]", FILTER)) %>%
    dplyr::filter(!is.na(seqnames)) %>%
    identity()

  return(variants)
}

#' Antijoin Mutect2 Reynolds Samples with Panel of Normals
#'
#' @param my_variants
#' @param my_pon
#'
#' @return
#' @export
#'
#' @examples
pon_subtract_m2_reynolds <- function(my_variants, my_pon){
  my_pon <- dplyr::filter(my_pon, AF > 0.05)
  variants <- my_variants %>%
    anti_join(my_pon, by = c("seqnames", "start", "end", "REF", "ALT", "GT" = "GT")) %>%
    dplyr::filter(!grepl("[a-z]", FILTER)) %>%
    dplyr::filter(!is.na(seqnames)) %>%
    identity()
  return(variants)
}


#' Filter Variant Calls
#'
#' @param my_variants
#'
#' @return
#' @export
#'
#' @examples
filter_calls  <- function (myv) {

  AF_columns = c("AF.TUMOR", "AF")
  AD_columns = c("AD.TUMOR.2", "AD.2")

  # myv <- myv %>%
  #   dplyr::arrange(desc(one_of(AF_columns))) %>%
  #   dplyr::group_by(sample, paramRangeID) %>%
  #   dplyr::filter(row_number() == (1))

  myv <- myv %>%
    # dplyr::filter(FILTER == "PASS" ) %>% # filter out variants which don't meet variant caller filtering criteria
    # dplyr::filter(!grepl(paste(c("intron_variant", "synonymous"), collapse="|"), Consequence)) %>% # filter out intronic or synonymous
    dplyr::filter_at(vars(one_of(AF_columns)), ~ . > 0.05) %>% # filter out tumor variant allele frequency less than 5 percent
    dplyr::filter_at(vars(one_of(AD_columns)), ~ . > 5) %>% # filter out tumor alternate allele read depth less than 5
    identity()

  myv <- myv %>%
    # dplyr::filter(!grepl("rs", paramRangeID)) %>% # filter out known snps
    dplyr::ungroup() %>%
    dplyr::mutate(SYMBOL = as.character(SYMBOL)) %>%
    dplyr::group_by(paramRangeID) %>%
    dplyr::mutate(recurrence = paste(as.character(sample), collapse = ";")) %>%
    dplyr::mutate(counts = dplyr::n()) %>%
    dplyr::mutate(path_meta_score = ifelse(SIFT_score >= 0.5, 1, 0) + ifelse(Polyphen2_score >= 0.5, 1, 0) + ifelse(MutationTaster_score >= 0.5, 1, 0) + ifelse(LRT_score >= 0.5, 1, 0)) %>%
    identity()

  return(myv)
}

#' Fitler Set of Variants
#'
#' @param my_variants
#' @param datatype
#'
#' @return
#' @export
#'
#' @examples
filter_variant_set  <- function (my_variants, datatype, gnomad_threshold = 0.0005) {
  if (datatype %in% c("tn")) {
    # browser()
    myv <- filter_calls(my_variants)
    myv <- myv %>%
      dplyr::filter((AD.NORMAL.1 != 0 | AD.NORMAL.2 != 0)) %>% # filter out variants with zero normal reads (ref or alt)
      dplyr::filter(gnomAD_AF < gnomad_threshold | is.na(gnomAD_AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
      identity()

    return(myv)

  } else if (datatype == "pon") {
    # browser()
    myv <- filter_calls(my_variants)
    myv <- myv %>%
      dplyr::filter(gnomAD_AF < gnomad_threshold | is.na(gnomAD_AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
      identity()

    return(myv)
  } else if (datatype == "strelka") {
    # browser()
    myv <- filter_calls(my_variants)
    myv <- myv %>%
      dplyr::filter(gnomAD_AF < gnomad_threshold | is.na(gnomAD_AF)) %>%  # filter out variants with a gnomad frequency greater than set threshold
      identity()

    return(myv)
  } else {
    stop("'datatype' describes the method of variant calling, either tumor/normal paired (tn) or panel of normals (pon)")
  }
}

#' Split Variant Set into Tumor and Cell Line Sets
#'
#' @param mydf
#'
#' @return
#' @export
#'
#' @examples
divide_t_cl <- function(mydf) {
  mydf <- dplyr::mutate(mydf, sample.Cell.Line = ifelse(grepl("C", sample), as.character(sample),"")) %>%
    dplyr::mutate(sample.Tumor = ifelse(grepl("T", sample), as.character(sample), "")) %>%
    dplyr::mutate(AF.CELL_LINE = ifelse(grepl("C", sample), as.character(AF.TUMOR), "")) %>%
    dplyr::mutate(AF.TUMOR = ifelse(grepl("T", sample), as.character(AF.TUMOR), "")) %>%
    dplyr::mutate(AD.TUMOR.1 = as.character(AD.TUMOR.1), AD.TUMOR.2 = as.character(AD.TUMOR.2)) %>%
    dplyr::mutate(AD.CELL_LINE.1 = ifelse(grepl("CL", sample), AD.TUMOR.1, "")) %>%
    dplyr::mutate(AD.CELL_LINE.2 = ifelse(grepl("CL", sample), AD.TUMOR.2, "")) %>%
    dplyr::mutate(AD.TUMOR.1 = ifelse(grepl("T", sample), AD.TUMOR.1, "")) %>%
    dplyr::mutate(AD.TUMOR.2 = ifelse(grepl("T", sample), AD.TUMOR.2, "")) %>%
    identity()

  return(mydf)
}

#' Compare variants with four prior studies
#'
#' @param cobrinik_vars
#' @param reference_vars
#' @param reference_genes
#' @param caller
#'
#' @return
#' @export
#'
#' @examples
five_author_intersect <- function(cobrinik_vars, reference_vars, reference_genes, caller){

  reference_vars <- mutate(reference_vars, start = as.numeric(start), end = as.numeric(end))

  stchlk_five_author_var_intxn <-  inner_join(cobrinik_vars, reference_vars, by = c("start", "end", "seqnames" = "chr", "SYMBOL" = "gene"))

  reference_vars <- dplyr::select(reference_vars, gene, reference_author = author, reference_sample = sample) %>%
    dplyr::group_by(reference_author, gene) %>%
    dplyr::mutate(reference_samples = paste0(reference_sample, collapse = "; ")) %>%
    dplyr::mutate(recurrence = dplyr::n()) %>%
    dplyr::select(-reference_sample) %>%
    dplyr::group_by(gene, reference_author, reference_samples) %>%
    dplyr::distinct() %>%
    identity()

  stchlk_four_author_gene_intxn <- inner_join(cobrinik_vars, reference_vars, by = c("SYMBOL" = "gene")) %>%
    dplyr::filter(SYMBOL != "RB1") %>%
    dplyr::select(sample, SYMBOL, reference_author, reference_samples, everything()) %>%
    identity()

  stchlk_four_author_gene_intxn <- stchlk_four_author_gene_intxn %>%
    dplyr::group_by(sample, SYMBOL) %>%
    dplyr::mutate(recurrence = length(str_split(reference_samples, ";"))) %>%
    dplyr::mutate(reference_samples = paste(reference_author, ":", reference_samples)) %>%
    dplyr::mutate(reference_samples = paste0(reference_samples, collapse = "; ")) %>%
    dplyr::mutate(recurrence = sum(recurrence.y)) %>%
    dplyr::select(-reference_author) %>%
    dplyr::select(recurrence, everything()) %>%
    dplyr::select(-recurrence.y) %>%
    dplyr::distinct() %>%
    identity()

  return(list(variants = stchlk_five_author_var_intxn, genes = stchlk_four_author_gene_intxn))

}

#' Compare Variants in Matched Tumor and Cell Line Samples
#'
#' @param cobrinik_vars
#' @param caller
#'
#' @return
#' @export
#'
#' @examples
compare_t_cl <- function(cobrinik_vars, caller){
  # browser()
  # find intersection between tumor and cell line variants ------------------
  t_vars <- dplyr::filter(cobrinik_vars, grepl("T", sample)) %>%
    mutate(prefix = substr(sample, 1,2)) %>%
    identity()

  cl_vars <- dplyr::filter(cobrinik_vars, grepl("CL", sample)) %>%
    mutate(prefix = substr(sample, 1,2)) %>%
    dplyr::rename(AF.CELL_LINE = AF.TUMOR) %>%
    identity()



  t_cl_vars <- dplyr::inner_join(t_vars, cl_vars, by = c("seqnames", "start", "end", "REF", "ALT", "prefix", "SYMBOL")) %>%
    # dplyr::filter(AF.TUMOR.x > 0.279 | AF.TUMOR.y > 0.279) %>%
    arrange(SYMBOL, prefix) %>%
    ungroup() %>%
    dplyr::rename(AD.CELL_LINE.1 = AD.TUMOR.1.y, AD.CELL_LINE.2 = AD.TUMOR.2.y) %>%
    dplyr::select(seqnames, start, end, SYMBOL, HGVSc = HGVSc.x, HGVSp = HGVSp.x, FILTER.TUMOR = FILTER.x, sample.Tumor = sample.x,
                  AF.TUMOR, AD.TUMOR.1 = AD.TUMOR.1.x, AD.TUMOR.2 = AD.TUMOR.2.x, FILTER.CELL_LINE = FILTER.y, sample.Cell.Line = sample.y, AF.CELL_LINE,
                  AD.CELL_LINE.1, AD.CELL_LINE.2, Consequence = Consequence.x) %>%
    mutate(AF.TUMOR = as.numeric(AF.TUMOR), AF.CELL_LINE = as.numeric(AF.CELL_LINE)) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    identity()

  t_cl_vars_besides_RB1 <- dplyr::filter(t_cl_vars, SYMBOL != "RB1")

  write_csv(t_cl_vars, paste0("~/rb_pipeline/doc/RB_exome_manuscript/SNV/", caller, "_t_cl_vars.csv"))
  write_csv(t_cl_vars_besides_RB1, paste0("~/rb_pipeline/doc/RB_exome_manuscript/", caller, "_t_cl_vars_besides_RB1.csv"))

  theme_set(theme_bw())

  cons_plot <- ggplot(t_cl_vars, aes(factor(Consequence))) + geom_bar(stat="count") +
    theme(axis.text.x = element_text(angle = 90,  hjust = 1),
          axis.ticks = element_blank(),
          axis.text.y = element_blank())
  print(cons_plot)

  return(t_cl_vars)

}

#' Load VCFs
#'
#' @param vcf_file_paths
#' @param sample_type
#'
#' @return
#' @export
#'
#' @examples
load_vcfs <- function(vcf_file_paths, sample_type){

  vcf_sample_names <- gsub("_.*$", "", basename(vcf_file_paths))
  vcf_sample_list <- purrr::map(vcf_file_paths, ~VariantAnnotation::readVcf(.x, "hg19"))
  names(vcf_sample_list) <- vcf_sample_names
  vcf_sample_list_path  <- paste0("~/rb_pipeline/results/SNV/", sample_type, "_list.rds")
  saveRDS(vcf_sample_list, file = vcf_sample_list_path)
  print(vcf_sample_list_path)
  return(vcf_sample_list)
}

#' Filter VCFs by gold standard variants
#'
#' @param vcf
#' @param gold_vars
#'
#' @return
#' @export
#'
#' @examples
filt_vcf_by_gold <- function(vcf, gold_vars){
  filt <- FilterRules(list(isgold = ~ rownames(x) %in% gold_vars$snp_id))
  ## Apply
  filt2 <- filterVcf(m2_tumor_filenames[[1]], "hg19", tempfile(), filters=filt)
  ## Filtered results
  filt_vcf <- readVcf(filt2, "hg19")

  return(filt_vcf)
}

#' Load Filtered VCFs
#'
#' @param vcf_file_paths
#' @param gold_vars
#'
#' @return
#' @export
#'
#' @examples
load_filt_vcfs <- function(vcf_file_paths, gold_vars){
  vcf_sample_names <- gsub("_.*$", "", basename(vcf_file_paths))
  vcf_sample_list <- purrr::map(vcf_file_paths, filt_vcf_by_gold, gold_vars)
  names(vcf_sample_list) <- vcf_sample_names
  return(vcf_sample_list)
}

#' Recalculate geo from webgestalt
#'
#' @param geo_output
#' @param gene_input
#' @param intgenelength
#'
#' @return
#' @export
#'
#' @examples
recalculate_geo <- function(geo_output, gene_input, intgenelength = 245){
  # browser()
  gene_input <- gene_input %>%
    janitor::tabyl(gene) %>%
    dplyr::mutate(userId = map2(.$gene, .$n, rep)) %>%
    dplyr::mutate(userId = map_chr(.$userId, paste0, collapse = ";")) %>%
    identity()

  old_enrichment <- geo_output

  enrichment <- old_enrichment %>%
    dplyr::mutate(symbols = str_split(userId, ";")) %>%
    select(-userId) %>%
    unnest(symbols) %>%
    dplyr::left_join(gene_input, by = c("symbols" = "gene")) %>%
    dplyr::group_by(geneSet) %>%
    dplyr::mutate(userId = paste0(userId, collapse = ";")) %>%
    dplyr::mutate(added_count = n-1)

  geneSet_adds <- enrichment %>%
    dplyr::select(geneSet, added_count) %>%
    dplyr::summarise(added_count = sum(added_count))

  enrichment2 <-
    enrichment %>%
    select(-added_count, -symbols, -n, -percent) %>%
    dplyr::distinct() %>%
    dplyr::right_join(geneSet_adds, by = "geneSet") %>%
    # dplyr::filter(added_count.y == max(added_count.y)) %>%
    identity()

  # FDR size 6306
  new_enrichment <-
    enrichment2 %>%
    dplyr::mutate(overlap = overlap + added_count) %>%
    dplyr::mutate(enrichmentRatio = overlap/expect) %>%
    dplyr::mutate(pValue = 1-phyper(overlap - 1, intgenelength, 16666 - intgenelength, size, lower.tail=TRUE, log.p=FALSE)) %>%
    # dplyr::select(-n, -percent) %>%
    dplyr::mutate(FDR = ifelse(added_count == 0, FDR, p.adjust(pValue, method = "BH"))) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::slice_head(n = 20) %>%
    identity()

  enrichment_compare <- dplyr::left_join(old_enrichment, new_enrichment, by = c("geneSet", "description", "link"), suffix = c(".old", ".new")) %>%
    select(order(colnames(.))) %>%
    dplyr::filter(!is.na(enrichmentRatio.new))

  go_column_names <- c("size", "overlap", "expect", "enrichmentRatio",
                       "pValue", "FDR", "overlapId", "userId")

  go_columns = paste0(go_column_names, ".new")
  names(go_columns) <- go_column_names

  go_columns <- c("geneSet", "description", go_columns)

  test2 <- enrichment_compare %>%
    dplyr::select(go_columns)
}


#' Title
#'
#' @param vars
#' @param antiseries
#'
#' @return
#' @export
#'
#' @examples
calc_recurrent_vars <- function(vars, antiseries = "CHLA-RB") {

  recurrent_vars <-
    vars %>%
    dplyr::filter(series != antiseries) %>%
    dplyr::group_by(gene, sample) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::mutate(sample_type = ifelse(str_detect(sample, "-T"), "Tumor", "Cell Type")) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(sample_types = paste0(sample_type, collapse = "; ")) %>%
    # dplyr::mutate(sample_types = strsplit(sample_types, split = "; ")) %>%
    # dplyr::filter(grepl("Tumor", sample_types)) %>%
    dplyr::select(-sample_types) %>%
    dplyr::mutate(recurrence = dplyr::n()) %>%
    dplyr::mutate(sample_number = str_extract(sample, "[0-9]+")) %>%
    dplyr::filter(n_distinct(sample_number) > 1) %>%
    dplyr::mutate(recurrence = paste0(author, ":", sample, collapse = "; ")) %>%
    identity()

  return(recurrent_vars)
}
