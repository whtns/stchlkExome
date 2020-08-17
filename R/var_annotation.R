#' #' Clean Strelka CSV
#'
#' @param annovar_csv
#'
#' @return
#' @export
#'
#' @examples
clean_strelka_csv <- function(annovar_csv){

  annovar_csv <- splitstackshape::cSplit(annovar_csv, "Otherinfo", sep="\t")
  annovar_csv <- splitstackshape::cSplit(annovar_csv, "Otherinfo_3", sep=";")
  annovar_csv <- tidyr::separate(annovar_csv, Otherinfo_5, c("TUMOR-DP", "TUMOR-FDP", "TUMOR-SDP", "TUMOR-SUBDP", "TUMOR-AU", "TUMOR-CU", "TUMOR-GU", "TUMOR-TU"), sep=":")
  annovar_csv <- tidyr::separate(annovar_csv, Otherinfo_6, c("NORMAL-DP", "NORMAL-FDP", "NORMAL-SDP", "NORMAL-SUBDP", "NORMAL-AU", "NORMAL-CU", "NORMAL-GU", "NORMAL-TU"), sep=":")
  colnames(annovar_csv)[30:42] <- c("SOMATIC", "QSS", "TQSS", "NT", "QSS_NT", "TQSS_NT", "SGT", "DP", "MQ", "MQ0", "ReadPosRankSum", "SNVSB", "SomaticEVS")
  annovar_csv[,31:42] <- data.frame(sapply(annovar_csv[,31:42], ~gsub("^.*=", "", x)))

  return(annovar_csv)
}


#' Calculate Strelka VAF
#'
#' @param nuc_counts
#'
#' @return
#' @export
#'
#' @examples
calc_strelka_VAF <- function(nuc_counts){
  a <- sort(c(unlist(nuc_counts)), decreasing = TRUE)
  a <- gsub(",.*", "", a)
  a <- as.numeric(a)
  a = a[2]/a[1]

}

#' Tidy T/N VCF from Mutect2
#'
#' @param my_vcf
#'
#' @return
#' @export
#'
#' @examples
mutect2_tn_tidy <- function(my_vcf){
  # browser()
  my_vcf@rowRanges$paramRangeID <- names(my_vcf)
  csq <- data.frame(mcols(ensemblVEP::parseCSQToGRanges(my_vcf))) %>%
    tibble::rownames_to_column("paramRangeID") %>%
    mutate(paramRangeID = gsub(".[0-9]$", "", paramRangeID)) %>%
    identity()

  evcf <- S4Vectors::expand(my_vcf)

  vcf_data <- list(
    rrs = as_tibble(SummarizedExperiment::rowRanges(evcf)),

    gtdf = as_tibble(VariantAnnotation::geno(evcf)$GT) %>%
      select(GT.TUMOR = contains("CL"), GT.TUMOR = contains("T"), GT.NORMAL = contains("N")),

    afdf = as_tibble(VariantAnnotation::geno(evcf)$AF) %>%
      select(AF.TUMOR = contains("CL"), AF.TUMOR = contains("T"), AF.NORMAL = contains("N")) %>%
      map_df(as.numeric),

    addf = as_tibble(VariantAnnotation::geno(evcf)$AD) %>%
      select(AD.TUMOR.1 = matches("*CL_1.1|*T_1.1"), AD.TUMOR.2 = matches("*CL_1.2|*T_1.2"), AD.NORMAL.1 = matches("*N_1.1"), AD.NORMAL.2 = matches("*N_1.2")) %>%
      # set_names(c("AD.NORMAL.1", "AD.TUMOR.1", "AD.NORMAL.2", "AD.TUMOR.2")) %>%
      map_df(as.numeric)
    # gnomad = as_tibble(as.numeric(info(evcf)$gnomad_AF)) %>%
    # set_names(c("gnomad.AF"))
  )

  vcf_data <- dplyr::bind_cols(vcf_data)

  vcf_data <- full_join(csq, vcf_data, by = c("paramRangeID"))

  return(vcf_data)

}

#' Tidy panel of normals VCF from Mutect2
#'
#' @param my_vcf
#'
#' @return
#' @export
#'
#' @examples
mutect2_pon_tidy <- function(my_vcf){
  # browser()
  my_vcf@rowRanges$paramRangeID <- names(my_vcf)
  csq <- data.frame(mcols(ensemblVEP::parseCSQToGRanges(my_vcf))) %>%
    tibble::rownames_to_column("paramRangeID") %>%
    mutate(paramRangeID = gsub(".[0-9]$", "", paramRangeID))

  evcf <- S4Vectors::expand(my_vcf)

  vcf_data <- list(
    rrs = as_tibble(SummarizedExperiment::rowRanges(evcf)),

    gtdf = as_tibble(VariantAnnotation::geno(evcf)$GT) %>%
      set_names(c("GT.NORMAL")),

    afdf = as_tibble(VariantAnnotation::geno(evcf)$AF) %>%
      set_names(c("AF.NORMAL")) %>%
      map_df(as.numeric),

    addf = as_tibble(VariantAnnotation::geno(evcf)$AD) %>%
      set_names(c("AD.NORMAL.1", "AD.NORMAL.2")) %>%
      map_df(as.numeric)
  )

  vcf_data <- dplyr::bind_cols(vcf_data)

  vcf_data <- full_join(csq, vcf_data, by = c("paramRangeID"))

  return(vcf_data)

}

#' Tidy Strelka snvs
#'
#' @param my_vcf
#'
#' @return
#' @export
#'
#' @examples
strelka_tidy_snvs <- function(my_vcf){
  # browser()
  my_vcf@rowRanges$paramRangeID <- names(my_vcf)
  csq <- data.frame(mcols(ensemblVEP::parseCSQToGRanges(my_vcf))) %>%
    tibble::rownames_to_column("paramRangeID") %>%
    mutate(paramRangeID = gsub(".[0-9]$", "", paramRangeID))

  evcf <- S4Vectors::expand(my_vcf)

  vcf_data <- list(
    rrs = as_tibble(SummarizedExperiment::rowRanges(evcf)),

    au = as_tibble(VariantAnnotation::geno(evcf)$AU) %>%
      set_names(levels(interaction(c("AU.NORMAL", "AU.TUMOR"), c("1", "2")))),

    cu = as_tibble(VariantAnnotation::geno(evcf)$CU) %>%
      set_names(levels(interaction(c("CU.NORMAL", "CU.TUMOR"), c("1", "2")))),

    gu = as_tibble(VariantAnnotation::geno(evcf)$GU) %>%
      set_names(levels(interaction(c("GU.NORMAL", "GU.TUMOR"), c("1", "2")))),

    tu = as_tibble(VariantAnnotation::geno(evcf)$TU) %>%
      set_names(levels(interaction(c("TU.NORMAL", "TU.TUMOR"), c("1", "2"))))

    # gnomad.AF = as_tibble(VariantAnnotation::info(evcf)$gno_af_all) %>%
    #   dplyr::select("gnomad.AF" = "value")
  )

  vcf_data <- dplyr::bind_cols(vcf_data)

  vcf_data <- full_join(csq, vcf_data, by = c("paramRangeID"))


  vcf_data[,"ALT.1"]  <- NULL

  vcf_data <- vcf_data %>%
    dplyr::mutate(AD.TUMOR.1 = case_when(
      REF == "C" ~ CU.TUMOR.1,
      REF == "G" ~ GU.TUMOR.1,
      REF == "T" ~ TU.TUMOR.1,
      REF == "A" ~ AU.TUMOR.1)) %>%
    dplyr::mutate(AD.NORMAL.1 = case_when(
      REF == "C" ~ CU.NORMAL.1,
      REF == "G" ~ GU.NORMAL.1,
      REF == "T" ~ TU.NORMAL.1,
      REF == "A" ~ AU.NORMAL.1)) %>%
    dplyr::mutate(AD.TUMOR.2 = case_when(
      ALT == "C" ~ CU.TUMOR.1,
      ALT == "G" ~ GU.TUMOR.1,
      ALT == "T" ~ TU.TUMOR.1,
      ALT == "A" ~ AU.TUMOR.1)) %>%
    dplyr::mutate(AD.NORMAL.2= case_when(
      ALT == "C" ~ CU.NORMAL.1,
      ALT == "G" ~ GU.NORMAL.1,
      ALT == "T" ~ TU.NORMAL.1,
      ALT == "A" ~ AU.NORMAL.1)) %>%
    dplyr::mutate(AF.TUMOR = AD.TUMOR.2 / (AD.TUMOR.1 + AD.TUMOR.2)) %>%
    dplyr::mutate(AF.NORMAL = AD.NORMAL.2 / (AD.NORMAL.1 + AD.NORMAL.2)) %>%
    identity()

  return(vcf_data)

}


# collate vcfs ------------------------------------------------------------

#' S3 Method for collating vcfs
#'
#' @param vcf_list
#' @param tidy_function
#'
#' @return
#' @export
#'
#' @examples
collate_vcfs <- function(vcf_list, tidy_function) {
  UseMethod("collate_vcfs")
}

#' Collate vcfs
#'
#' @param vcf_list
#' @param tidy_function
#'
#' @return
#' @export
#'
#' @examples
collate_vcfs.pon <- function(vcf_list, tidy_function){

  evcf_list <- purrr::map(vcf_list, ~tidy_function(.x))
  evcf_list <- purrr::map(evcf_list, standardize_vcf_cols)

  tidy_vcfs <- dplyr::bind_rows(evcf_list, .id = "sample")

  # rfpred ------------------------------------------------------------------

  rfp_input <- dplyr::select(data.frame(tidy_vcfs), chrom = seqnames, pos = start, ref = REF, alt = ALT) %>%
    dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
    drop_na()

  rfp0 <- rfPred::rfPred_scores(variant_list=rfp_input,
                                data="~/rb_pipeline/bin/all_chr_rfPred.txtz",
                                index="~/rb_pipeline/bin/all_chr_rfPred.txtz.tbi", all.col = TRUE)

  rfp0 <- dplyr::mutate(rfp0, chromosome = paste0("chr", chromosome))

  tidy_vcfs <- dplyr::full_join(tidy_vcfs, rfp0, by = c("seqnames" = "chromosome", "start" = "position_hg19", "REF" = "reference", "ALT" = "alteration"))

}

#' Collate vcfs
#'
#' @param vcf_list
#' @param tidy_function
#'
#' @return
#' @export
#'
#' @examples
collate_vcfs.tn <- function(vcf_list, tidy_function){

  evcf_list <- purrr::map(vcf_list, ~tidy_function(.x))

  tidy_vcfs <- dplyr::bind_rows(evcf_list, .id = "sample")

  # rfpred ------------------------------------------------------------------

  rfp_input <- dplyr::select(data.frame(tidy_vcfs), chrom = seqnames, pos = start, ref = REF, alt = ALT) %>%
    dplyr::mutate(chrom = gsub("chr", "", chrom)) %>%
    drop_na()

  rfp0 <- rfPred::rfPred_scores(variant_list=rfp_input,
                                data="~/rb_pipeline/bin/all_chr_rfPred.txtz",
                                index="~/rb_pipeline/bin/all_chr_rfPred.txtz.tbi", all.col = TRUE)

  rfp0 <- dplyr::mutate(rfp0, chromosome = paste0("chr", chromosome))

  tidy_vcfs <- dplyr::full_join(tidy_vcfs, rfp0, by = c("seqnames" = "chromosome", "start" = "position_hg19", "REF" = "reference", "ALT" = "alteration"))

}

#' Standardize Vcf column names for PON samples
#'
#' @param vcf_df
#'
#' @return
#' @export
#'
#' @examples
standardize_vcf_cols <- function(vcf_df){

  standard_col_names <- c("FILTER", "GT", "AF", "AD.1", "AD.2")
  firstcol = grep("FILTER", colnames(vcf_df))
  lastcol = grep("^AD.*2$", colnames(vcf_df))

  names(vcf_df)[c(firstcol:lastcol)] <- standard_col_names
  return(vcf_df)
}


#' Retidy vcfs
#'
#' @param my_vcfs
#' @param my_pon
#'
#' @return
#' @export
#'
#' @examples
retidy_vcfs <- function(my_vcfs, my_pon, datatype){

  if (datatype %in% c("tn","cl")){
    my_vcfs <- dplyr::arrange(my_vcfs, desc(AF.TUMOR)) %>%
      dplyr::group_by(sample, snp_id) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::filter(!ExonicFunc.refGene %in% c("synonymous_SNV", "nonframeshift_insertion", "nonframeshift_deletion")) %>%
      dplyr::filter(FILTER == "PASS" | FILTER == "alt_allele_in_normal" | FILTER == "germline_risk") %>%  #need to evaluate suitablity of this filter!
      dplyr::filter((AD.NORMAL.1 != 0 | AD.NORMAL.2 != 0)) %>%
      dplyr::mutate(path_meta_score = ifelse(SIFT_score >=0.5, 1, 0) +
                      ifelse(Polyphen2_score >=0.5, 1, 0) +
                      ifelse(MutationTaster_score >=0.5, 1, 0) +
                      ifelse(LRT_score >=0.5, 1, 0)) %>%
      dplyr::filter((ExonicFunc.refGene == "nonsynonymous_SNV" & path_meta_score > 1) | ExonicFunc.refGene != "nonsynonymous_SNV") %>%
      dplyr::filter(gnomad.AF.value < 0.01 | is.na(gnomad.AF.value)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Gene.refGene = as.character(Gene.refGene))

    variants_w_o_pon <- my_vcfs %>%
      dplyr::group_by(sample, snp_id) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::filter(AF.TUMOR > 0.1) %>%
      dplyr::filter(AD.TUMOR.2 > 10) %>%
      dplyr::filter(!grepl("rs", snp_id)) %>%
      dplyr::group_by(snp_id) %>%
      dplyr::mutate(recurrence = paste(as.character(sample), collapse = ";")) %>%
      dplyr::mutate(counts = dplyr::n())

    my_pon <- dplyr::filter(my_pon, AF > 0.35)

    variants <- variants_w_o_pon %>%
      anti_join(my_pon, by = c("seqnames", "start", "end", "REF", "ALT", "GT.TUMOR" = "GT"))

  } else if (datatype == "pon"){
    my_vcfs <- dplyr::arrange(my_vcfs, desc(AF)) %>%
      dplyr::group_by(sample, snp_id) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::filter(!ExonicFunc.refGene %in% c("synonymous_SNV", "nonframeshift_insertion", "nonframeshift_deletion")) %>%
      dplyr::filter(FILTER == "PASS" | FILTER == "alt_allele_in_normal" | FILTER == "germline_risk") %>%  #need to evaluate suitablity of this filter!
      dplyr::filter((AD.1 != 0 | AD.2 != 0)) %>%
      dplyr::mutate(path_meta_score = ifelse(SIFT_score >=0.5, 1, 0) +
                      ifelse(Polyphen2_score >=0.5, 1, 0) +
                      ifelse(MutationTaster_score >=0.5, 1, 0) +
                      ifelse(LRT_score >=0.5, 1, 0)) %>%
      dplyr::filter((ExonicFunc.refGene == "nonsynonymous_SNV" & path_meta_score > 1) | ExonicFunc.refGene != "nonsynonymous_SNV") %>%
      dplyr::filter(gnomad.AF.value < 0.01 | is.na(gnomad.AF.value)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Gene.refGene = as.character(Gene.refGene))

    variants_w_o_pon <- my_vcfs %>%
      dplyr::group_by(sample, snp_id) %>%
      dplyr::filter(row_number() == 1) %>%
      dplyr::filter(AF > 0.1) %>%
      dplyr::filter(AD.2 > 10) %>%
      dplyr::filter(!grepl("rs", snp_id)) %>%
      dplyr::group_by(snp_id) %>%
      dplyr::mutate(recurrence = paste(as.character(sample), collapse = ";")) %>%
      dplyr::mutate(counts = dplyr::n())

    my_pon <- dplyr::filter(my_pon, AF > 0.35)

    variants <- variants_w_o_pon %>%
      anti_join(my_pon, by = c("seqnames", "start", "end", "REF", "ALT", "GT"))

  } else{
    stop("'datatype' describes the method of variant calling, either tumor/normal paired (tn) or panel of normals (pon)")
  }

  genes <- variants %>%
    dplyr::group_by(snp_id) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::group_by(Gene.refGene) %>%
    dplyr::mutate(gene_counts = sum(counts)) %>%
    dplyr::mutate(gene_recurrence = paste(recurrence, collapse=";")) %>%
    dplyr::mutate(gene_recurrence = purrr::map_chr(strsplit(as.character(gene_recurrence) ,";"), ~ paste(unique(.x), collapse=";"))) %>%
    dplyr::mutate(gene_recurrence_counts = purrr::map_int(strsplit(as.character(gene_recurrence) ,";"), ~ length(unique(.x)))) %>%
    dplyr::arrange(desc(counts)) %>%
    dplyr::ungroup()


  samples <- variants %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, Gene.refGene) %>%
    dplyr::arrange(desc(sample))


  samples <- aggregate(Gene.refGene ~ sample, data = unique(samples), paste, collapse = ",")

  my_list <-  list("variants_w_o_pon" = variants_w_o_pon, "variants" = variants, "genes" = genes, "samples" = samples)

}




# save tidy datasets as csvs ----------------------------------------------

#' Save Tidy csvs
#'
#' @param df_list
#' @param base_path
#'
#' @return
#' @export
#'
#' @examples
save_tidy_csvs <- function(df_list, base_path){
  dir.create(base_path)

  subpath_names <- names(df_list)

  subpaths <- purrr::map(subpath_names, ~ paste0(base_path, .x, "_tidy_table.csv"))

  clean_list_cols <- function(dataset2) {
    dataset2 <- dplyr::mutate(dataset2, recurrence = purrr::map(recurrence, paste, collapse = "; "), recurrence = unlist(recurrence))
  }

  df_list[1:2] <- purrr::map(df_list[1:2], clean_list_cols)

  purrr::map2(df_list, subpaths, function(x,y) write.table(x, y, sep = ",", row.names = FALSE))

}

# read back in tidy csvs --------------------------------------------------

#' Read Tidy csvs
#'
#' @param basepath
#'
#' @return
#' @export
#'
#' @examples
read_tidy_csvs <- function(basepath){
  my_path <- list.files(basepath, full.names = TRUE)
  tidy_vcfs <- purrr::map(my_path, read.table, sep=",", header  = TRUE)
  tidy_vcfs_names <- gsub("_.*$", "", basename(my_path))
  tidy_vcfs <- setNames(tidy_vcfs, tidy_vcfs_names)

}


