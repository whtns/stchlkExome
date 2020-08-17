
#' Sort Tumor before Cell Line
#'
#' @param karyo_granges
#'
#' @return
#' @export
#'
#' @examples
sort_T_b4_CL <- function(karyo_granges){

  karyo_names <- split(names(karyo_granges), gsub("\\..*", "", names(karyo_granges)))
  karyo_names <- lapply(karyo_names, sort, decreasing = TRUE)

  karyo_granges <- karyo_granges[unlist(karyo_names)]
}

#' Load Segmentation Objects
#'
#' @param segmentation_files
#'
#' @return
#' @export
#'
#' @examples
load_seg_objs <- function(segmentation_files){
  if (length(segmentation_files) == 1){
    segmentation_objects <- get(load(segmentation_files))
    segmentation_objects <- split(segmentation_objects$output, segmentation_objects$output$ID)
    return(segmentation_objects)
  } else{   # if loading rdata files containing single sample per file ---------------
    segment_names <- strsplit(segmentation_files,'/')
    segment_names <- sapply(segment_names, '[[', 8)

    segmentation_objects <- lapply(segmentation_files, function(x) get(load(x)))
    names(segmentation_objects) <- c(segment_names)
  }
}

#' Load Bin Objects
#'
#' @param segmentation_files
#'
#' @return
#' @export
#'
#' @examples
load_bin_objs <- function(segmentation_files){
  if (length(segmentation_files) == 1){
    bin_objects <- get(load(segmentation_files))
    bod <- bin_objects$data
    bod <- bod[, !grepl("none|N_1.*N_1", colnames(bod))]
    bod <- bod %>%
      dplyr::mutate(n_pos = rowSums(. > 0)) %>%
      dplyr::mutate(n_neg = rowSums(. < 0)) %>%
      dplyr::select(chrom, maploc, n_pos, n_neg) %>%
      identity()
    return(bod)
  } else{   # if loading rdata files containing single sample per file ---------------
    segment_names <- strsplit(segmentation_files,'/')
    segment_names <- sapply(segment_names, '[[', 8)

    bod <- lapply(segmentation_files, function(x) get(load(x)))
    names(bod) <- c(segment_names)
  }
}

#' Load All Objects
#'
#' @param segmentation_files
#'
#' @return
#' @export
#'
#' @examples
load_all_objs <- function(segmentation_files){
  if (length(segmentation_files) == 1){
    SCNA_objects <- get(load(segmentation_files))
    return(SCNA_objects)
  } else{   # if loading rdata files containing single sample per file ---------------
    segment_names <- strsplit(segmentation_files,'/')
    segment_names <- sapply(segment_names, '[[', 8)

    SCNA_objects <- lapply(segmentation_files, function(x) get(load(x)))
    names(SCNA_objects) <- c(segment_names)
  }
}


#' Bind Segmentation Plots
#'
#' @param SCNA_obj_list
#' @param sample_type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
RbindSegPlot <- function(SCNA_obj_list, sample_type = NULL, ...){
  SingleSegPlot <- function(seg_output, range.CNA = c(-2, 2),
                            color.palette = colorRampPalette(c("blue", "white", "red"))) {


    ## Use only tumor sample_ids
    if(sample_type == "tumor"){
      sample_ids <- unique(grep("^.*T.*N.*$", seg_output$ID, value = TRUE))
    } else if(sample_type == "cell_line"){   ## Use only cell line sample_ids
      sample_ids <- unique(grep("^.*CL.*N.*$", seg_output$ID, value = TRUE))
    } else
      # sample_ids <- c(unique(grep("^.*T.*N.*$", seg_output$ID, value = TRUE)), unique(grep("^.*CL.*N.*$", seg_output$ID, value = TRUE)))

      ## Use all sample_ids by default
      if(!exists("sample_ids")) {
        sample_ids <- unique(seg_output$ID)
      }

    ## Select sample_ids
    seg_output <- seg_output[seg_output$ID %in% sample_ids, ]

    # dna_copy_object$data <- dna_copy_object$data[, c("chrom", "maploc", sample_ids)]
    # names(dna_copy_object$data) <- c("chrom", "maploc", "seg.mean")

    order.sample_ids <- unique(seg_output$ID)
    seg_output$ID <- factor(seg_output$ID, levels = order.sample_ids)
    # return(list(seg_output, dna_copy_object$data))
    return(seg_output)

  }

  SCNA_obj_list <- lapply(SCNA_obj_list, SingleSegPlot)
  # seg_obj_list <- lapply(SCNA_obj_list, "[[", 1)
  # point_obj_list <- lapply(SCNA_obj_list, "[[", 2)
  seg_obj_df <- dplyr::bind_rows(SCNA_obj_list, .id = "sample_id")
  # point_obj_df <- rbindlist(point_obj_list, idcol = "sample_id")
  # return(list("segments" = seg_obj_df, "bins" = point_obj_df))
  return(list("segments" = seg_obj_df))
}

#' Create Segmentation Granges
#'
#' @param karyo_seg
#'
#' @return
#' @export
#'
#' @examples
create_seg_granges <- function(karyo_seg){

  chroms <- unique(karyo_seg$chrom)

  # prep data for karyoplots ----------------------------------------

  # rescale max gain of each sample to shrink small amp areas and increase dynamic range of gains
  ceiling_gain <- 3
  total_max_gain <- max(karyo_seg$seg.mean)
  avg_max_gain <- group_by(karyo_seg, sample_id) %>%
    summarize(top = max(seg.mean))
  median_max_gain <- median(avg_max_gain$top)

  floor_loss <- -5
  total_min_loss <- min(karyo_seg$seg.mean)
  avg_min_loss <- group_by(karyo_seg, sample_id) %>%
    summarize(top = min(seg.mean))
  median_min_loss <- median(avg_min_loss$top)

  karyo_seg_gains <- karyo_seg %>%
    dplyr::filter(seg.mean > 0) %>%
    mutate(seg.mean.scaled = ifelse(seg.mean > ceiling_gain, ceiling_gain, seg.mean)) %>%
    mutate(seg.mean.scaled = rescale(seg.mean.scaled, to=c(0,1))) %>%
    identity()

  karyo_seg_losses <- karyo_seg %>%
    dplyr::filter(seg.mean <= 0) %>%
    mutate(seg.mean.scaled = ifelse(seg.mean < floor_loss, floor_loss, seg.mean)) %>%
    mutate(seg.mean.scaled = rescale(seg.mean.scaled, to=c(-1,0))) %>%
    identity()

  karyo_seg <- bind_rows(karyo_seg_gains, karyo_seg_losses)

  # karyo_seg <- mutate(karyo_seg, seg.mean = ifelse(seg.mean > avg_max_gain, avg_max_gain, seg.mean)) %>%
  #   mutate(seg.mean = rescale(seg.mean, to=c(-1,1))) %>%  # rescale gains and losses to unit scale
  #   # mutate(seg.mean = ifelse(seg.mean < 0, rescale(seg.mean, to=c(-1,0)), seg.mean)) %>%  # rescale gains and losses to unit scale
  #   # mutate(seg.mean = ifelse(seg.mean > 0, rescale(seg.mean, to=c(0,1)), seg.mean)) %>%  # rescale gains and losses to unit scale
  #   identity()

  seg_obj_list <- split(karyo_seg, karyo_seg$sample_id)
  seg_granges <- map(seg_obj_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)

  seg_names <- split(names(seg_granges), gsub("-.*", "", names(seg_granges)))
  seg_names <- lapply(seg_names, sort, decreasing = TRUE)
  seg_granges <- seg_granges[unlist(seg_names)]

  # seg_granges <- seg_granges[!grepl("N", names(seg_granges))]


}

#' Plot SCNA group
#'
#' @param raw_cov_list
#' @param num_seg_granges
#' @param group_space
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_scna_group <- function(raw_cov_list, num_seg_granges, group_space = 0.1, kp, tr.i, tr.o, ...){
  # num_seg_granges[1] <- num_seg_granges[1]+group_space

  num_seg_granges <- num_seg_granges + c(1, rep(-1, length(num_seg_granges)-1))*group_space*seq(length(num_seg_granges))

  # num_seg_granges <- num_seg_granges + rep(c(1,-1), length(num_seg_granges))[1:length(num_seg_granges)]*group_space

  # kpDataBackground(kp, r0=tr.o*num_seg_granges[1], r1=tr.o*num_seg_granges[1]+tr.i, data.panel = 2) %>%
  #   kpAbline(h=0, r0=tr.o*num_seg_granges[1], r1=tr.o*num_seg_granges[1]+tr.i, data.panel = 2)
  map2(raw_cov_list, num_seg_granges, plot_base_seg, kp = kp, tr.i = tr.i, tr.o = tr.o, ...)
}


#' Plot Base Segmentation
#'
#' @param cov_granges
#' @param tn
#' @param kp
#' @param chr_select
#' @param cex
#' @param tr.i
#' @param tr.o
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_base_seg <- function(cov_granges, tn, kp, chr_select, cex, tr.i, tr.o, ...){

  # kp <- plotKaryotype(genome="hg19", plot.type = 1, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select) %>%
  #   kpAddChromosomeNames(col="red", cex = 0.3) %>%
  #   kpAddBaseNumbers()

  # print(unique(cov_granges$sample_id))
  kp <- kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i, data.panel = 2, color = "white") %>%
    # kpAxis(ymin =0, y = 1, cex = 0.3, r0=(tr.o*tn), r1=(tr.o*tn+tr.i)) %>%
    kpHeatmap(cov_granges, y = cov_granges$seg.mean.scaled, ymin=-1, ymax=1, r0=tr.o*tn, r1=tr.o*tn+tr.i, colors = c("blue", "white", "red"), data.panel = 2) %>%
    kpAddLabels(labels=unique(cov_granges$sample_id), r0=tr.o*tn, r1=tr.o*tn+tr.i,  pos=1, label.margin = 0.04, col="black", cex=cex, data.panel = 2) %>%
    identity()

}


#' Plot All Segmentation to Heatmap
#'
#' @param raw_cov_list
#' @param chr_select
#' @param marker_granges
#' @param marker_labels
#' @param cex
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_scna_granges <- function(raw_cov_list, chr_select = "auto", marker_granges = NULL, marker_labels = "Peak Sites", cex = 0.6, group = FALSE, group_space = 0.1, ...) {

  sample_groups <-
    names(raw_cov_list) %>%
    stringr::str_extract("^[0-9]*")

  num_seg_granges <- seq(0, length(raw_cov_list)-1)

  tn = max(num_seg_granges)+1
  tr.o = 0.99/tn
  tr.i = tr.o - 0.002
  # cex = tr.o*cex_scale
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- 3
  plot.params$topmargin <- 10
  plot.params$data1height <- 1
  plot.params$data1inmargin <- 1
  plot.params$data2inmargin <- 1
  plot.params$data2height <- 150
  plot.params$data1outmargin <- 1
  plot.params$bottommargin <- 0.1
  plot.params$leftmargin <- 0.1

  dots = list(...)
  plot.params[names(dots)] <- dots

    kp <- plotKaryotype(genome="hg19", plot.type = 3, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select, ideogram.plotter = NULL, ...) %>%
      kpAddChromosomeNames(col="black", cex = cex)
    if (group){
      raw_cov_list <- split(raw_cov_list, sample_groups)
      num_seg_granges <- split(num_seg_granges, sample_groups)
      purrr::map2(raw_cov_list, num_seg_granges, ~plot_scna_group(.x, .y, group_space = group_space, kp = kp, chr_select = chr_select, cex = cex, tr.i = tr.i, tr.o = tr.o, ...))
    } else {
      map2(raw_cov_list, num_seg_granges, plot_base_seg, kp = kp, cex = cex, tr.i = tr.i, tr.o = tr.o, ...)
    }

    # map2(raw_cov_list, num_seg_granges, plot_base_seg, kp = kp, chr_select = chr_select, cex = cex, tr.i = tr.i, tr.o = tr.o, ...)
    if(!is.null(marker_granges)){
      kpPlotRegions(kp, data = marker_granges, r0=tr.o*tn, r1=tr.o*tn+tr.i, col = "black", lty=1, lwd=0.25, border="black", data.panel=2) %>%
      kpAddLabels(labels=marker_labels, r0=tr.o*tn, r1=tr.o*tn+tr.i,  pos=1, label.margin = 0.04, col="black", cex=cex, data.panel = 2)
    }
}

#' Load Segmentation Files
#'
#' @param seg_file
#'
#' @return
#' @export
#'
#' @examples
load_seg_files <- function(seg_file){
  fs::path_filter(seg_file, "*segment.Rdata$")
}

#' Title
#'
#' @param baf_segment
#' @param chr
#' @param feature_start
#' @param feature_end
#'
#' @return
#' @export
#'
#' @examples
seg_span_gene <- function(baf_segment, chr, feature_start, feature_end){
  baf_segment <- dplyr::filter(baf_segment, chrom == chr & (!(start < feature_start & end < feature_end) & !(start > feature_end & end > feature_end)))
  # baf_segment <- baf_segment[(baf_segment$chrom == chr & (!(baf_segment$start < feature_start & baf_segment$end < feature_end) & !(baf_segment$start > feature_end & baf_segment$end > feature_end)))]
}

#' Title
#'
#' @param segmentation_files
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
clean_bin_objs <- function(segmentation_files, ...){
  bin_objects <- lapply(segmentation_files, load_bin_objs)
  bin_objects <- reduce(bin_objects, left_join, by = c("chrom", "maploc")) %>%
    dplyr::mutate(n_pos = n_pos.x + n_pos.y) %>%
    dplyr::mutate(n_neg = n_neg.x + n_neg.y) %>%
    dplyr::mutate(gain_loss_diff = n_pos - n_neg) %>%
    # dplyr::select(chrom, maploc, gain_loss_diff) %>%
    identity()

  return(bin_objects)
}

#' Title
#'
#' @param bin_objects
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
plot_bin_objects <- function(bin_objects, filename){
  ggplot(bin_objects, aes(x=maploc, y=gain_loss_diff, color=chrom)) +
    geom_line()
}

#' Title
#'
#' @param dna_copy_object
#' @param samples
#' @param remove_pattern
#' @param sample_type
#'
#' @return
#' @export
#'
#' @examples
format_SCNA_data <- function(dna_copy_object, samples, remove_pattern, sample_type = NULL) {

  if(sample_type=="proband"){
    samples <- unique(grep("[[:digit:]]_1", dna_copy_object$output$ID, value = TRUE ))
  } else if (sample_type == "all"){
    samples <- unique(dna_copy_object$output$ID)
  }

  # none_refs <- unique(grep("none", samples, value = TRUE))
  # normal_normal <- unique(grep("N_1.*N_1", samples, value = TRUE))

  remove_samples <- unlist(lapply(remove_pattern, function(x) unique(grep(x, samples, value = TRUE))))
  ## Use all samples by default
  if(missing(samples)) {
    samples <- unique(dna_copy_object$output$ID)
  }

  ## Select samples
  dna_copy_object$output <- dna_copy_object$output[!dna_copy_object$output$ID %in% remove_samples, ]

  order.samples <- unique(dna_copy_object$output$ID)
  dna_copy_object$output$ID <- factor(dna_copy_object$output$ID, levels = order.samples)
  return(dna_copy_object$output)

}

#' Title
#'
#' @param my_data
#' @param outfile
#'
#' @return
#' @export
#'
#' @examples
find_seg_genes <- function(my_data, outfile = NULL) {
  my_data <- my_data %>%
    dplyr::rename(sample = ID, start = loc.start, end = loc.end, seg.mean = seg.mean, chrom = chrom)

  peak <- my_data %>%
    dplyr::filter(!is.na(seg.mean))

  peak <- peak %>%
    as.data.frame %>%
    dplyr::select(chrom, start, end, sample, seg.mean) %>%
    mutate(chrom=paste0('chr', chrom)) %>%
    mutate(chrom=case_when(
      chrom == "chr23" ~ "chrX",
      chrom == "chr24" ~ "chrY",
      TRUE ~ as.character(chrom)
    )) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    identity()

  peak2 <- mergeByOverlaps(genes(txdb), peak)

  gene_overlaps <- data.frame(peak2$`genes(txdb)`)[,c("seqnames", "start", "end")]
  peaks <- data.frame(peak2$peak)[,c("start", "end")]
  colnames(peaks) <- c("seg.start", "seg.end")

  peaks <- cbind(gene_overlaps, peaks, "gene_id" = peak2$gene_id, "sample" = peak2$sample, "seg.mean" = peak2$seg.mean)

  entrez_to_symbol <- as.data.frame(org.Hs.egSYMBOL)

  peaks <- dplyr::right_join(entrez_to_symbol, peaks, by = "gene_id")

  # write.table(peaks, outfile, sep = ",", row.names = FALSE)
  return(peaks)

}

#' Title
#'
#' @param peak_genes_dat
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
clean_gene_objs <- function(peak_genes_dat, ...){

  gene_objects <- group_by(peak_genes_dat, gene_id, seqnames, symbol, start, end) %>%
    dplyr::count(sign(seg.mean)) %>%
    tidyr::spread(`sign(seg.mean)`, n) %>%
    dplyr::rename(n_neg = `-1`, n_pos = `1`) %>%
    dplyr::mutate_all(funs(replace(., is.na(.), 0))) %>%
    dplyr::mutate(gain_loss_diff = n_pos - n_neg) %>%
    dplyr::select(gene_id, symbol, seqnames, start, end, gain_loss_diff) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::filter(seqnames %in% c("chr1", "chr2", "chr6", "chr7", "chr13", "chr16", "chr19")) %>%
    identity()

  return(gene_objects)
}

#' Title
#'
#' @param gene_objects
#' @param title
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
plot_gene_objects <- function(gene_objects, title, filename){

  gene_objects <- clean_gene_objs(gene_objects)
  gene_plot <- ggplot(gene_objects, aes(x=start, y=gain_loss_diff, color=seqnames)) +
    geom_line() +
    ggtitle(title)

  print(gene_plot)
  return(gene_plot)
  #
  # ggplot(gene_objects, aes(x=start, y=gain_loss_diff, color=seqnames)) +
  #   geom_segment(aes(x=start, xend=end, y=gain_loss_diff, yend=gain_loss_diff)) +
  #   ggtitle(title)

}

#' Title
#'
#' @param segmentation_files
#' @param cohort
#'
#' @return
#' @export
#'
#' @examples
create_scna_df <- function(segmentation_files, cohort){
  segmentation_objects <- do.call(c, lapply(segmentation_files, load_seg_objs))

  # tidy segmentation data ---------------------------------------------------------------
  segment_df <- RbindSegPlot(segmentation_objects, sample_type = "all")

  karyo_seg <- segment_df$segments %>%
    dplyr::rename("start" = "loc.start", "end" = "loc.end") %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    mutate(chrom = gsub("chr23", "chrX", chrom)) %>%
    mutate(chrom = gsub("chr24", "chrY", chrom)) %>%
    mutate(sample_id = gsub("log2.", "", gsub("_.*", "", ID))) %>%
    identity()

  if (cohort == "vc"){
    # karyo_seg <- karyo_seg[!grepl("none", karyo_seg$ID),]
    # karyo_seg <- karyo_seg[!grepl("N", karyo_seg$sample_id),]
    karyo_seg
  }

  if(cohort == "reynolds") karyo_seg <- karyo_seg[grepl("none", karyo_seg$ID),]

  return(karyo_seg)

}

#' Title
#'
#' @param cohort
#' @param segmentation_files
#' @param as_grange
#'
#' @return
#' @export
#'
#' @examples
collate_scna_segments <- function(cohort, segmentation_files, as_grange = TRUE) {
  # browser()

  segmentation_files = segmentation_files[[cohort]]

  segmentation_files <- purrr::map(segmentation_files, load_seg_files)
  karyo_seg <- create_scna_df(segmentation_files, cohort = cohort)

  if (!as_grange) return(karyo_seg)

  seg_granges <- create_seg_granges(karyo_seg)

  if (cohort == "vc"){
    seg_granges <- sort_T_b4_CL(seg_granges)
  }

  return(seg_granges)


}


