
#' Title
#'
#' @param cov_granges
#' @param name_cov_grange
#' @param tn
#' @param kp
#' @param chr_select
#' @param mid.regs
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_base_loh <- function(cov_granges, name_cov_grange, tn, kp, chr_select, mid.regs, ...){

  tr.i <- 0.037
  tr.o <- 0.040
  kp <- kpDataBackground(kp, r0=tr.o*tn, r1=tr.o*tn+tr.i) %>%
    # kpAxis(ymin =0, y = 1, cex = 0.3, r0=(tr.o*tn), r1=(tr.o*tn+tr.i)) %>%
    # kpPlotRegions(data = mid.regs, r0=tr.o*tn, r1=tr.o*tn+tr.i, col = NA, lty=1, lwd=0.5, border="blue", data.panel=2) %>%
    kpHeatmap(cov_granges, y = cov_granges$mBAF, ymin=0, ymax=1, r0=tr.o*tn, r1=tr.o*tn+tr.i, colors = c("blue", "white", "red")) %>%
    kpAddLabels(labels=name_cov_grange, r0=tr.o*tn, r1=tr.o*tn+tr.i,  pos=1, label.margin = 0.035, col="red", cex=1.3) %>%
    identity()
}

#' Title
#'
#' @param raw_cov_list
#' @param file_name
#' @param chr_select
#' @param mid.regs
#' @param marker_granges
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_loh_granges_to_file <- function(raw_cov_list, file_name, chr_select = "auto", mid.regs, marker_granges, plotheight = 10, plotwidth = 50, ...) {
  tr.i <- 0.037
  tr.o <- 0.040
  num_seg_granges <- seq(0, length(raw_cov_list)-1)
  tn = max(num_seg_granges)+1
  pdf(file_name, width = plotwidth, height = plotheight)
  plot.params <- getDefaultPlotParams(plot.type=5)

  plot.params$ideogramheight <- 3
  plot.params$topmargin <- 9
  plot.params$data1height <- 1
  plot.params$data1inmargin <- 1
  plot.params$data2inmargin <- 1
  plot.params$data2height <- 150
  plot.params$data1outmargin <- 1
  plot.params$bottommargin <- 1
  plot.params$leftmargin <- 0.1
  kp <- plotKaryotype(genome="hg19", plot.type = 5, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select) %>%
    kpAddChromosomeNames(col="red", cex=1.3, srt=65)
  num_baf_granges <- seq(0, length(raw_cov_list)-1)
  for (i in names(raw_cov_list)) {
    num_baf <- which(names(raw_cov_list) == i) - 1
    plot_base_loh(raw_cov_list[[i]], i, num_baf, kp, chr_select, mid.regs)
  }
  kpPlotRegions(kp, data = marker_granges, r0=tr.o*tn, r1=tr.o*tn+tr.i, col = "black", lty=1, lwd=0.25, border="black", data.panel=2) %>%
    kpAddLabels(labels="Peak Sites", r0=tr.o*tn, r1=tr.o*tn+tr.i,  pos=1, label.margin = 0.04, col="red", cex=1.3, data.panel = 2)
  dev.off()
}

#' Plot Chromosome to File
#'
#' @param raw_cov_list
#' @param file_name
#' @param chr_select
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_chrom_to_file <- function(raw_cov_list, file_name, chr_select, ...) {
  png("Plot3.png", width = 4, height = 4, units = 'in', res = 800)
  plot.params <- getDefaultPlotParams(plot.type=1)
  plot.params <- getDefaultPlotParams(plot.type=1)
  plot.params$ideogramheight <- 5
  plot.params$data1height <- 150
  plot.params$data1inmargin <- 1
  plot.params$bottommargin <- 20
  plot.params$topmargin <- 20
  kp <- plotKaryotype(genome="hg19", plot.type = 1, plot.params = plot.params, labels.plotter = NULL, chromosomes = chr_select) %>%
    kpAddChromosomeNames(col="red", cex = 0.3) %>%
    kpAddBaseNumbers()

  map2(raw_cov_list, num_baf_granges, plot_base_seg, kp, chr_select)
  dev.off()
}


#' Title
#'
#' @param raw_cov_data
#' @param plot_title
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_bafs_stacked <- function(raw_cov_data, plot_title, ...){
  file_dest <- paste0(plot_path, plot_title, ".png")
  plot.params <- getDefaultPlotParams(plot.type=1)
  plot.params$ideogramheight <- 5
  plot.params$data1height <- 10
  plot.params$data1inmargin <- 1
  plot.params$bottommargin <- 20
  plot.params$topmargin <- 20
  kp <- plotKaryotype(genome="hg19", plot.type = 1, main = plot_title, plot.params = plot.params, labels.plotter = NULL, chromosomes = "auto") %>%
    kpDataBackground(r0=0.68, r1=0.97) %>%
    kpAddChromosomeNames(col="red", srt=30) %>%
    kpDataBackground() %>%
    kpAxis(ymin =0, y = 1, cex = 1.0) %>%
    kpHeatmap(raw_cov_data, y = raw_cov_data$mBAF, ymin=0, ymax=1)

}


#' Title
#'
#' @param raw_cov_data
#' @param plot_title
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_bafs_contiguous <- function(raw_cov_data, plot_title, ...){
  file_dest <- paste0(plot_path, plot_title, ".png")
  plot.params <- getDefaultPlotParams(plot.type=4)
  plot.params$ideogramheight <- 5
  plot.params$data1height <- 10
  plot.params$data1inmargin <- 1
  plot.params$bottommargin <- 20
  plot.params$topmargin <- 20
  kp <- plotKaryotype(genome="hg19", plot.type = 4, main = plot_title, plot.params = plot.params, labels.plotter = NULL, chromosomes = "auto") %>%
    kpDataBackground(r0=0.68, r1=0.97) %>%
    kpAddChromosomeNames(col="red", srt=30) %>%
    kpDataBackground() %>%
    kpAxis(ymin =0, y = 1, cex = 1.0) %>%
    kpHeatmap(raw_cov_data, y = raw_cov_data$mBAF, ymin=0, ymax=1)

}

#chromosome names need to be in format "chrX"
#' Title
#'
#' @param baf_segment
#'
#' @return
#' @export
#'
#' @examples
clean_baf_segment <- function(baf_segment){
  baf_segment <- mutate(baf_segment, Chr = paste0("chr", Chr)) %>%
    mutate(Chr = ifelse(Chr == "chr23", "chrX", ifelse(Chr == "chr24", "chrY", Chr))) %>%
    # dplyr::filter(mBAF > 0.7 & mBAF <= 1.00) %>%
    dplyr::filter(!grepl("N", Assay)) %>%
    identity()
}
