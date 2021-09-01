#!R

args = commandArgs(trailingOnly=TRUE)

dupradar_run <- function(gtf, paired, stranded, threads) {

  ##convert stranded if needed
  if(is.numeric(stranded)){
    stranded <- as.character(stranded)
  }

  #define BAMs
  in_bams <- dir(pattern = "bam")
  dup_out <- as.list(in_bams)

  dir.create("dupradar", showWarnings = FALSE)

  for (x in 1:length(in_bams)){

  dup_out <- lapply(seq_along(in_bams), function(x){
    dupo <- dupRadar::analyzeDuprates(in_bams[x], gtf, stranded, paired, threads)
    pdf(paste0("dupradar/", in_bams[[x]], ".duprate_exp_densplot.pdf"))
      dupRadar::duprateExpDensPlot(DupMat = dupo)
    dev.off()

    pdf(paste0("dupradar/", in_bams[[x]], ".duprate_exp_boxplot.pdf"))
      dupRadar::duprateExpBoxplot(DupMat = dupo)
    dev.off()
    return(dupo)
  })
  names(dup_out) <- in_bams
  save(dup_out, file="dupradar/dupradar.analyzeDuprates.RData")
}

dupradar_run(args[1], args[2], args[3], args[4])
