#' @name sc_atac_create_qc_per_bc_file
#' @title generating a file useful for producing the qc plots
#' @description uses the peak file and annotation files for features
#' 
#' @param frags_file The fragment file
#' @param peaks_file The peak file
#' @param promoters_file The path of the promoter annotation file 
#' @param tss_file The path of the tss annotation file 
#' @param enhs_file The path of the enhs annotation file 
#' @param output_folder

#' @importFrom data.table fread setkey copy :=
#' 
#' @export
#' 
sc_atac_create_qc_per_bc_file <- function(frags_file,
                                          peaks_file,
                                          promoters_file,
                                          tss_file,
                                          enhs_file,
                                          output_folder) {
  
  
  # https://github.com/wbaopaul/scATAC-pro/blob/master/scripts/src/get_qc_per_barcode.R
  
  # chr <- .N <- bc <- total_frags <- isMito <- NULL
  # 
  # out.frag.overlap.file <- file.path(output_folder, "qc_per_bc_file.txt")
  # 
  # frags <- fread(frags_file, select=1:4, header = FALSE)
  # names(frags) <- c('chr', 'start', 'end', 'bc')
  # setkey(frags, chr, start)
  # 
  # frags[, 'total_frags' := .N, by = bc]
  # frags <- frags[total_frags > 5]
  # 
  # frags <- frags[!grepl(chr, pattern = 'random', ignore.case = TRUE)]
  # frags <- frags[!grepl(chr, pattern ='un', ignore.case = TRUE)]
  # 
  # peaks <- fread(peaks_file, select=1:3, header = FALSE)
  # promoters <- fread(promoters_file, select=1:3, header = FALSE)
  # tss <- fread(tss_file, select=1:3, header = FALSE)
  # enhs <- fread(enhs_file, select=1:3, header = FALSE)
  # names(peaks) = names(promoters) = names(tss) =
  #   names(enhs) = c('chr', 'start', 'end')
  # 
  # 
  # setkey(peaks, chr, start)
  # setkey(promoters, chr, start)
  # setkey(tss, chr, start)
  # setkey(enhs, chr, start)
  # 
  # chrs <- unique(frags$chr)
  # 
  # tss[, 'start' := start - 1000]
  # tss[, 'end' := end + 1000]
  # fragsInRegion <- NULL
  # 
  # for(chr0 in chrs){
  #   peaks0 <- peaks[chr == chr0]
  #   promoters0 <- promoters[chr == chr0]
  #   
  #   tss0 <- tss[chr == chr0]
  #   enhs0 <- enhs[chr == chr0]
  #   frags0 <- frags[chr == chr0]
  #   frags <- frags[chr != chr0]
  #   if(nrow(peaks0) == 0){
  #     frags0[, 'peaks' := 0]
  #   }else{
  #     frags0[, 'peaks' := sc_atac_getOverlaps_read2AnyRegion(frags0, peaks0)]
  #   }
  #   
  #   if(nrow(promoters0) == 0){
  #     frags0[, 'promoters' := 0]
  #   }else{
  #     frags0[, 'promoters' := sc_atac_getOverlaps_read2AnyRegion(frags0, promoters0)]
  #   }
  #   
  #   if(nrow(tss0) == 0){
  #     frags0[, 'tss' := 0]
  #   }else{
  #     frags0[, 'tss' := sc_atac_getOverlaps_read2AnyRegion(frags0, tss0)]
  #   }
  #   
  #   if(nrow(enhs0) == 0){
  #     frags0[, 'enhs' := 0]
  #   }else{
  #     frags0[, 'enhs' := sc_atac_getOverlaps_read2AnyRegion(frags0, enhs0)]
  #   }
  #   
  #   
  #   fragsInRegion <- rbind(fragsInRegion, frags0)
  #   cat(paste(chr0, 'Done! \n'))
  # }
  # rm(frags)
  # 
  # # get barcode region count matrix
  # fragsInRegion[, 'isMito' := ifelse(chr %in% c('chrM'), 1, 0)]
  # fragsInRegion[, c('chr', 'start', 'end') := NULL]
  # fragsInRegion[, 'frac_mito' := sum(isMito)/total_frags, by = bc]
  # fragsInRegion[, 'isMito' := NULL]
  # fragsInRegion[, 'frac_peak' := sum(peaks)/total_frags, by = bc]
  # fragsInRegion[, 'peaks' := NULL]
  # fragsInRegion[, 'frac_promoter' := sum(promoters)/total_frags, by = bc]
  # fragsInRegion[, 'promoters' := NULL]
  # fragsInRegion[, 'frac_tss' := sum(tss)/total_frags, by = bc]
  # fragsInRegion[, 'tss' := NULL]
  # fragsInRegion[, 'frac_enhancer' := sum(enhs)/total_frags, by = bc]
  # fragsInRegion[, 'enhs' := NULL]
  # fragsInRegion <- unique(fragsInRegion)
  # 
  # utils::write.table(fragsInRegion, file = out.frag.overlap.file, sep = '\t',
  #                    row.names = FALSE, quote = FALSE)
  
  pro.gr <- rtracklayer::import(promoters_file)
  tss_df <- data.table::fread(tss_file, select=c(1:3), header = F, col.names = c("chr", "start", "end"))
  tss.gr <- GenomicRanges::makeGRangesFromDataFrame(tss_df)
  enhs.gr <- rtracklayer::import(enhs_file)
  
  # Peak file
  peaks.gr <-  rtracklayer::import(peaks_file)
  unique_peak.gr <- data.frame(peaks.gr) %>% distinct(seqnames, start, end) %>% as.data.frame() %>% GenomicRanges::makeGRangesFromDataFrame()
  
  
  # Sinto fragments
  fragments <- data.table::fread(frags_file, select=1:5, header = FALSE, col.names = c("seqnames", "start", "end", "barcode", "count")) 
  fragments[]
  fragments.gr <- fragments %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # Get all barcodes
  barcodes <- unique(fragments$barcode)
  
  mito_counts <- fragments %>% group_by(barcode) %>% mutate(count = sum(seqnames == "chrM")) %>% select(barcode, count)
  
  # Compute overlaps of fragments with each feature
  peak.overlaps <- GenomicRanges::findOverlaps(query = fragments.gr,
                                               subject = unique_peak.gr,
                                               type = "any",
                                               ignore.strand = TRUE)
  pro.overlaps <- GenomicRanges::findOverlaps(query = fragments.gr,
                                              subject = pro.gr,
                                              type = "any",
                                              ignore.strand = TRUE)
  enhs.overlaps <- GenomicRanges::findOverlaps(query = fragments.gr,
                                               subject = enhs.gr,
                                               type = "any",
                                               ignore.strand = TRUE)
  tss.overlaps <- GenomicRanges::findOverlaps(query = fragments.gr,
                                              subject = tss.gr,
                                              type = "any",
                                              maxgap = 1000,
                                              ignore.strand = TRUE)
  
  total_counts <- table(fragments$barcode)
  peak_counts <- table(fragments.gr[unique(queryHits(peak.overlaps))]$barcode)
  promoter_counts <- table(fragments.gr[unique(queryHits(pro.overlaps))]$barcode)
  enhancer_counts <- table(fragments.gr[unique(queryHits(enhs.overlaps))]$barcode)
  tss_counts <- table(fragments.gr[unique(queryHits(tss.overlaps))]$barcode)
  
  qc_table <- data.frame(total_frags = rep(0, length(barcodes)),
                         frac_peak = rep(0, length(barcodes)),
                         frac_promoter = rep(0, length(barcodes)),
                         frac_tss = rep(0, length(barcodes)),
                         frac_enhancer = rep(0, length(barcodes)),
                         frac_mito = rep(0, length(barcodes)))
  
  rownames(qc_table) <- barcodes
  
  qc_table[names(total_counts), "total_frags"] <- total_counts
  qc_table[names(peak_counts), "frac_peak"] <- peak_counts
  qc_table[names(promoter_counts), "frac_promoter"] <- promoter_counts  
  qc_table[names(enhancer_counts), "frac_enhancer"] <- enhancer_counts
  qc_table[names(tss_counts), "frac_tss"] <- tss_counts
  qc_table[mito_counts$barcode, "frac_mito"] <- mito_counts$count
  
  qc_table$frac_peak <- qc_table$frac_peak/qc_table$total_frags
  qc_table$frac_promoter <- qc_table$frac_promoter/qc_table$total_frags
  qc_table$frac_enhancer <- qc_table$frac_enhancer/qc_table$total_frags
  qc_table$frac_tss <- qc_table$frac_tss/qc_table$total_frags
  qc_table$frac_mito <- qc_table$frac_mito/qc_table$total_frags
  
  qc_table <- qc_table[qc_table$total_frags > 5, ] %>% rownames_to_column("bc")
  
  write.csv(qc_table, file = "/stornext/Home/data/allstaff/y/yang.p/scPipe_testing/scPipe_atac_output/cell_qc_metrics.csv", row.names = FALSE)
}
