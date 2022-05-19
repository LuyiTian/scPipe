#' @name sc_atac_create_qc_per_bc_file
#' @title generating a file useful for producing the qc plots
#' @description uses the peak file and annotation files for features
#' @param inbam The input bam file
#' @param frags_file The fragment file
#' @param peaks_file The peak file
#' @param promoters_file The path of the promoter annotation file 
#' @param tss_file The path of the tss annotation file 
#' @param enhs_file The path of the enhs annotation file 
#' @param output_folder
#' 
#' @param lower the lower threshold for the data if using the \code{emptydrops} function for cell calling.
#' 
#' @param min_uniq_frags The minimum number of required unique fragments required for a cell (used for \code{filter} cell calling)
#' @param max_uniq_frags The maximum number of required unique fragments required for a cell (used for \code{filter} cell calling)
#' @param min_frac_peak The minimum proportion of fragments in a cell to overlap with a peak (used for \code{filter} cell calling)
#' @param min_frac_tss The minimum proportion of fragments in a cell to overlap with a tss (used for \code{filter} cell calling)
#' @param min_frac_enhancer The minimum proportion of fragments in a cell to overlap with a enhancer sequence (used for \code{filter} cell calling)
#' @param min_frac_promoter The minimum proportion of fragments in a cell to overlap with a promoter sequence (used for \code{filter} cell calling)
#' @param max_frac_mito The maximum proportion of fragments in a cell that are mitochondrial (used for \code{filter} cell calling)
#' @param output_folder A string indicating the path of the output folder
#' 
#' @importFrom data.table fread setkey copy :=
#' 
#' @export
#' 
sc_atac_create_qc_per_bc_file <- function(inbam, 
                                          frags_file,
                                          peaks_file,
                                          promoters_file,
                                          tss_file,
                                          enhs_file,
                                          output_folder,
                                          lower = NULL) {
  
  
  # https://github.com/wbaopaul/scATAC-pro/blob/master/scripts/src/get_qc_per_barcode.R
  
  chr <- .N <- bc <- total_frags <- isMito <- NULL
  
  out.frag.overlap.file <- file.path(output_folder, "qc_per_bc_file.txt")
  
  frags <- fread(frags_file, select=1:4, header = FALSE)
  names(frags) <- c('chr', 'start', 'end', 'bc')
  setkey(frags, chr, start)
  
  frags[, 'total_frags' := .N, by = bc]
  frags <- frags[total_frags > 5]
  
  frags <- frags[!grepl(chr, pattern = 'random', ignore.case = TRUE)]
  frags <- frags[!grepl(chr, pattern ='un', ignore.case = TRUE)]
  
  peaks <- fread(peaks_file, select=1:3, header = FALSE)
  promoters <- fread(promoters_file, select=1:3, header = FALSE)
  tss <- fread(tss_file, select=1:3, header = FALSE)
  enhs <- fread(enhs_file, select=1:3, header = FALSE)
  names(peaks) = names(promoters) = names(tss) =
    names(enhs) = c('chr', 'start', 'end')
  
  
  setkey(peaks, chr, start)
  setkey(promoters, chr, start)
  setkey(tss, chr, start)
  setkey(enhs, chr, start)
  
  chrs <- unique(frags$chr)
  
  tss[, 'start' := start - 1000]
  tss[, 'end' := end + 1000]
  fragsInRegion <- NULL
  
  for(chr0 in chrs){
    peaks0 <- peaks[chr == chr0]
    promoters0 <- promoters[chr == chr0]
    
    tss0 <- tss[chr == chr0]
    enhs0 <- enhs[chr == chr0]
    frags0 <- frags[chr == chr0]
    frags <- frags[chr != chr0]
    if(nrow(peaks0) == 0){
      frags0[, 'peaks' := 0]
    }else{
      frags0[, 'peaks' := sc_atac_getOverlaps_read2AnyRegion(frags0, peaks0)]
    }
    
    if(nrow(promoters0) == 0){
      frags0[, 'promoters' := 0]
    }else{
      frags0[, 'promoters' := sc_atac_getOverlaps_read2AnyRegion(frags0, promoters0)]
    }
    
    if(nrow(tss0) == 0){
      frags0[, 'tss' := 0]
    }else{
      frags0[, 'tss' := sc_atac_getOverlaps_read2AnyRegion(frags0, tss0)]
    }
    
    if(nrow(enhs0) == 0){
      frags0[, 'enhs' := 0]
    }else{
      frags0[, 'enhs' := sc_atac_getOverlaps_read2AnyRegion(frags0, enhs0)]
    }
    
    
    fragsInRegion <- rbind(fragsInRegion, frags0)
    message(paste(chr0, 'Done!'))
  }
  rm(frags)
  
  # get barcode region count matrix
  fragsInRegion[, 'isMito' := ifelse(chr %in% c('chrM'), 1, 0)]
  fragsInRegion[, c('chr', 'start', 'end') := NULL]
  fragsInRegion[, 'frac_mito' := sum(isMito)/total_frags, by = bc]
  fragsInRegion[, 'isMito' := NULL]
  fragsInRegion[, 'frac_peak' := sum(peaks)/total_frags, by = bc]
  fragsInRegion[, 'peaks' := NULL]
  fragsInRegion[, 'frac_promoter' := sum(promoters)/total_frags, by = bc]
  fragsInRegion[, 'promoters' := NULL]
  fragsInRegion[, 'frac_tss' := sum(tss)/total_frags, by = bc]
  fragsInRegion[, 'tss' := NULL]
  fragsInRegion[, 'frac_enhancer' := sum(enhs)/total_frags, by = bc]
  fragsInRegion[, 'enhs' := NULL]
  fragsInRegion <- unique(fragsInRegion)
  
  utils::write.table(fragsInRegion, file = out.frag.overlap.file, sep = '\t',
                     row.names = FALSE, quote = FALSE)
  
  
}
