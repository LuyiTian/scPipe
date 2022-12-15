#' @name sc_atac_create_cell_qc_metrics
#' @title generating a file useful for producing the qc plots
#' @description uses the peak file and annotation files for features
#' 
#' @param frags_file The fragment file
#' @param peaks_file The peak file
#' @param promoters_file The path of the promoter annotation file 
#' @param tss_file The path of the tss annotation file 
#' @param enhs_file The path of the enhs annotation file 
#' @param output_folder The path of the output folder for resultant files
#' 
#' @importFrom data.table fread setkey copy :=

#' 
#' @returns Nothing (Invisible 'NULL')
#' @export
#' 
sc_atac_create_cell_qc_metrics <- function(frags_file,
                                            peaks_file,
                                            promoters_file,
                                            tss_file,
                                            enhs_file,
                                            output_folder) {
    
    pro.gr <- rtracklayer::import(promoters_file)
    tss_df <- data.table::fread(tss_file, select=c(seq_len(3)), header = FALSE, col.names = c("chr", "start", "end"))
    tss.gr <- GenomicRanges::makeGRangesFromDataFrame(tss_df)
    enhs.gr <- rtracklayer::import(enhs_file)
    
    # Peak file
    peaks.gr <-  rtracklayer::import(peaks_file)
    unique_peak.gr <- data.frame(peaks.gr) %>% dplyr::distinct(seqnames, start, end) %>% as.data.frame() %>% GenomicRanges::makeGRangesFromDataFrame()
    
    
    # Sinto fragments
    fragments <- data.table::fread(frags_file, select=seq_len(5), header = FALSE, col.names = c("seqnames", "start", "end", "barcode", "count")) 
    fragments.gr <- fragments %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # Get all barcodes
    barcodes <- unique(fragments$barcode)
    
    mito_counts <- fragments %>% dplyr::group_by(barcode) %>% dplyr::mutate(count = sum(seqnames == "chrM")) %>% dplyr::select(barcode, count)
    
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
    peak_counts <- table(fragments.gr[unique(S4Vectors::queryHits(peak.overlaps))]$barcode)
    promoter_counts <- table(fragments.gr[unique(S4Vectors::queryHits(pro.overlaps))]$barcode)
    enhancer_counts <- table(fragments.gr[unique(S4Vectors::queryHits(enhs.overlaps))]$barcode)
    tss_counts <- table(fragments.gr[unique(S4Vectors::queryHits(tss.overlaps))]$barcode)
    
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
    
    qc_table <- qc_table[qc_table$total_frags > 5, ] %>% tibble::rownames_to_column("bc")
    
    utils::write.csv(qc_table, file = file.path(output_folder, "cell_qc_metrics.csv"), row.names = FALSE)
}
