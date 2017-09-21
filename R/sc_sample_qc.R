#' @name sc_sample_qc
#' @title quality control information for a small sample scRNA-seq dataset to 
#' demonstrate capabilities of scPipe
#' @description This data.frame contains cell quality control information for the 100 
#' cells.
#' @return NULL, but makes a data frame with cell quality control data.frame
#' @docType data
#' @usage sc_sample_qc
#' @format a data.frame instance, 1 row per cell.
#' @source Christin Biben, WEHI
#' @author Luyi Tian
NULL

#' @name UMI_duplication
#' @title UMI dulication statistics for a small sample scRNA-seq dataset to 
#' demonstrate capabilities of scPipe
#' @description This data.frame contains UMI dulication statistics, where the
#' first column is the number of duplication, and the second column is the
#' count of UMIs.
#' 
#' @return NULL, but makes a data frame with UMI dulication statistics
#' @docType data
#' @usage sc_sample_qc
#' @format a data.frame instance.
#' @source Christin Biben, WEHI
#' @author Luyi Tian
NULL

#' @name cell_barcode_matching
#' @title cell barcode demultiplex statistics for a small sample scRNA-seq dataset to 
#' demonstrate capabilities of scPipe
#' @description This data.frame contains cell barcode demultiplex statistics with 
#' several rows:
#' * barcode_unmatch_ambiguous_mapping
#' * barcode_unmatch_mapped_to_intron
#' * barcode_match
#' * barcode_unmatch_unaligned
#' * barcode_unmatch_aligned
#' * barcode_unmatch_mapped_to_exon
#' 
#' @return NULL, but makes a data frame with cell barcode demultiplex statistics
#' @docType data
#' @usage sc_sample_qc
#' @format a data.frame instance.
#' @source Christin Biben, WEHI
#' @author Luyi Tian
NULL