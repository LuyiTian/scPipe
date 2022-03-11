#' @name sc_sample_qc
#' @title quality control information for a small sample scRNA-seq dataset to 
#' demonstrate capabilities of scPipe.
#' @description This data.frame contains cell quality control information 
#' for the 100 cells. For each cell it has: \itemize{
#'  \item{unaligned} the number of unaligned reads.
#'  \item{aligned_unmapped} the number of reads that aligned to genome but 
#'  fail to map to any features. 
#'  \item{mapped_to_exon} is the number of reads that mapped
#'  to exon.
#'  \item{mapped_to_intron} is the number of reads that mapped to intron.
#'  \item{ambiguous_mapping} is the number of reads that mapped to multiple
#'  features. They are not considered in the following analysis.
#'  \item{mapped_to_ERCC} is the number of reads that mapped to ERCC 
#'  spike-in controls.
#'  \item{mapped_to_MT} is the number of reads that mapped to mitochondrial
#'  genes.
#'  \item{total_count_per_cell} is the number of reads that mapped to exon
#'  after UMI deduplication. In contrast, `mapped_to_exon` is the 
#'  number of reads mapped to exon before UMI deduplication.
#'  \item{number_of_genes} is the number of genes detected for each cells
#'  \item{non_ERCC_percent} is 1 - (percentage of ERCC reads). Reads are
#'  UMI deduplicated.
#'  \item{non_mt_percent} is 1 - (percentage of mitochondrial reads).
#'  Reads are UMI deduplicated.
#'  \item{non_ribo_percent} is 1- (percentage of ribosomal reads).
#'  Reads are UMI deduplicated.
#'  }
#' @return NULL, but makes a data frame with cell quality control data.frame
#' @docType data
#' @format a data.frame instance, one row per cell.
#' @source Christin Biben (WEHI). She FACS sorted cells from several immune
#' cell types including B cells, granulocyte and some early progenitors. 
#' @author Luyi Tian
#' @examples 
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' head(QC_metrics(sce))
#' plot_mapping(sce,percentage=TRUE,dataname="sc_sample")
#' 
NULL

#' @name UMI_duplication
#' @title UMI duplication statistics for a small sample scRNA-seq dataset to 
#' demonstrate capabilities of scPipe
#' @description This data.frame contains UMI duplication statistics, where the
#' first column is the number of duplication, and the second column is the
#' count of UMIs.
#' 
#' @return NULL, but makes a data frame with UMI dulication statistics
#' @docType data
#' @format a data.frame instance, one row per cell.
#' @source Christin Biben (WEHI). She FACS sorted cells from several immune
#' cell types including B cells, granulocyte and some early progenitors. 
#' @author Luyi Tian
#' @examples 
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' sce = SingleCellExperiment(assays = list(counts =as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' demultiplex_info(sce) = cell_barcode_matching
#' UMI_dup_info(sce) = UMI_duplication
#' 
#' head(UMI_dup_info(sce))
#' 
NULL

#' @name cell_barcode_matching
#' @title cell barcode demultiplex statistics for a small sample scRNA-seq 
#' dataset to demonstrate capabilities of scPipe
#' @description This data.frame contains cell barcode demultiplex statistics 
#' with several rows:\itemize{
#' \item{barcode_unmatch_ambiguous_mapping} is the number of reads that do not
#' match any barcode, but aligned to the genome and mapped to multiple 
#' features.
#' \item{barcode_unmatch_mapped_to_intron} is the number of reads that do
#' not match any barcode, but aligned to the genome and mapped to intron.
#' \item{barcode_match} is the number of reads that match the cell barcodes
#' \item{barcode_unmatch_unaligned} is the number of reads that do not 
#' match any barcode, and not aligned to the genome
#' \item{barcode_unmatch_aligned} is the number of reads that do not 
#' match any barcode, but aligned to the genome and do not mapped to any 
#' feature
#' \item{barcode_unmatch_mapped_to_exon} is the number of reads that do not
#' match any barcode, but aligned to the genome and mapped to the exon
#' }
#' @return NULL, but makes a data frame with cell barcode demultiplex statistics
#' @docType data
#' @format a data.frame instance, one row per cell.
#' @source Christin Biben (WEHI). She FACS sorted cells from several immune
#' cell types including B cells, granulocyte and some early progenitors. 
#' @author Luyi Tian
#' @examples 
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' sce = SingleCellExperiment(assays = list(counts =as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' demultiplex_info(sce) = cell_barcode_matching
#' UMI_dup_info(sce) = UMI_duplication
#' 
#' demultiplex_info(sce)
#' 
NULL