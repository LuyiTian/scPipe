#####################################################
# Integrate multi-omic scRNA-Seq and scATAC-Seq data
#####################################################

#' @name sc_integrate
#' @title Integrate multi-omic scRNA-Seq and scATAC-Seq data
#' @description Generates an integrated SCE object with scRNA-Seq and scATAC-Seq data produced by the scPipe pipelines
#' @param sce_rna SCE object for scRNA-Seq data
#' @param sce_atac SCE object for scATAC-Seq data
#' @param output_folder The path of the output folder to store the integrated SCE object in
#' @param barcode_match_file A .csv file with scRNA-Seq based barcodes in column 1 and scATAC-Seq barcodes in column 2 with a header row
#' @param rna_barcode A .csv file (w/ header row) that is output from scPipe RNA-Seq preprocessing pipeline that has cell name and barcode assignment (cell name should be the first column and barcode sequence in the second column). If the colnames of the inputing SCE for RNA-Seq is different to what is in the `barcode_match_file`, ths file is needed (optional)
#' @param atac_barcode A .csv file (w/ header row) that has assignment of barcode to different cell names assigned in the ATAC-Seq based SCE object (cell name should be the first column and barcode sequence in the second column). If the colnames of the inputing SCE for ATAC-Seq is different to what is in the `barcode_match_file`, ths file is needed (optional)
#' @param atac_revcomp A logical parameter if TRUE will convert the ATAC-Seq barcode list to revcomp(); Default: FALSE
#' @export
#'

sc_integrate <- function(sce_rna, 
                         sce_atac,
                         barcode_match_file,
                         rna_barcode     = NULL,
                         atac_barcode    = NULL,
                         atac_revcomp    = FALSE,
                         output_folder   = NULL){
  
  # sanity check to know whether all files are avaialble
  if(!all(file.exists(sce_rna, sce_atac, barcode_match_file))){
    stop("At least one of the input files for R1 does not exist")    
  }
  
  if(is.null(output_folder)) {
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory is not provided. Created directory: ", output_folder, "\n")
  }
  
  # load the SCE objects
  sce_rna  <- readRDS(sce_rna)
  sce_atac <- readRDS(sce_atac) 
  
  # first rename the assay for ATAC
  assay(sce_atac, "ATAC") <- counts(sce_atac)
  assay(sce_rna, "RNA")   <- counts(sce_rna)
  
  # read in the barcode file
  barcode_match_file <- read.csv(barcode_match_file, header = TRUE)
  
  # >>> sanity check to see whether the info in the barcode file is the same as column names of the SCE object here <<<
  
  if(!is.null(rna_barcode)) {
    # >>>  match the barcodes in the SCE object here <<<
    rna_bc_anno                                           <- read.csv(rna_barcode, header = TRUE)
    colnames(rna_bc_anno)[1:2]                            <- c("cell_name", "barcode_sequence")
    colnames(assay(sce_rna, withDimnames = FALSE, "RNA")) <- rna_bc_anno$barcode_sequence[match(colnames(assay(sce_rna, "RNA")), rna_bc_anno$cell_name)]
    rownames(colData(sce_rna))                            <- rna_bc_anno$barcode_sequence[match(rownames(colData(sce_rna)), rna_bc_anno$cell_name)]
  }
  
  
  if(!is.null(atac_barcode)) {
    # >>>  match the barcodes in the SCE object here <<<
    atac_bc_anno                                            <- read.csv(atac_barcode, header = TRUE)
    colnames(atac_bc_anno)[1:2]                             <- c("cell_name", "barcode_sequence")
    colnames(assay(sce_atac, withDimnames = FALSE, "ATAC")) <- atac_bc_anno$barcode_sequence[match(colnames(assay(sce_atac, "ATAC")), atac_bc_anno$cell_name)]
    rownames(colData(sce_atac))                             <- atac_bc_anno$barcode_sequence[match(rownames(colData(sce_atac)), atac_bc_anno$cell_name)]
  }
  
  if (atac_revcomp) {
    atac_bc_anno@barcode_sequence <- Biostrings::reverseComplement(atac_bc_anno@barcode_sequence)
  }
  
  # now we can use this to join the RNA-Seq and ATAC-Seq
  
  # first create a sep. data frame with matching RNA and ATAC barcodes and filter the two SCE objects based on that
  rna_vector  <- barcode_match_file$rna[colnames(assay(sce_rna, withDimnames = FALSE, "RNA"))]
  atac_vector <- barcode_match_file$atac[colnames(assay(sce_atac, withDimnames = FALSE, "ATAC"))]
  
  
  
}