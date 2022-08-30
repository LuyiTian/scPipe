
###########################################################
# scATAC-Se and scRNA-Se data integration for multi-omic data
###########################################################

#' @name sc_multiomic_integration
#'
#' @title Inegrating modalities from multi-omic data using barcode matching
#' @description Using a barcode file match modalities of multi-omic datasets 
#' 
#' @param sce.atac SCE object pre-processed for scATAC-Seq data
#' @param sce.rna SCE object pre-processed for scRNA-Seq data
#' @param barcode.match.file A .csv file that contains the barcode matching form each of the modalities; currently supporting only RNA and ATAC
#' @param output_folder The path of the output folder 
#'
#' @examples
#' \dontrun{
#' sc_multiomic_integration(sce.atac = sceatac, sce.rna = scerna, barcode.match.file = barcode_file) 
#' }
#' 
#' @export
#' 

sc_multiomic_integration <- function(sce.atac, sce.rna, output_folder = NULL){
  
  # Check if output directory exists
  if (is.null(output_folder)) {
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory does not exist. Created directory: ", output_folder, "\n")
  }
  
  #Create log folder/file
  log_and_stats_folder <- file.path(output_folder, "scPipe_atac_stats")
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file            <- file.path(log_and_stats_folder, "log_file.txt")
  if(!file.exists(log_file)) file.create(log_file)
  cat(
    paste0(
      "sc_atac_tfidf starts at ",
      as.character(Sys.time()),
      "\n"
    ),
    file = log_file, append = TRUE)
  
  # sanity check to know whether both are SCE objects
  if ((class(x = sce.atac) != "SingleCellExperiment") && (class(x = sce.rna) != "SingleCellExperiment")) { 
    stop("One of the SCE objects are corrupted or don' exist. Please check them again before running this fucnction \n") 
  }
  
  
  

}