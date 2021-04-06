###########################################################
# Create a SingleCell Experiment Object for scATAC-Seq data
###########################################################

#' sc_atac_create_sce()
#'
#' @return 
#'
#' @examples
#' \dontrun{
#' 
#' 
#' }
#'
#' @export
#'

sc_atac_create_sce <- function(input_folder = NULL, 
                               organism     = NULL, 
                               feature_type = NULL, 
                               pheno_data   = NULL, 
                               report       = FALSE) {
  
  if(is.null(input_folder)){
    input_folder <- file.path(getwd(), "scPipe-atac-output")
    input_stats_folder <- file.path(getwd(), "scPipe-atac-output/scPipe_atac_stats")
  } else {
    cat("input location ", input_folder, " is not valid. Please enter the full path to proceed. \n")
    break;
  }
  
  if (!dir.exists(input_folder)){
    cat("Default input folder could not be found at " , file.path(getwd()),  "Please enter the full input path to proceed \n");
    break;
  } else {
    input_stats_folder <- file.path(input_folder, "scPipe_atac_stats")
  }
  
  #feature_cnt   <- readMM(file.path(input_folder, "sparse_matrix.mtx"))
  feature_cnt   <- readRDS(file.path(input_folder, "sparse_matrix.rds"))
  cell_stats    <- read.csv(file.path(input_stats_folder, "filtered_stats_per_cell.csv"), row.names=1)
  feature_stats <- read.csv(file.path(input_stats_folder, "filtered_stats_per_feature.csv"))
  
  # need to change from here.... (check whether I need to filter before saving to the SCE object)
  
  # can I order a matrix like a csv file like below? test...
  feature_cnt      <- feature_cnt[, order(colnames(feature_cnt))]
  cell_stats       <- cell_stats[order(rownames(cell_stats)), ]
  
  # generating the SCE object
  sce                         <- SingleCellExperiment(assays = list(counts = as.matrix(feature_cnt)))
  sce@metadata$scPipe$version <- packageVersion("scPipe")  # set version information
  
  if(!is.null(organism)){
    organism(sce) <- organism
  }
  if(!is.null(feature_type)){
    feature_type(sce) <- feature_type
  }
  
  QC_metrics(sce) <- cell_stats
  
  if(!is.null(pheno_data)){
    colData(sce) <- cbind(colData(sce), pheno_data[order(rownames(pheno_data)),])
  }
  
  #feature_info(sce) <- feature_stats
  
  saveRDS(sce, file = paste(input_folder,"/scPipe_atac_SCEobject.rds",sep = ""))
  
  if(report){
    sc_atac_create_report(input_folder = file.path(getwd(), "scPipe-atac-output/scPipe_atac_stats"),
                          output_folder= input_folder,
                          sample_name  = NULL,
                          organism     = organism,
                          feature_type = feature_type)
  }
  
  
  return(sce)
  
}
