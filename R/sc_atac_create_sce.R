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

sc_atac_create_sce <- function(input_folder, organism=NULL, feature_type=NULL, pheno_data=NULL, report=FALSE) {
  
  if(input_folder == ''){
    input_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(input_folder)){
    cat("Default input folder could not be found at " , file.path(getwd()),  "Please enter the full input path to proceed \n");
    break;
  }
  
  feature_cnt <- readMM(file.path(input_folder, "sparse_matrix.mtx"))
  cell_stat   <- read.csv(file.path(input_folder, "scPipe_atac_stat", "filtered_stats_per_cell.csv"), row.names=1)
  
  # need to change from here....
  demultiplex_stat <- read.csv(file.path(input_folder, "stat", "overall_stat.csv"))
  UMI_dup_stat <- read.csv(file.path(input_folder, "stat", "UMI_duplication_count.csv"))
  
  gene_cnt = gene_cnt[, order(colnames(gene_cnt))]
  cell_stat = cell_stat[order(rownames(cell_stat)), ]
  
  
  sce = SingleCellExperiment(assays = list(counts =as.matrix(gene_cnt)))
  sce@metadata$scPipe$version = packageVersion("scPipe")  # set version information
  if(!is.null(organism)){
    organism(sce) = organism
  }
  if(!is.null(gene_id_type)){
    gene_id_type(sce) = gene_id_type
  }
  QC_metrics(sce) = cell_stat
  if(!is.null(pheno_data)){
    colData(sce) = cbind(colData(sce), pheno_data[order(rownames(pheno_data)),])
  }
  
  demultiplex_info(sce) = demultiplex_stat
  UMI_dup_info(sce) = UMI_dup_stat
  #if(any(grepl("^ERCC-", rownames(sce)))){
  #  isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
  #}
  
  
  if(report){
    create_report(sample_name=basename(datadir),
                  outdir=datadir,
                  organism=organism,
                  gene_id_type=gene_id_type)
  }
  
  
  return(sce)

}
