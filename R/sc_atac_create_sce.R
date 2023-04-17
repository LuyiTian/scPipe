###########################################################
# Create a SingleCell Experiment Object for scATAC-Seq data
###########################################################

#' sc_atac_create_sce()
#'
#' @param input_folder The output folder produced by the pipeline
#' @param organism The type of the organism
#' @param sample_name The name of the sample
#' @param feature_type The type of the feature
#' @param pheno_data The pheno data
#' @param report Whether or not a HTML report should be produced
#'
#' @returns a SingleCellExperiment object created from the scATAC-Seq data provided
#'
#' @importFrom utils read.csv
#' @importFrom tibble column_to_rownames
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @examples
#' \dontrun{
#' sc_atac_create_sce(
#'    input_folder = input_folder,
#'    organism = "hg38",
#'    feature_type = "peak",
#'    report = TRUE)
#' }    
#' 
#' @export
#'
sc_atac_create_sce <- function(input_folder  = NULL, 
                                organism      = NULL, 
                                sample_name   = NULL,
                                feature_type  = NULL, 
                                pheno_data    = NULL, 
                                report        = FALSE) {
    
    # if(is.null(input_folder)){
    #   input_folder <- file.path(getwd(), "scPipe-atac-output")
    #   input_stats_folder <- file.path(getwd(), "scPipe-atac-output/scPipe_atac_stats")
    # } else {
    #   cat("input location ", input_folder, " is not valid. Please enter the full path to proceed. \n")
    #   break;
    # }
    
    if(is.null(input_folder)){
        input_folder <- file.path(getwd(), "scPipe-atac-output")
        input_stats_folder <- file.path(getwd(), "scPipe-atac-output/scPipe_atac_stats")
    }
    
    if (!dir.exists(input_folder)){
        message("Default input folder could not be found at " , input_folder,  "\nPlease enter the full input path to proceed");
    } else {
        input_stats_folder <- file.path(input_folder, "scPipe_atac_stats")
    }
    
    #feature_cnt   <- readMM(file.path(input_folder, "sparse_matrix.mtx"))
    feature_cnt   <- readRDS(file.path(input_folder, "sparse_matrix.rds"))
    cell_stats    <- utils::read.csv(file.path(input_stats_folder, "filtered_stats_per_cell.csv"), row.names=1)
    feature_stats <- utils::read.csv(file.path(input_stats_folder, "filtered_stats_per_feature.csv"))
    

    # need to change from here.... (check whether I need to filter before saving to the SCE object)
    
    # can I order a matrix like a csv file like below? test...
    feature_cnt      <- feature_cnt[, order(c(colnames(feature_cnt)))]
    
    
    qc <- utils::read.csv(file.path(input_folder, "cell_qc_metrics.csv")) %>% tibble::column_to_rownames(var = "bc")
    cell_stats <- merge(x = cell_stats, y = qc, by = 0, all.x = TRUE) %>% tibble::column_to_rownames(var = "Row.names")
    cell_stats       <- cell_stats[order(c(rownames(cell_stats))), ]
    

    # generating the SCE object
    sce                         <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = feature_cnt))

    sce@metadata$scPipe$version <- packageVersion("scPipe")  # set version information
    
    if(!is.null(organism)){
        organism(sce) <- organism
    }
    
    if(!is.null(feature_type)){
        feature_type(sce) <- feature_type
    }
    
    # Saving demultiplexing stats to sce object
    stats_file <- file.path(input_stats_folder, "mapping_stats_per_barcode.csv")
    raw <- read.csv(stats_file, row.names = "barcode")
    
    
    QC_metrics(sce) <- cell_stats
    demultiplex_info(sce) <- raw
    
    
    
    if(!is.null(pheno_data)){
        colData(sce) <- cbind(colData(sce), pheno_data[order(rownames(pheno_data)),])
    }

    feature_info(sce) <- feature_stats
    saveRDS(sce, file = file.path(input_folder, "scPipe_atac_SCEobject.rds"))
    
    message("SCE object is saved in:\n\t", input_folder,"/scPipe_atac_SCEobject.rds")
    
    if(report){
        sc_atac_create_report(input_folder  = file.path(input_folder),
                            output_folder = file.path(input_folder, "scPipe_atac_stats"),
                            sample_name   = sample_name,
                            organism      = organism,
                            feature_type  = feature_type)
    }
    
    
    return(sce)
    
}

