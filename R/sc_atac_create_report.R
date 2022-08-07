
###############################################
# Create the final report for scATAC-Seq data
###############################################

#' @name sc_atac_create_report
#' @title HTML report generation
#' @description Generates a HTML report using the output folder produced by the pipeline
#'
#' @param input_folder The path of the folder produced by the pipeline
#' @param feature_file is the peak or genome_bin file that was generated via ScPipe workflow or outside the workflow (.bed format)
#' @param output_folder The path of the output folder to store the HTML report in
#' @param organism A string indicating the name of the organism being analysed
#' @param sample_name A string indicating the name of the sample
#' @param feature_type A string indicating the type of the feature (`genome_bin` or `peak`)
#' @param n_barcode_subset if you require only to visualise stats for a sample of barcodes to improve processing time (integer)
#' @param tss_file bed file of the TSS regions
#' @param promoter_file bed file of the promoter regions
#' @param enhancer_file bed file fo the enhancer regions
#' @export
#'


sc_atac_create_report <- function(input_folder, 
                                  feature_file,
                                  output_folder    = NULL, 
                                  organism         = NULL,
                                  sample_name      = NULL,
                                  feature_type     = NULL,
                                  n_barcode_subset = NULL,
                                  tss_file         = NULL,
                                  promoter_file    = NULL,
                                  enhancer_file    = NULL){
  
  if (!dir.exists(input_folder)){
    stop("The input folder could not be found at ", input_folder);
  }
  
  if (!dir.exists(file.path(input_folder, "scPipe_atac_stats"))) {
    stop("The input stats folder could not be found at ", 
         file.path(input_folder, "scPipe_atac_stats"), 
         ". Either run the pipeline to generate it or create a `scPipe_atac_stats` folder that contains the required files.");
  }
  
  available_organisms = c("hg19",
                          "hg38",
                          "mm10")
  
  if(!is.null(organism) && (organism %in% available_organisms) && (is.null(promoter_file))) {
    message(organism, " is a recognized organism. Using promoter files in repository ...\n")
    anno_paths    <- system.file("extdata/annotations/", package = "scPipe", mustWork = TRUE)
    promoter_file <- file.path(anno_paths, paste0(organism, "_promoter.bed.gz"))
  }
  
  if(!is.null(organism) && (organism %in% available_organisms) && (is.null(tss_file))) {
    message(organism, " is a recognized organism. Using TSS files in repository ...\n")
    anno_paths    <- system.file("extdata/annotations/", package = "scPipe", mustWork = TRUE)
    tss_file      <- file.path(anno_paths, paste0(organism, "_tss.bed.gz"))
  }
  
  if(!is.null(organism) && (organism %in% available_organisms) && (is.null(enhancer_file))) {
    message(organism, " is a recognized organism. Using enhancer files in repository ...\n")
    anno_paths    <- system.file("extdata/annotations/", package = "scPipe", mustWork = TRUE)
    enhancer_file <- file.path(anno_paths, paste0(organism, "_enhancer.bed.gz"))
  }

  
  if(is.null(output_folder)) {
    output_folder <- file.path(input_folder, "scPipe_atac_stats")
  }
  
  if (!dir.exists(output_folder)){
    cat("Output folder could not be found at " , output_folder,  "Creating it now... \n");
    dir.create(output_folder, recursive=TRUE)
  }
  
  report.file <- system.file("extdata/rmd_report_skeleton.Rmd", package = "scPipe") 
  if (report.file == "") # Used if installed via the repo
    report.file <- system.file("inst/extdata/rmd_report_skeleton.Rmd", package = "scPipe", mustWork = TRUE)
  
  rmarkdown::render(
    input  = report.file,
    params = list(
      input_folder     = input_folder,
      organism         = organism,
      sample           = sample_name,
      feature_type     = feature_type,
      feature_file     = feature_file,
      n_barcode_subset = n_barcode_subset,
      tss_file         = tss_file,
      promoter_file    = promoter_file,
      enhancer_file    = enhancer_file
    ),
    output_file = file.path(output_folder, "scPipe_atac_report.html")  
  )
}
