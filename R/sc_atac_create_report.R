#' @name sc_atac_create_report
#' @title HTML report generation
#' @description Generates a HTML report using the output folder produced by the pipeline
#'
#' @param input_folder The path of the folder produced by the pipeline
#' @param output_folder The path of the output folder to store the HTML report in
#' @param sample_name A string indicating the name of the sample
#' @param organism A string indicating the name of the organism being analysed
#' @param feature_type A string indicating the type of the feature (`genome_bin` or `peak`)
#'
#' @export
#'


sc_atac_create_report <- function(input_folder, 
                                  output_folder = NULL, 
                                  sample_name  = NULL, 
                                  organism     = NULL, 
                                  feature_type = NULL){
  
  if (!dir.exists(input_folder)){
    stop("The input folder could not be found at ", input_folder);
  }
  
  if (!dir.exists(file.path(input_folder, "scPipe_atac_stats"))) {
    stop("The input stats folder could not be found at ", 
         file.path(input_folder, "scPipe_atac_stats"), 
         ". Either run the pipeline to generate it or create a `scPipe_atac_stats` folder that contains the required files.");
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
      input_folder = input_folder,
      organism = organism
    ),
    output_file = file.path(output_folder, "scPipe_atac_report.html")  
  )
}
