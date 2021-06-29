#' sc_atac_create_report()
#'
#' @return 
#'
#' @examples
#' sc_atac_create_report(log_and_stats_folder = here::here("out/test_out_bam/log_and_stats/"))
# 
#  sc_atac_create_report(
#    log_and_stats_folder = here::here("out/test_out_bam/log_and_stats/"),
#    output_folder = here::here("out/test_out_bam/log_and_stats/"))
#' \dontrun{
#' 
#' 
#' }
#'
#' @import data.table
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
  
  rmarkdown::render(
    input  = here::here("R", "rmd_report_skeleton.Rmd"),
    params = list(
      input_folder = input_folder,
      organism = organism
    ),
    output_file = file.path(output_folder, "scPipe_atac_report.html")  
  )
}

# # Example:
# 
# sc_atac_create_report(log_and_stats_folder = here::here("out/test_out_bam/log_and_stats/"))
# 
# sc_atac_create_report(
#   log_and_stats_folder = here::here("out/test_out_bam/log_and_stats/"),
#   output_folder = here::here("out/test_out_bam/log_and_stats/"))
