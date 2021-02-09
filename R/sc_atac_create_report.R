#' sc_atac_create_report()
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


sc_atac_create_report <- function(input_folder, 
                                  output_folder, 
                                  sample_name  = NULL, 
                                  organism     = NULL, 
                                  feature_type = NULL){
  
  # log_and_stats_folder and output_folder must have a full path ? not really..
  
  if(input_folder == ''){
    input_folder <- file.path(getwd(), "scPipe-atac-output/scPipe_atac_stats")
  }
  
  if (!dir.exists(input_folder)){
    cat("Default input folder could not be found at " , file.path(getwd()),  "Please enter the full input path to proceed \n");
    break;
  }
  
  if(output_folder == ''){
    output_folder <- file.path(getwd(), "scPipe-atac-output/scPipe_atac_stats")
  }
  
  if (!dir.exists(output_folder)){
    cat("Output folder could not be found at " , output_folder,  "Creating one now... \n");
    dir.create(output_folder, recursive=TRUE)
  }
  
  rmarkdown::render(
    input  <- "R/rmd_report_skeleton.Rmd",
    params <- list(
      log_and_stats_folder = input_folder
    ),
    output_file <- paste0(output_folder, "/scPipe_atac_report.html")
  )
}

# # Example:
# 
# sc_atac_create_report(log_and_stats_folder = here::here("out/test_out_bam/log_and_stats/"))
# 
# sc_atac_create_report(
#   log_and_stats_folder = here::here("out/test_out_bam/log_and_stats/"),
#   output_folder = here::here("out/test_out_bam/log_and_stats/"))
