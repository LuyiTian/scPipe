#' sc_atac_call_peaks()
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

sc_atac_call_peaks = function(inbam, output_folder = ""){
  
  if(output_folder == ''){
    output_folder = file.path(getwd(), "output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  system2("macs2", c("callpeak", "-t", inbam , "--outdir", output_folder)) 
  
}