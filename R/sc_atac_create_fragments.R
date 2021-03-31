#' @name sc_atac_create_fragments
#' @title generating the popular fragments for scATAC-Seq data using sinto
#' @description Takes in a tagged and sorted BAM file and outputs the associated fragments in a .bed file
#' @param
#' @param 
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

# 
sc_atac_create_fragments <- function(inbam, output_folder = "") {
  # Need to use sinto to generate fragment file
  reticulate::use_virtualenv("scPipe_env") 
  reticulate::import("sinto")
  # Check if output directory exists
  if(output_folder == ''){
    output_folder = file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  system2("sinto", c("fragments", "-b", inbam, "-f", paste0(output_folder, "/fragments.bed")))
  
}