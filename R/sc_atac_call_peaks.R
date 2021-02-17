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

sc_atac_call_peaks = function(inbam, ref = NULL, output_folder = NULL, genome_size = NULL){
  if (is.null(ref) && is.null(genome_size)) {
    stop("No genome or genome size was specified. Must specify at least one!")
  } else if (!file.exists(ref)) {
    stop("Specified reference file doesn't exist!")
  } else if (is.null(genome_size)) {
    genome_size = signif(file.info(ref)$size, 2)
    cat("Using genome size of", genome_size, "\n")
  }
  
  
  if(is.null(output_folder)) {
    output_folder = file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  library(MACSr)
  MACSr::callpeak(inbam, nomodel = TRUE, shift = 100, extsize = 200, gsize=genome_size, outdir = output_folder)
}