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

sc_atac_peak_calling <- function(inbam, output_folder = ""){
  
  if(output_folder == ''){
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
  
  system2("macs2", c("callpeak", "-t", inbam , "--outdir", output_folder)) 
  
  # add the following to extend peaks
  
  # dd = fread(peakFile)
  # 
  # names(dd)[2:3] = c('start', 'end')
  # dd[, 'ss' := (end - start)]
  # dd[, 'midp' := floor(end/2 + start/2)]
  # dd[, 'start' := ifelse(ss < 500, midp - 250, start)]
  # dd[, 'end' := ifelse(ss < 500, midp + 250, end)]
  # dd[, c('ss', 'midp') := NULL]
  # 
  # write.table(dd, file = peakFile, row.names = F, col.names = F, quote = F,
  #             sep = '\t')
  
}