#' sc_atac_align()
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

sc_atac_aligning <- function (ref, 
                              readFile1, 
                              readFile2     = NULL, 
                              readDir       = NULL, 
                              output_folder = "", 
                              output_file   = '',
                              input_format  = "FASTQ",
                              output_format = "BAM",
                              type          = "dna",
                              nthreads      = 1){
  
  if(output_folder == ''){
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory is not provided. Created directory: ", output_folder, "\n")
  }
  
  if(!file.exists(readFile1)){
    stop("Input File readFile1 does not exist")    
  }
  
  if(!is.null(readFile2) && !file.exists(readFile2)){
    stop("Input File readFile2 does not exist")    
  }
  
  log_and_stats_folder = paste0(output_folder, "/log_and_stats/")
  dir.create(log_and_stats_folder, showWarnings = F)
  log_file = paste0(log_and_stats_folder, "log_file.txt")
  stats_file = paste0(log_and_stats_folder, "stats_file_align.txt")
  if(!file.exists(log_file)) file.create(log_file)
  # file.create(stats_file)
  
  cat(
    paste0(
      "sc_atac_aligning starts at ",
      as.character(Sys.time()),
      "\n"
    ), 
    file = log_file, append = T)
  
  
  
  # creating an index
  indexPath <-  file.path(output_folder, "genome_index") 
  buildindex (basename=indexPath,
              reference=ref)
  
  if (!is.null(output_file)) {
    fileNameWithoutExtension <- strsplit(readFile1, "\\.")[[1]][1]
    outbam                   <- paste(fileNameWithoutExtension, "_aligned.bam", sep = "")
    cat("Output file name is not provided. Aligned reads are saved in ", outbam, "\n")
  }
  else{
    fileNameWithoutExtension <- strsplit(output_file, "\\.")[[1]][1]
    outbam                   <- paste(output_folder, "/", output_file, sep="")
  }
  
  #execute Rsubread align()
  align_output_df = Rsubread::align(
    index = indexPath,
    readfile1 = readFile1,
    readfile2 = readFile2, 
    output_file = outbam)
  
  write.csv(align_output_df, file = stats_file)
  
  #generating the bam index
  Rsamtools::sortBam(outbam, paste0(fileNameWithoutExtension, "_aligned_sorted"))
  #Rsamtools::indexBam(outbam)
  Rsamtools::indexBam(paste0(fileNameWithoutExtension, "_aligned_sorted.bam"))
  
  # get the unmapped mapped stats to be output and stored in a log file
  #can use Rsamtools::idxstatsBam()
  bamstats <- Rsamtools::idxstatsBam(paste0(fileNameWithoutExtension, "_aligned_sorted.bam"))
  write.csv2(bamstats, paste(output_folder, "/", fileNameWithoutExtension, "_alignment_stats.csv"), sep = "")
  
  cat(
    paste0(
      "sc_atac_aligning finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = T)
  
  return(align_output_df)
} 





