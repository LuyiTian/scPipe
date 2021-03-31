##########################################################
# Aligning Demultiplxed FASTQ Reads to a Known Reference
##########################################################

#' @name sc_atac_aligning
#' @title aligning the demultiplexed FASTQ reads
#' @description after we run the \code{sc_atac_trim_barcode} to demultiplex the fastq files, we are using this
#' function to align those fastq files to a known reference.
#' @param ref the reference genome file (.fasta, .fa format)
#' @param readFile1 the first fastq file which is mandatory
#' @param readFile2 the second fastq file, which is required if the data is paired-end
#' @param index_path if the Rsubread genome build is available user can enter the path here
#' @return 
#'
#' @examples
#' \dontrun{
#' sc_atac_aligning(ref
#'     readFile1  
#'     readFile2 
#'     nthreads  = 6) 
#' }
#'@export

sc_atac_aligning <- function (ref, 
                              readFile1, 
                              readFile2     = NULL, 
                              output_folder = NULL, 
                              output_file   = NULL,
                              input_format  = "FASTQ",
                              output_format = "BAM",
                              index_path    = NULL,
                              type          = "dna",
                              nthreads      = 1){
  
  if(is.null(output_folder)) {
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
  
  log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file             <- paste0(log_and_stats_folder, "log_file.txt")
  stats_file           <- paste0(log_and_stats_folder, "stats_file_align.txt")
  if(!file.exists(log_file)) file.create(log_file)
  # file.create(stats_file)
  
  cat(
    paste0(
      "sc_atac_aligning starts at ",
      as.character(Sys.time()),
      "\n"
    ), 
    file = log_file, append = TRUE)
  
  
  # creating an index if not available
  if (!is.null(index_path)) {
    indexPath <- index_path
    if (!file.exists(paste0(indexPath, ".log"))) {
      stop("Genome index does not exist in the specificed location. Please check the full index path again.")   
      }
  } else {
    cat("Genome index location not specified. Looking for the index in", output_folder, "\n")
    indexPath <-  file.path(output_folder, "genome_index") 
    if (file.exists(paste0(indexPath, ".log"))) {
      cat("Genome index foound in ", output_folder, "...\n")
    } else {
      cat("Genome index not foound. Creating one in ", output_folder, ". This will take a while ...\n")
      Rsubread::buildindex(basename=indexPath, reference=ref)
    }
  }
  
  # Generate the output filename
  if (is.null(output_file)) {
    # Only exception is if filename (excluding directory name) contains '.' then will only extract the first part
    fileNameWithoutExtension <- paste0(output_folder, "/", strsplit(basename(readFile1), "\\.")[[1]][1])
    outbam                   <- paste0(fileNameWithoutExtension, "_aligned.bam")
    cat("Output file name is not provided. Aligned reads are saved in ", outbam, "\n")
  }
  else {
    fileNameWithoutExtension <- paste(output_folder, strsplit(output_file, "\\.")[[1]][1], sep = "/")
    outbam                   <- paste0(output_folder, "/", output_file)
  }
  
  #execute Rsubread align()
  align_output_df <- Rsubread::align(
    index       = indexPath,
    readfile1   = readFile1,
    readfile2   = readFile2,
    sortReadsByCoordinates = TRUE,
    output_file = outbam)
  
  write.csv(align_output_df, file = stats_file, row.names = TRUE, quote = FALSE)
  
  #generating the bam index
  #Rsamtools::indexBam(paste0(output_folder, "/", fileNameWithoutExtension, "_aligned.bam"))
  Rsamtools::indexBam(paste0(fileNameWithoutExtension, "_aligned.bam"))
  
  # get the unmapped mapped stats to be output and stored in a log file
  #bamstats <- Rsamtools::idxstatsBam(paste0(output_folder, "/",fileNameWithoutExtension, "_aligned.bam"))
  bamstats <- Rsamtools::idxstatsBam(paste0(fileNameWithoutExtension, "_aligned.bam"))
  write.csv(bamstats, file = paste0(log_and_stats_folder, "stats_file_align_per_chrom.csv"), row.names = FALSE, quote = FALSE)
  
  cat(
    paste0(
      "sc_atac_aligning finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = TRUE)
  
  return(align_output_df)
} 