##########################################################
# Aligning Demultiplxed FASTQ Reads to a Known Reference
##########################################################

#' @name sc_aligning
#' @title aligning the demultiplexed FASTQ reads using the Rsubread:align()
#' @description after we run the \code{sc_trim_barcode} or \code{sc_atac_trim_barcode} to demultiplex the fastq files, we are using this
#' function to align those fastq files to a known reference.
#' @param ref a character string specifying the path to reference genome file (.fasta, .fa format)
#' @param index_path character string specifying the path/basename of the index files, if the Rsubread genome build is available 
#' @param tech a character string giving the sequencing technology. Possible value includes "atac" or "rna"
#' @param R1 a mandatory character vector including names of files that include sequence reads to be aligned. For paired-end reads, this gives the list of files including first reads in each library. File format is FASTQ/FASTA by default.
#' @param R2 a character vector, the second fastq file, which is required if the data is paired-end
#' @param output_folder a character string, the name of the output folder
#' @param output_file a character vector specifying names of output files. By default, names of output files are set as the file names provided in R1 added with an suffix string
#' @param input_format a string indicating the input format
#' @param output_format a string indicating the output format
#' @param nthreads numeric value giving the number of threads used for mapping. 
#'
#' @examples
#' \dontrun{
#' sc_aligning(index_path,
#'     tech = 'atac',
#'     R1,  
#'     R2, 
#'     nthreads  = 6) 
#' }
#'@export

sc_aligning <- function (
  R1, 
  R2            = NULL,
  tech          = "atac",
  index_path    = NULL,
  ref           = NULL,
  output_folder = NULL, 
  output_file   = NULL,
  input_format  = "FASTQ",
  output_format = "BAM",
  nthreads      = 1){
  
  
  if(!all(file.exists(R1))){
    stop("At least one of the input files for R1 does not exist")    
  }
  
  if(!is.null(R2) && !all(file.exists(R2))){
    stop("At least one of input file for R2 does not exist")
  }
  
  if(tech == "atac") {
    cat("ATAC-Seq mode is selected...\n")
    
    if(is.null(output_folder)) {
      output_folder        <- file.path(getwd(), "scPipe-atac-output")
    }
    
    if (!dir.exists(output_folder)){
      dir.create(output_folder,recursive=TRUE)
      cat("Output directory is not provided. Created directory: ", output_folder, "\n")
    }
    log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
    type                 <- "dna"
  } else if(tech == "rna") {
    cat("RNA-Seq mode is selected...")
    
    if(is.null(output_folder)) {
      stop("output_folder cannot be NULL for rna mode. Aborting...\n")
    }
    log_and_stats_folder <- output_folder
    type                 <- "rna"  
  } else
  
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file             <- paste0(log_and_stats_folder, "log_file.txt")
  stats_file           <- paste0(log_and_stats_folder, "stats_file_align.txt")
  if(!file.exists(log_file)) file.create(log_file)
  
  cat(
    paste0(
      "sc_aligning starts at ",
      as.character(Sys.time()),
      "\n"
    ), 
    file = log_file, append = TRUE
  )
  
  
  # creating an index if not available
  if (is.null(index_path) && is.null(ref)) {
    stop("either a subread index path or a reference.fa path needs to be added \n")
  } else {
    if (!is.null(index_path)) {
      indexPath <- index_path
      if (!file.exists(paste0(indexPath, ".log"))) {
        stop("Genome index does not exist in the specificed location. Please check the full index path again.\n")   
      }
    } else {
      cat("Genome index location not specified. Looking for the index in", output_folder, "\n")
      indexPath <-  file.path(output_folder, "genome_index") 
      if (file.exists(paste0(indexPath, ".log"))) {
        cat("Genome index found in ", output_folder, "...\n")
      } else {
        cat("Genome index not found. Creating one in ", output_folder, " ...\n")
        if(file.exists(ref)){
          Rsubread::buildindex(basename=indexPath, reference=ref)
        } else {
          stop("reference file does not exist. Please check the path and retry. \n")
        }
        
      }
    }
    
    
  }
  
  
  # Generate the output filename
  if (is.null(output_file)) {
    # Only exception is if filename (excluding directory name) contains '.' then will only extract the first part
    fileNameWithoutExtension <- paste0(output_folder, "/", strsplit(basename(R1), "\\.")[[1]][1])
    outbam                   <- paste0(fileNameWithoutExtension, "_aligned.bam")
    cat("Output file name is not provided. Aligned reads are saved in ", outbam, "\n")
  }
  else {
    fileNameWithoutExtension <- paste(output_folder, strsplit(output_file, "\\.")[[1]][1], sep = "/")
    outbam                   <- paste0(output_folder, "/", output_file)
  }
  
  #execute Rsubread align()
  
  if(!is.null(R2) && !all(file.exists(R2))){ # paired-end
    align_output_df <- Rsubread::align(
      index       = indexPath,
      readfile1          = R1,
      readfile2          = R2,
      sortReadsByCoordinates = TRUE,
      output_file = outbam)
  } else {                                   # single-end
    align_output_df <- Rsubread::align(
      index       = indexPath,
      readfile1          = R1,
      sortReadsByCoordinates = TRUE,
      output_file = outbam)
    
  }
  
  utils::write.csv(align_output_df, file = stats_file, row.names = TRUE, quote = FALSE)
  
  #generating the bam index
  Rsamtools::indexBam(paste0(fileNameWithoutExtension, "_aligned.bam"))
  
  # get the unmapped mapped stats to be output and stored in a log file
  bamstats <- Rsamtools::idxstatsBam(paste0(fileNameWithoutExtension, "_aligned.bam"))
  utils::write.csv(bamstats, file = paste0(log_and_stats_folder, "stats_file_align_per_chrom.csv"), row.names = FALSE, quote = FALSE)
  
  cat(
    paste0(
      "sc_aligning finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = TRUE)
  
  return(align_output_df)
}
