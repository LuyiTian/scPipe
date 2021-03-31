###############################
# Demultiplexing FASTQ Reads
###############################

#' @name sc_atac_trim_barcode
#' @title demultiplex raw single-cell ATAC-Seq fastq reads
#' @description  single-cell data need to be demultiplexed in order to retain the information of the cell barcodes
#' the data belong to. Here we reformat fastq files so barcode/s (and if available the UMI sequences) are moved from
#' the sequence into the read name. Since scATAC-Seq data are mostly paired-end, both `r1` and `r2` are demultiplexed in this function.
#' @param r1 read one for pair-end reads.
#' @param r2 read two for pair-end reads, NULL if single read.
#' @param bc_file the barcode information, can be either in a \code{fastq} format (e.g. from 10x-ATAC) or
#' from a \core{.csv} file (here the barcode is expected to be on the second column). 
#' Currently, for the fastq approach, this can be a list of barcode files.
#' @param output_folder the output dir for the demultiplexed fastq file, which will contain 
#' fastq files with reformatted barcode and UMI into the read name. 
#' Files ending in \code{.gz} will be automatically compressed.
#' @param bc_start barcode start position (0-indexed), which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param bc_length barcode length, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param umi_start if available, the start position of the molecular identifier.
#' @param umi_length if available, the start position of the molecular identifier.  
#' @param rmN logical, whether to remove reads that contains N in UMI or cell barcode.
#' @param rmlow logical, whether to remove the low quality reads.
#' @param min_qual the minimum base pair quality that is allowed.
#' @param num_below_min the maximum number of base pairs below the quality threshold.
#' @examples
#' \dontrun{
#' using a barcode fastq file
#' sc_atac_trim_barcode (
#' r1            = r1, 
#' r2            = r2, 
#' bc_file       = barcode_fastq,
#' rmN           = TRUE,
#' rmlow         = TRUE,
#' output_folder = "")
#' 
#' using a barcode csv file
#' sc_atac_trim_barcode (
#' r1            = r1, 
#' r2            = r2, 
#' bc_file       = barcode_1000, 
#' bc_start      = 3, 
#' bc_length     = 16, 
#' rmN           = TRUE,
#' rmlow         = TRUE,
#' output_folder = "")
#' }
#'@export

sc_atac_trim_barcode = function(
  r1,
  r2,
  bc_file,
  output_folder = "",
  bc_start=-1,
  bc_length=-1,
  umi_start=0,
  umi_length=0,
  umi_in = "both",
  rmN = FALSE,
  rmlow = FALSE,
  min_qual = 20,
  num_below_min = 2,
  id1_st = -1,
  id1_len = -1,
  id2_st = -1,
  id2_len = -10) {

  if(output_folder == ''){
    output_folder = file.path(getwd(), "scPipe-atac-output")
  }
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  }
 
  log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
  dir.create(log_and_stats_folder, showWarnings = F)
  
  log_file             <- paste0(log_and_stats_folder, "log_file.txt")
  stats_file           <- paste0(log_and_stats_folder, "stats_file_trimbarcode.txt")
  if(!file.exists(log_file)) file.create(log_file)
  file.create(stats_file)

  cat(
    paste0( "trimbarcode starts at ", as.character(Sys.time()),"\n"), file = log_file, append = TRUE)
  
  if (substr(r1, nchar(r1) - 2, nchar(r1)) == ".gz") {
    write_gz = TRUE
  }
  else {
    write_gz = FALSE
  }

  if (!is.null(bc_file)) {
    if (!file.exists(r1)) {stop("read1 fastq file does not exist.")}
    i=1;
    for (bc in bc_file) {
      if (!file.exists(bc)) {stop("Barcode file does not exist.")}
      bc_file[i] = path.expand(bc)
      i = i+1;
    }

    if(umi_start != 0){
      if(umi_in %in% c("both", "R1", "R2")){
        cat("UMI Present in: ", umi_in, "\n")
      }else{
        stop("Invalid value of umi_in. Possible values are both, R1 and R2")
      }
    }

    # expand tilde to home path for downstream gzopen() call
    r1 = path.expand(r1)

    if(!is.null(r2)){
      if (!file.exists(r2)) {stop("read2 file does not exist.")}
      r2 = path.expand(r2)
    }else{
      r2=""
    }
    cat("Saving the output at location: ")
    cat(output_folder)
    cat("\n")

    if(file_ext(bc_file) != 'csv'){
      out_vec = rcpp_sc_atac_trim_barcode_paired(
        output_folder,
        r1,
        bc_file,
        r2,write_gz,
        rmN,
        rmlow,
        min_qual,
        num_below_min,
        id1_st,
        id1_len,
        id2_st,
        id2_len,
        umi_start,
        umi_length)

      cat("Total Reads: ", out_vec[1], 
          "\nTotal N's removed: ", out_vec[2], 
          "\nremoved_low_qual: ", out_vec[3], 
          "\nUnique sequences read in barcode file: ", out_vec[4],
          "\n",
          file = stats_file, append = TRUE)
      
    }
    else {
      cat("Using barcode CSV file, since barcode FastQ file is not passed \n")
      if(bc_start == -1 || bc_length == -1 ){
        stop("Please pass bc_start and bc_length values")
      }

      # Check if given barcode start position is valid
      if (!check_barcode_start_position(r1, bc_file, bc_start, bc_length, 10000, .8)) {
          stop("Please change bc_start and try again")
      }

      out_vec = rcpp_sc_atac_trim_barcode(
        output_folder,
        r1,
        r2,
        bc_file,
        bc_start,
        bc_length,
        umi_start,
        umi_length,
        umi_in,
        write_gz,
        rmN,
        rmlow,
        min_qual,
        num_below_min,
        id1_st,
        id1_len,
        id2_st,
        id2_len)

      bc <- data.table::fread("../new_project_collab/seqATAC/data/barcode.csv", select = 2, col.names = "bc")
      
      cat("Total Reads: ", out_vec[1], 
          "\nTotal N's removed: ", out_vec[2], 
          "\nremoved_low_qual: ", out_vec[3], 
          "\nExact match Reads: ", out_vec[4], 
          "\nApprox Match Reads: ", out_vec[5], 
          "\nTotal barcodes: ", length(unique(bc$bc)),
          "\n",
          file = stats_file, append = TRUE)
      
    }
  }else{
    stop("Barcode file is mandatory")
  }


  cat(
    paste0(
      "trimbarcode finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = TRUE)
  
  # return(out_vec)
}



