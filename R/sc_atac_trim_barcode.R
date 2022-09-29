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
#' @param barcodeFastq optional list of FASTQ files containing barcode sequences as sequences (from 10x-ATAC). Each FASTQ must be of the same length
#' as r1 and r2.
#' @param valid_barcode_file optional CSV file path of the valid (expected) barcode sequences to be found in either the reads or supplied \code{barcodeFastq}.
#' Must contain one barcode per line on the second column (comma separated)
#' @param output_folder the output dir for the demultiplexed fastq file, which will contain 
#' fastq files with reformatted barcode and UMI into the read name. 
#' Files ending in \code{.gz} will be automatically compressed.
#' @param id1_st barcode start position (0-indexed) for read 1, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param id2_st barcode start position (0-indexed) for read 2, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param id1_len barcode length for read 1, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param id2_len barcode length for read 2, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param umi_start if available, the start position of the molecular identifier.
#' @param umi_length if available, the start position of the molecular identifier.  
#' @param umi_in umi_in
#' @param rmN logical, whether to remove reads that contains N in UMI or cell barcode.
#' @param rmlow logical, whether to remove reads that have low quality barcode sequences
#' @param min_qual the minimum base pair quality that is allowed (default = 20).
#' @param num_below_min the maximum number of base pairs below the quality threshold.
#' @param no_reverse_complement specifies if the reverse complement of the barcode sequence should be 
#' used for barcode error correction (only when barcode sequences are provided as fastq files). FALSE (default)
#' lets the function decide whether to use reverse complement, and TRUE forces the function to
#' use the forward barcode sequences.
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
#' bc_file       = barcode_csv, 
#' id1_st        = 0,
#` id1_len       = 16,
#` id2_st        = 0,
#` id2_len       = 16 
#' rmN           = TRUE,
#' rmlow         = TRUE,
#' output_folder = "")
#' }
#'@export

sc_atac_trim_barcode <- function(
	r1,
	r2,
	barcodeFastq = NULL,
	valid_barcode_file = "",
	output_folder = "",
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
	id2_len = -10,
	no_reverse_complement=FALSE) {
  
	if(output_folder == ''){
	output_folder <- file.path(getwd(), "scPipe-atac-output")
	}
	if (!dir.exists(output_folder)){
	dir.create(output_folder,recursive=TRUE)
	cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
	}

	log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
	dir.create(log_and_stats_folder, showWarnings = FALSE)

	log_file             <- paste0(log_and_stats_folder, "log_file.txt")
	stats_file           <- paste0(log_and_stats_folder, "stats_file_trimbarcode.txt")
	if(!file.exists(log_file)) file.create(log_file)
	file.create(stats_file)

	cat(
	paste0( "trimbarcode starts at ", as.character(Sys.time()),"\n"), file = log_file, append = TRUE)

	if (substr(r1, nchar(r1) - 2, nchar(r1)) == ".gz") {
	write_gz <- TRUE
	}
	else {
	write_gz <- FALSE
	}

	if (!is.null(barcodeFastq)) {
		i=1;
		for (bc in barcodeFastq) {
			if (!file.exists(bc)) {stop("Barcode file does not exist.")}
			barcodeFastq[i] <- path.expand(bc)
			i          <- i+1;
		}
	}

	if (!file.exists(r1)) {stop("read1 fastq file does not exist.")}
	
	if ((valid_barcode_file != "") && file_ext(valid_barcode_file) != 'csv') {
		stop("Valid Barcode File must be a CSV")
	}

	if(umi_start != 0){
		if(umi_in %in% c("both", "R1", "R2")){
		cat("UMI Present in: ", umi_in, "\n")
		}else{
		stop("Invalid value of umi_in. Possible values are both, R1 and R2")
		}
	}

	# expand tilde to home path for downstream gzopen() call
	r1 <- path.expand(r1)

	if(!is.null(r2)){
		if (!file.exists(r2)) {stop("read2 file does not exist.")}
		r2 <- path.expand(r2)
	}else{
		r2 <- ""
	}
	cat("Saving the output at location: ")
	cat(output_folder)
	cat("\n")

	if(!is.null(barcodeFastq)){
		# fastq barcode files are provided, run the FASTQ method
		out_vec <- rcpp_sc_atac_trim_barcode_paired(
		output_folder,
		r1,
		barcodeFastq,
		r2,
		valid_barcode_file,
		write_gz,
		rmN,
		rmlow,
		min_qual,
		num_below_min,
		no_reverse_complement)
		
		cat("Total Reads: ", out_vec[1],
			"\nTotal N's removed: ", out_vec[2],
			"\nremoved_low_qual: ", out_vec[3],
			"\nUnique sequences read in barcode file: ", out_vec[4],
			"\n",
			file = stats_file, append = TRUE)
		
	} else {
		cat("Using barcode CSV file, since barcode FastQ files are not passed \n")
		if(id1_st < 0 || id2_st < 0 || id1_len < 0 || id2_len < 0 ) {
			stop("Please pass positive integer values for id1_st, id2_st, id1_len, and id2_len")
		}
		
		if (valid_barcode_file != "") {
			# trim the barcode csv file (which contains the actual barcodes in the second column)
			# into a file with barcodes on each line and no whitespace
			temp_barcode_file <- paste0(output_folder, "/tempbarcode.csv")
			on.exit(if(file.exists(temp_barcode_file)) {file.remove(temp_barcode_file)})
			
			# change this to handle multiple barcode files!! TODO
			barcodes <- read.csv(valid_barcode_file, header=FALSE, strip.white=TRUE)
			write(barcodes$V2, temp_barcode_file)
			
			# Check if given barcode start position is valid
			# check_barcode_start_position is expecting a single barcode, of only the barcode sequences, no commas
			cat("Checking if id1_st is valid\n")
			if (!check_barcode_start_position(r1, temp_barcode_file, valid_barcode_file, id1_st, id1_len, 10000, .8)) {
				if (tolower(readline(prompt="Continue anyway? (y/n) ")) != "y") {
					stop("Please change id1_st and try again")
				}
				cat("Continuing...")
			}
			cat("Checking if id2_st is valid\n")
			if (!check_barcode_start_position(r2, temp_barcode_file, valid_barcode_file, id2_st, id2_len, 10000, .8)) {
				if (tolower(readline(prompt="Continue anyway? (y/n) ")) != "y") {
					stop("Please change id2_st and try again")
				}
				cat("Continuing...")
			}
		} else {
			temp_barcode_file = ""
		}
		
		out_vec <- rcpp_sc_atac_trim_barcode(
			output_folder,
			r1,
			r2,
			temp_barcode_file,
			#bc_start,
			#bc_length,
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
		
		# concatenate results to stats_file
		cat("Total Reads: ", out_vec[1],
			"\nTotal N's removed: ", out_vec[2],
			"\nremoved_low_qual: ", out_vec[3],
			"\nExact match Reads: ", out_vec[4],
			"\nReads Matched After Correction: ", out_vec[5],
			"\nTotal barcodes: ", out_vec[6],
			"\n",
			file = stats_file, append = TRUE)
		
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