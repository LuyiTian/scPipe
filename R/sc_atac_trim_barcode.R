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
#' from a \code{.csv} file (here the barcode is expected to be on the second column). 
#' Currently, for the fastq approach, this can be a list of barcode files.
#' @param valid_barcode_file optional file path of the valid (expected) barcode sequences to be found in the bc_file (.txt, can be txt.gz). 
#' Must contain one barcode per line on the second column separated by a comma (default ="").
#' If given, each barcode from bc_file is matched against the barcode of
#' best fit (allowing a hamming distance of 1). If a FASTQ \code{bc_file} is provided, barcodes with a higher mapping quality, as given by
#' the fastq reads quality score are prioritised.
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
#'
#' @returns None (invisible `NULL`)
#'
#' @examples
#' data.folder <- system.file("extdata", package = "scPipe", mustWork = TRUE)
#' r1      <- file.path(data.folder, "small_chr21_R1.fastq.gz") 
#' r2      <- file.path(data.folder, "small_chr21_R3.fastq.gz") 
#' 
#' # Using a barcode fastq file:
#'
#' # barcodes in fastq format
#' barcode_fastq      <- file.path(data.folder, "small_chr21_R2.fastq.gz") 
#' 
#' sc_atac_trim_barcode (
#' r1            = r1, 
#' r2            = r2, 
#' bc_file       = barcode_fastq,
#' rmN           = TRUE,
#' rmlow         = TRUE,
#' output_folder = tempdir())
#'
#' # Using a barcode csv file:
#'
#' # barcodes in .csv format
#' barcode_1000       <- file.path(data.folder, "chr21_modified_barcode_1000.csv")
#' 
#' \dontrun{
#' sc_atac_trim_barcode (
#' r1            = r1, 
#' r2            = r2, 
#' bc_file       = barcode_1000, 
#' id1_st        = 0,
#` id1_len       = 16,
#` id2_st        = 0,
#` id2_len       = 16 
#' rmN           = TRUE,
#' rmlow         = TRUE,
#' output_folder = tempdir())
#' }
#'@export
sc_atac_trim_barcode <- function(
    r1,
    r2,
    bc_file = NULL,
    valid_barcode_file = "",
    output_folder = "",
    umi_start=0,
    umi_length=0,
    umi_in = "both",
    rmN = FALSE,
    rmlow = FALSE,
    min_qual = 20,
    num_below_min = 2,
    id1_st = -0,
    id1_len = 16,
    id2_st = 0,
    id2_len = 16,
    no_reverse_complement=FALSE) {
    
    if(output_folder == ''){
    output_folder <- file.path(getwd(), "scPipe-atac-output")
    }
    if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    message("Output Directory Does Not Exist. Created Directory: ", output_folder)
    }

    log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
    dir.create(log_and_stats_folder, showWarnings = FALSE)

    log_file             <- paste0(log_and_stats_folder, "log_file.txt")
    stats_file           <- paste0(log_and_stats_folder, "stats_file_trimbarcode.txt")
    if(!file.exists(log_file)) file.create(log_file)
    file.create(stats_file)

    write(
        c("trimbarcode starts at ", as.character(Sys.time()), "\n"), file = log_file, append = TRUE
    )

    if (substr(r1, nchar(r1) - 2, nchar(r1)) == ".gz") {
    	write_gz <- TRUE
    }
    else {
    	write_gz <- FALSE
    }

    if (is.null(bc_file)) {
	stop("Barcode file is mandatory")
	}

	i <- 1;
	for (bc in bc_file) {
		if (!file.exists(bc)) {stop("Barcode file does not exist.")}
		bc_file[i] <- path.expand(bc)
		i          <- i+1;
	}

	if (!file.exists(r1)) {stop("read1 fastq file does not exist.")}
	
	if ((valid_barcode_file != "") && tools::file_ext(valid_barcode_file) != 'csv') {
		stop("Valid Barcode File must be a CSV")
	}

	if(umi_start != 0) {
		if(umi_in %in% c("both", "R1", "R2")) {
			message("UMI Present in: ", umi_in)
		}else{
			stop("Invalid value of umi_in. Possible values are both, R1 and R2")
		}
	}

	# expand tilde to home path for downstream gzopen() call
	r1 <- path.expand(r1)

	if(!is.null(r2)) {
		if (!file.exists(r2)) {stop("read2 file does not exist.")}
		r2 <- path.expand(r2)
	} else{
		r2 <- ""
	}
	message("Saving the output at location: \n", output_folder)

	if(tools::file_ext(bc_file) != "csv"){
		# fastq barcode files are provided, run the FASTQ method
		out_vec <- rcpp_sc_atac_trim_barcode_paired(
			output_folder,
			r1,
			bc_file,
			r2,
			valid_barcode_file,
			write_gz,
			rmN,
			rmlow,
			min_qual,
			num_below_min,
			no_reverse_complement)
	
		write(
			c(
				"Total Reads: ", out_vec[1],
				"\nTotal N's removed: ", out_vec[2],
				"\nremoved_low_qual: ", out_vec[3],
				"\nUnique sequences read in barcode file: ", out_vec[4],
				"\n"
			),
			file = stats_file, append = TRUE)
	
	} else {
		message("Using barcode CSV file, since barcode FastQ file is not passed")
		if(id1_st < 0 || id2_st < 0 || id1_len < 0 || id2_len < 0 ){
			stop("Please pass positive integer values for id1_st, id2_st, id1_len, and id2_len")
		}
		
		
		# trim the barcode csv file (which contains the actual barcodes in the second column)
		# into a file with barcodes on each line and no whitespace
		temp_barcode_file <- paste0(output_folder, "/tempbarcode.csv")
		on.exit(if(file.exists(temp_barcode_file)) {file.remove(temp_barcode_file)})
		
		# change this to handle multiple barcode files!! TODO
		barcodes <- read.csv(bc_file, header=FALSE, strip.white=TRUE)
		write(barcodes$V2, temp_barcode_file)

		# perform the same manipulation for valid
		
		# Check if given barcode start position is valid
		# check_barcode_start_position is expecting a single barcode, of only the barcode sequences, no commas
		message("Checking if id1_st is valid")
		if (!check_barcode_start_position(r1, temp_barcode_file, bc_file, id1_st, id1_len, 10000, .8)) {
			if (tolower(readline(prompt="Continue anyway? (y/n) ")) != "y") {
				stop("Please change id1_st and try again")
			}
			message("Continuing...")
		}
		message("Checking if id2_st is valid")
		if (!check_barcode_start_position(r2, temp_barcode_file, bc_file, id2_st, id2_len, 10000, .8)) {
			if (tolower(readline(prompt="Continue anyway? (y/n) ")) != "y") {
				stop("Please change id2_st and try again")
			}
			message("Continuing...")
		}
		
		out_vec <- rcpp_sc_atac_trim_barcode(
			output_folder,
			r1,
			r2,
			temp_barcode_file, # not urgent but this needs to be changed to bc_file or barcode_file
			valid_barcode_file,
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
		write(
			c(
				"Total Reads: ", out_vec[1],
				"\nTotal N's removed: ", out_vec[2],
				"\nremoved_low_qual: ", out_vec[3],
				"\nExact match Reads: ", out_vec[4],
				"\nReads Matched After Correction: ", out_vec[5],
				"\nTotal barcodes: ", out_vec[6],
				"\n"
			),
			file = stats_file, append = TRUE)
	}
    
    
    write(
        c(
            "trimbarcode finishes at ",
            as.character(Sys.time()),
            "\n\n"
        ),
        file = log_file, append = TRUE)
    
    # return(out_vec)
}