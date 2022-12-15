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
#' @param type type of sequencing data (`RNA` or `DNA`)
#' @param input_format a string indicating the input format
#' @param output_format a string indicating the output format
#' @param nthreads numeric value giving the number of threads used for mapping. 
#'
#' @returns the file path of the output aligned BAM file
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
    type          = "dna",
    nthreads      = 1){

    if(!all(file.exists(R1))){
        stop("At least one of the input files for R1 does not exist")    
    }
    
    if(!is.null(R2) && !all(file.exists(R2))){
        stop("At least one of input file for R2 does not exist")
    }
    
    if(tech == "atac") {
        message("ATAC-Seq mode is selected...")
    
        if(is.null(output_folder)) {
            output_folder        <- file.path(getwd(), "scPipe-atac-output")
        }
        
        if (!dir.exists(output_folder)){
            dir.create(output_folder,recursive=TRUE)
            message("Output directory is not provided. Created directory: ", output_folder)
        }
        
        log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
        type                 <- "dna"
        
    } else if(tech == "rna") {
        message("RNA-Seq mode is selected...")
        
        if(is.null(output_folder)) {
            stop("output_folder cannot be NULL for rna mode. Aborting...\n")
        }
        log_and_stats_folder <- output_folder
        type                 <- "rna"  
    } #else

    dir.create(log_and_stats_folder, showWarnings = FALSE)
    log_file             <- paste0(log_and_stats_folder, "log_file.txt")
    stats_file           <- paste0(log_and_stats_folder, "stats_file_align.txt")
    if(!file.exists(log_file)) file.create(log_file)

    write(
        c(
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
            message("Genome index location not specified. Looking for the index in", output_folder)
            indexPath <-  file.path(output_folder, "genome_index") 
            if (file.exists(paste0(indexPath, ".log"))) {
                message("Genome index found in ", output_folder, "...")
            } else {
                message("Genome index not found. Creating one in ", output_folder, " ...")
                if(file.exists(ref)){
                    Rsubread::buildindex(basename=indexPath, reference=ref)
                } else {
                    stop("reference file does not exist. Please check the path and retry. \n")
                }
                
            }
        }
    }

    # Check for partial/nomatch files
    if(tech == "atac") {
        containing_folder <- dirname(R1) # Assume partial and nomatch files are also in the same directory as supplied input fastq files
        input_folder_files <- list.files(containing_folder)

        # Initialise demultiplexing stats
        barcode_completematch_count <- length(readLines(R1))/2
        demux_stats <- data.frame(status = c("barcode_completematch_count"),
                                count = c(barcode_completematch_count))
        # Concatenate the complete and partial matches

        partial_matches_R1 <- file.path(containing_folder, input_folder_files[grep("dem.+partialmatch.+R1.+fastq", input_folder_files)])
        partial_matches_R3 <- file.path(containing_folder, input_folder_files[grep("dem.+partialmatch.+R3.+fastq", input_folder_files)])

        if (all(file.exists(partial_matches_R1, partial_matches_R3)) && !identical(partial_matches_R1, character(0)) && !identical(partial_matches_R3, character(0))) {
            if (length(readLines(partial_matches_R1)) > 0 && length(readLines(partial_matches_R3)) > 0) {
                message("Found partial match fastq files, proceeding to concatenate with complete match fastq files.")
                barcode_partialmatch_count <- length(readLines(partial_matches_R1))/2
                demux_stats <- demux_stats %>% tibble::add_row(status = "barcode_partialmatch_count", count = barcode_partialmatch_count)
                concat_filename_R1 <- paste0("demultiplexed_complete_partialmatch_", stringr::str_remove(basename(R1), stringr::regex("dem.+completematch_")))
                concat_file_R1 <- file.path(containing_folder, concat_filename_R1)
                concat_filename_R3 <- paste0("demultiplexed_complete_partialmatch_", stringr::str_remove(basename(R2), stringr::regex("dem.+completematch_")))
                concat_file_R3 <- file.path(containing_folder, concat_filename_R3)
                system2("zcat", c(R1, partial_matches_R1, "|", "gzip", "-c", ">", concat_file_R1))
                system2("zcat", c(R2, partial_matches_R3, "|", "gzip", "-c", ">", concat_file_R3))

                if (!all(file.exists(concat_file_R1, concat_file_R3))) {
                    stop("Couldn't concatenate files!\n")
                }

                message("Output concatenated read files to:")
                message("R1:", concat_file_R1)
                message("R3:", concat_file_R3)

                # Replace original fastq files with concatenated files for aligning
                R1 <- concat_file_R1
                R2 <- concat_file_R3

            } else {
                message("No partial matches, checking for reads with non-matched barcodes.")
            }

        }
        # ------------ Align the nomatch file -------

        no_matches_R1 <- file.path(containing_folder, input_folder_files[grep("nomatch.+R1.+fastq", input_folder_files)])
        no_matches_R3 <- file.path(containing_folder, input_folder_files[grep("nomatch.+R3.+fastq", input_folder_files)])
        if (all(file.exists(no_matches_R1, no_matches_R3)) && !identical(no_matches_R1, character(0)) && !identical(no_matches_R3, character(0))) {
            if (length(readLines(no_matches_R1)) > 0 && length(readLines(no_matches_R3)) > 0) {
                message("Found barcode non-matched demultiplexed FASTQ files. Proceeding to align them.")

                fileNameWithoutExtension <- paste0(output_folder, "/", strsplit(basename(no_matches_R1), "\\.")[[1]][1])
                nomatch_bam <- paste0(fileNameWithoutExtension, "_aligned.bam")
                Rsubread::align(
                index = file.path(output_folder, "genome_index"),
                readfile1 = no_matches_R1,
                readfile2 = no_matches_R3,
                sortReadsByCoordinates = TRUE,
                type = "DNA",
                nthreads = 12,
                output_file = nomatch_bam)

                # Extract columns
                bam_tags <-list(bc="CB", mb="OX")
                param <- Rsamtools::ScanBamParam(tag = as.character(bam_tags),  mapqFilter=20)
                bamfl <- open(Rsamtools::BamFile(nomatch_bam))
                params <- Rsamtools::ScanBamParam(what=c("flag"), tag=c("CB"))
                bam0 <- Rsamtools::scanBam(bamfl, param = params)

                flag_defs <- tibble::tibble(
                type =
                    paste0("barcode_unmatch_", c("one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "both_reads_unmapped", "both_reads_unmapped", "mapped", "mapped", "mapped", "mapped", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously"))
                ,
                flag =
                    c(73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77, 141, 99, 147, 83, 163, 67, 131, 115, 179, 81, 161, 97, 145, 65, 129, 113, 177))

                # Create stats data frame
                demux_stats <- rbind(demux_stats, as.data.frame(table((data.frame(flag = bam0[[1]]$flag) %>% dplyr::left_join(flag_defs, by = "flag"))[,c('type')])) %>%
                                    dplyr::rename(status = Var1, count = Freq))

            } else {
                message("No reads found with non-matching barcodes.")
            }
        }
        utils::write.csv(demux_stats, file.path(log_and_stats_folder, "demultiplexing_stats.csv"), row.names = FALSE)
        message("Outputted demultiplexing stats file to", file.path(log_and_stats_folder, "demultiplexing_stats.csv"), "\n")

    } 

    # Generate the output filename
    if (is.null(output_file)) {
        # Only exception is if filename (excluding directory name) contains '.' then will only extract the first part
        fileNameWithoutExtension <- paste0(output_folder, "/", strsplit(basename(R1), "\\.")[[1]][1])
        outbam                   <- paste0(fileNameWithoutExtension, "_aligned.bam")
        message("Output file name is not provided. Aligned reads are saved in ", outbam)
    }
    else {
        fileNameWithoutExtension <- paste(output_folder, strsplit(output_file, "\\.")[[1]][1], sep = "/")
        outbam                   <- paste0(output_folder, "/", output_file)
    }
    
    #execute Rsubread align()
    
    if(!is.null(R2) && file.exists(R2)){       # paired-end
        align_output_df <- Rsubread::align(
        index       = indexPath,
        readfile1          = R1,
        readfile2          = R2,
        sortReadsByCoordinates = TRUE,
        type = type,
        nthreads = nthreads,
        output_file = outbam)
    } else {                                   # single-end
        align_output_df <- Rsubread::align(
        index       = indexPath,
        readfile1          = R1,
        sortReadsByCoordinates = TRUE,
        type = type,
        nthreads = nthreads,
        output_file = outbam)
        
    }
    
    utils::write.csv(align_output_df, file = stats_file, row.names = TRUE, quote = FALSE)
    
    #generating the bam index
    Rsamtools::indexBam(paste0(fileNameWithoutExtension, "_aligned.bam"))
    
    # get the unmapped mapped stats to be output and stored in a log file
    bamstats <- Rsamtools::idxstatsBam(paste0(fileNameWithoutExtension, "_aligned.bam"))
    utils::write.csv(bamstats, file = paste0(log_and_stats_folder, "stats_file_align_per_chrom.csv"), row.names = FALSE, quote = FALSE)

    write(
        c(
        "sc_aligning finishes at ",
        as.character(Sys.time()),
        "\n\n"
        ), 
        file = log_file, append = TRUE)
    
    return(outbam)
}
