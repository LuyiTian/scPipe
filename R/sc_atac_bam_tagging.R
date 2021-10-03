#############################
# Demultiplxing FASTQ Reads
#############################

#' sc_atac_bam_tagging()
#' @name sc_atac_bam_tagging
#' @title BAM tagging
#' @description Demultiplexes the reads
#' 
#' @param inbam The input BAM file
#' @param output_folder The path of the output folder
#' @param bc_length The length of the cellular barcodes
#' @param bam_tags The BAM tags
#' @param nthreads The number of threads
#' 
#' @export
#' 
sc_atac_bam_tagging <- function(inbam,
                                output_folder = NULL,
                                bc_length = NULL,
                                bam_tags = list(bc="CB", mb="OX"),
                                nthreads = 1) {

  if (any(!file.exists(inbam))) {
    stop("At least one input bam file should be present")
  } else {
    inbam <- path.expand(inbam)
  }

  if(is.null(output_folder)){
    output_folder <- file.path(getwd(), "scPipe-atac-output")
    # output_folder <- paste0(sub('[/][^/]+$', '', inbam)) # same as output_folder <- basename(inbam)
  }

  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory does not exist. Creating ", output_folder, "\n")
  }

  get_filename_without_extension <- function(path) {
    vec <- strsplit(basename(path), "\\.")[[1]]
    return(paste(vec[1:length(vec)-1], collapse = "."))
  }
  fileNameWithoutExtension <- get_filename_without_extension(inbam)
  
  outbam                   <- file.path(output_folder, paste0(fileNameWithoutExtension, "_tagged_sorted.bam"))
  outsortedbam             <- file.path(output_folder, paste0(fileNameWithoutExtension, "_tagged_sorted")) # don't want extension

  log_and_stats_folder       <- file.path(output_folder, "scPipe_atac_stats")
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file                   <- file.path(log_and_stats_folder, "log_file.txt")
  stats_file                 <- file.path(log_and_stats_folder, "stats_file_bam_tagging.txt")
  if(!file.exists(log_file)) file.create(log_file)

  cat(
    paste0(
      "sc_atac_tagging starts at ",
      as.character(Sys.time()),
      "\n"
    ),
    file = log_file, append = TRUE)

  outbam <- path.expand(outbam)
  if(!file.exists(outbam)){
    file.create(outbam)

  rcpp_sc_atac_bam_tagging(inbam, outbam, bam_tags$bc, bam_tags$mb, nthreads)

  cat("Tagged BAM file is located in: \n")
  cat(outbam)
  cat("\n")

  # Rsamtools::sortBam(outbam, outsortedbam, indexDestination = TRUE, maxMemory = 1024/32*max_memory)
  Rsamtools::indexBam(paste0(outsortedbam, ".bam"))

  cat("Sorted & indexed tagged BAM file is located in: \n")
  }
  cat(paste0(outsortedbam, ".bam"))
  cat("\n")


  if(is.null(bc_length)){
    cat("Using default value for barcode length (bc_length = 16) \n")
    bc_length <- 16
  }


  flag_defs <- tibble::tibble(
    type =
      c("one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "both_reads_unmapped", "both_reads_unmapped", "mapped", "mapped", "mapped", "mapped", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously")
    ,
    flag =
      c(73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77, 141, 99, 147, 83, 163, 67, 131, 115, 179, 81, 161, 97, 145, 65, 129, 113, 177))


  bam0 <- Rsamtools::scanBam(inbam)
  barcodes <- substr(bam0[[1]]$qname, 1, bc_length)

  barcode_info <- tibble::tibble(
    barcode = barcodes,
    flag    = bam0[[1]]$flag) %>%
    dplyr::left_join(flag_defs, by = "flag") %>%
    dplyr::left_join(as.data.frame(table(barcodes)) %>%
                purrr::set_names(c("barcode", "number_of_reads")),
              by = "barcode")

  barcode_stats_filename <- file.path(log_and_stats_folder, "mapping_stats_per_barcode.csv")
  cat("Saving csv file with barcode stats in", barcode_stats_filename, "\n")
  utils::write.csv(barcode_info, barcode_stats_filename, row.names = FALSE)

  cat(
    paste0(
      "sc_atac_tagging finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ),
    file = log_file, append = TRUE)
}