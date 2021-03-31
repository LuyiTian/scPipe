#############################
# Demultiplxing FASTQ Reads
#############################

#' sc_atac_bam_tagging()
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

sc_atac_bam_tagging <- function(inbam, 
                                output_folder = "",
                                bc_length = NULL,
                                bam_tags = list(bc="CB", mb="OX"),
                                nthreads = 1
) {
  
  if (any(!file.exists(inbam))) {
    stop("At least one input bam file should be present")
  } else {
    cat("sc_atac_bam_tagging1\n")
    inbam = path.expand(inbam)
  }
  
  if(output_folder == ''){
    #output_folder <- file.path(getwd(), "scPipe-atac-output")
    output_folder <- paste0(sub('[/][^/]+$', '', inbam))
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory does not exist. Created at: ", output_folder, "\n")
  }
  
  #if(output_folder == ""){
    # fileNameWithoutExtension <- strsplit(basename(inbam), "\\.")[[1]][1]
    # outbam                   <- paste(fileNameWithoutExtension, "_tagged.bam", sep = "")
    # outsortedbam             <- paste(fileNameWithoutExtension, "_tagged_sorted", sep = "")
  #} else{
    # if(!dir.exists(output_folder)){
    #   cat(output_folder, "does not exist.\nCreating folder...")
    #   dir.create(output_folder)
    #   cat("Created.\n")
    # }
    fileNameWithoutExtension <- strsplit(basename(inbam), "\\.")[[1]][1]
    outbam                   <- paste(output_folder, "/", fileNameWithoutExtension, "_tagged.bam", sep = "")
    outsortedbam             <- paste(output_folder, "/", fileNameWithoutExtension, "_tagged_sorted", sep = "")
  #}
  
  
  log_and_stats_folder       <- paste0(output_folder, "/scPipe_atac_stats/")
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file                   <- paste0(log_and_stats_folder, "log_file.txt")
  stats_file                 <- paste0(log_and_stats_folder, "stats_file_bam_tagging.txt")
  if(!file.exists(log_file)) file.create(log_file)
  # file.create(stats_file)
  
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
  }
  
  # if(output_folder == ''){
  #   output_folder <- file.path(getwd(), "scPipe-atac-output")
  # }
  # 
  # if (!dir.exists(output_folder)){
  #   dir.create(output_folder,recursive=TRUE)
  #   cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  # }
  
  
  rcpp_sc_atac_bam_tagging(inbam, outbam, bam_tags$bc, bam_tags$mb,nthreads)
  
  cat("Tagged BAM file is located in: \n")
  cat(outbam)
  cat("\n")
  
  Rsamtools::sortBam(outbam, outsortedbam, indexDestination = TRUE)
  Rsamtools::indexBam(paste0(outsortedbam, ".bam"))
  
  cat("Sorted & indexed tagged BAM file is located in: \n")
  cat(paste0(outsortedbam, ".bam"))
  cat("\n")
  
  
  if(is.null(bc_length)){
    cat("Using default value for barcode length (bc_length = 16) \n")
    bc_length <- 16
  }
  
  
  flag_defs = tibble(
    type = 
      c("one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "both_reads_unmapped", "both_reads_unmapped", "mapped", "mapped", "mapped", "mapped", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously")
    ,
    flag = 
      c(73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77, 141, 99, 147, 83, 163, 67, 131, 115, 179, 81, 161, 97, 145, 65, 129, 113, 177))
  
  
  bam0 = scanBam(inbam)
  barcodes = substr(bam0[[1]]$qname, 1, bc_length)
  
  barcode_info <- tibble(
    barcode = barcodes,
    flag    = bam0[[1]]$flag) %>%
    left_join(flag_defs, by = "flag") %>%
    left_join(as.data.frame(table(barcodes)) %>% 
                purrr::set_names(c("barcode", "number_of_reads")),
              by = "barcode")
  
  barcode_stats_filename <- paste(log_and_stats_folder, "mapping_stats_per_barcode.csv", sep = "")
  cat("Saving csv file with barcode stats in", barcode_stats_filename)
  write.csv(barcode_info, barcode_stats_filename, row.names = FALSE)
  
  
  # generate the fragment file for the BAM file 
  # need bedtools v2.26.0 or later
  # system2("bedtools", c("bamToBed", "i", outsortedbam), "|", "awk", c(-F"#" '{print $1"\t"$2}'), stdout = paste(output_folder,"/fragments.bed",sep = ""))
  
  cat(
    paste0(
      "sc_atac_tagging finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = TRUE)
  
}
