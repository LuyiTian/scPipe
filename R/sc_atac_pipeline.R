
#' @name sc_atac_pipeline
#' @title A convenient function for running the whole pipeline
#' 
#' @param r1 The first read fastq file
#' @param r2 The second read fastq file
#' @param barcode_fastq The barcode fastq file
#' @param organism The name of the organism e.g. hg38
#' @param reference The reference genome file
#' @param feature_type The feature type (either `genome_bin` or `peak`)
#' @param remove_duplicates Whether or not to remove duplicates (samtools is required) 
#' @param samtools_path A custom path of samtools to use for duplicate removal
#' @param cell_calling The desired cell calling method (either `emptydrops`, `filter`, or `cellranger`)
#' @param output_folder The path of the output folder
#'
#' @export
#' 
sc_atac_pipeline <- function(r1,
                             r2,
                             barcode_fastq,
                             organism,
                             reference,
                             feature_type,
                             remove_duplicates = FALSE,
                             samtools_path = NULL,
                             cell_calling = FALSE,
                             output_folder = NULL) {
  
  get_filename_without_extension <- function(path, extension_length = 1) {
    vec <- strsplit(basename(path), "\\.")[[1]]
    name.size <- length(vec) - extension_length
    return(paste(vec[1:name.size], collapse = "."))
  }

  r1_name <- get_filename_without_extension(r1, extension_length = 2)
  r2_name <- get_filename_without_extension(r2, extension_length = 2)

  sc_atac_trim_barcode (r1            = r1,
                        r2            = r2,
                        bc_file       = barcode_fastq,
                        rmN           = TRUE,
                        rmlow         = TRUE,
                        output_folder = output_folder)

  demux_r1        <- file.path(output_folder, paste0("demux_", r1_name, ".fastq.gz"))
  demux_r2        <- file.path(output_folder, paste0("demux_", r2_name, ".fastq.gz"))

  reference       <- reference
  
  sc_aligning(ref = reference,
              R1 = demux_r1,
              R2 = demux_r2,
              nthreads  = 12,
              output_folder = output_folder)

  bam_to_tag  <- file.path(output_folder, paste0("demux_", r1_name, "_aligned.bam"))

  sc_atac_bam_tagging(inbam         = bam_to_tag,
                      output_folder = output_folder,
                      bam_tags      = list(bc="CB", mb="OX"),
                      nthreads      =  12)

  sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted.bam"))

  if (isTRUE(remove_duplicates)) {
    sc_atac_remove_duplicates(inbam = sorted_tagged_bam,
                              samtools_path = samtools_path,
                              output_folder = output_folder)
    sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted_markdup.bam"))
  }

  sc_atac_create_fragments(inbam = sorted_tagged_bam,
                           output_folder = output_folder)

  features <- NULL
  if (feature_type == "peak") {
    sc_atac_peak_calling(inbam = sorted_tagged_bam,
                         output_folder = output_folder,
                         reference)

    features <- file.path(output_folder, "NA_peaks.narrowPeak")
  } else { # otherwise is genome bins
    features <- reference
  }

  sc_atac_feature_counting (insortedbam   = sorted_tagged_bam,
                            feature_input = features,
                            bam_tags      = list(bc="CB", mb="OX"),
                            feature_type  = feature_type,
                            organism      = organism,
                            cell_calling  = cell_calling,
                            genome_size   = NULL,
                            bin_size      = NULL,
                            yieldsize     = 1000000,
                            mapq          = 30,
                            exclude_regions = TRUE,
                            output_folder = output_folder,
                            fix_chr       = "none")
  
  sce <- sc_atac_create_sce(input_folder = output_folder,
                            organism     = "hg38",
                            feature_type = "peak",
                            pheno_data   = NULL,
                            report       = TRUE)
  return(sce)
}

#' @name sc_atac_pipeline_quick_test
#' @title A quick test for running the pipeline
#'
sc_atac_pipeline_quick_test <- function() {
  data.folder <- system.file("extdata", package = "scPipe", mustWork = TRUE)
  out <- tryCatch(
    {
      sce <- sc_atac_pipeline(r1 = file.path(data.folder, "small_R1.fastq.gz"),
                              r2 = file.path(data.folder, "small_R3.fastq.gz"),
                              barcode_fastq = file.path(data.folder, "small_R2.fastq.gz"),
                              organism = "hg38",
                              reference = file.path(data.folder, "genome.fa"),
                              remove_duplicates = FALSE,
                              feature_type = "peak",
                              cell_calling = "filter",
                              output_folder = file.path(tempdir(), "scPipe-atac-output"))
      cat("Successfully ran pipeline.\n")
    },
    finally = {
      system2("rm", c("-rf", file.path(data.folder, "scPipe-atac-output")))
    }
  )
}
