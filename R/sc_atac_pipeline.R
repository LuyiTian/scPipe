#' @name sc_atac_pipeline
#' 
#' @title A convenient function for running the entire pipeline
#' 
#' @param r1 The first read fastq file
#' @param r2 The second read fastq file
#' @param barcode_fastq The barcode fastq file (need either this or `barcode_csv`)
#' @param barcode_csv The barcode csv file (need either this or `barcode_fastq`)
#' @param organism The name of the organism e.g. hg38
#' @param reference The reference genome file
#' @param feature_type The feature type (either `genome_bin` or `peak`)
#' @param remove_duplicates Whether or not to remove duplicates (samtools is required) 
#' @param samtools_path A custom path of samtools to use for duplicate removal
#' @param bin_size The size of the bins for feature counting with the `genome_bin` feature type
#' @param yieldsize The number of reads to read in for feature counting
#' @param mapq The minimum MAPQ score
#' @param exclude_regions Whether or not the regions should be excluded
#' @param excluded_regions_filename The filename of the file containing the regions to be excluded
#' @param cell_calling The desired cell calling method either \code{cellranger}, \code{emptydrops} or  \code{filter}
#' @param promoters_file The path of the promoter annotation file (if the specified organism isn't recognised)
#' @param tss_file The path of the tss annotation file (if the specified organism isn't recognised)
#' @param enhs_file The path of the enhs annotation file (if the specified organism isn't recognised)
#' @param gene_anno_file The path of the gene annotation file (gtf or gff3 format)
#' @param fix_chr Specify `none`, `exclude_regions`, `feature` or `both` to prepend the string "chr" to the start of the associated file
#' @param lower the lower threshold for the data if using the \code{emptydrops} function for cell calling.
#' @param genome_size The size of the genome (used for the \code{cellranger} cell calling method)
#' @param min_uniq_frags The minimum number of required unique fragments required for a cell (used for \code{filter} cell calling)
#' @param max_uniq_frags The maximum number of required unique fragments required for a cell (used for \code{filter} cell calling)
#' @param min_frac_peak The minimum proportion of fragments in a cell to overlap with a peak (used for \code{filter} cell calling)
#' @param min_frac_tss The minimum proportion of fragments in a cell to overlap with a tss (used for \code{filter} cell calling)
#' @param min_frac_enhancer The minimum proportion of fragments in a cell to overlap with a enhancer sequence (used for \code{filter} cell calling)
#' @param min_frac_promoter The minimum proportion of fragments in a cell to overlap with a promoter sequence (used for \code{filter} cell calling)
#' @param max_frac_mito The maximum proportion of fragments in a cell that are mitochondrial (used for \code{filter} cell calling)
#' @param report Whether or not a HTML report should be produced
#' @param nthreads The number of threads to use for alignment (sc_align) and demultiplexing (sc_atac_bam_tagging)
#' @param output_folder The path of the output folder
#'
#' @export
#' 
sc_atac_pipeline <- function(r1,
                             r2,
                             barcode_fastq = NULL,
                             barcode_csv = NULL,
                             organism,
                             reference,
                             feature_type,
                             remove_duplicates = FALSE,
                             samtools_path = NULL,
                             genome_size   = NULL,
                             bin_size      = NULL,
                             yieldsize     = 10000000,
                             mapq          = 30,
                             exclude_regions = TRUE,
                             excluded_regions_filename = NULL,
                             fix_chr = "none",
                             lower = NULL,
                             cell_calling = "filter",
                             promoters_file = NULL,
                             tss_file       = NULL,
                             enhs_file      = NULL,
                             gene_anno_file = NULL,
                             min_uniq_frags = 3000,
                             max_uniq_frags = 50000,
                             min_frac_peak = 0.3,
                             min_frac_tss = 0,
                             min_frac_enhancer = 0,
                             min_frac_promoter = 0.1,
                             max_frac_mito = 0.15,
                             report = TRUE,
                             nthreads = 12,
                             output_folder = NULL) {
  
  get_filename_without_extension <- function(path, extension_length = 1) {
    vec <- strsplit(basename(path), "\\.")[[1]]
    name.size <- length(vec) - extension_length
    return(paste(vec[1:name.size], collapse = "."))
  }

  r1_name <- get_filename_without_extension(r1, extension_length = 2)
  r2_name <- get_filename_without_extension(r2, extension_length = 2)
  
  if (!is.null(barcode_fastq)) {
    sc_atac_trim_barcode (r1            = r1,
                          r2            = r2,
                          bc_file       = barcode_fastq,
                          rmN           = TRUE,
                          rmlow         = TRUE,
                          output_folder = output_folder)
  } else if (!is.null(barcode_csv)) {
    sc_atac_trim_barcode (r1            = r1,
                          r2            = r2,
                          bc_file       = barcode_csv,
                          id1_st = 0,
                          id1_len = 16,
                          id2_st = 0,
                          id2_len = 16,
                          rmN           = TRUE,
                          rmlow         = TRUE,
                          output_folder = output_folder)
  } else {
    return()
  }
  
  if (!is.null(barcode_fastq)) {
    demux_r1        <- file.path(output_folder, paste0("demux_", r1_name, ".fastq.gz"))
    demux_r2        <- file.path(output_folder, paste0("demux_", r2_name, ".fastq.gz"))
  } else {
    demux_r1        <- file.path(output_folder, paste0("demultiplexed_completematch_", r1_name, ".fastq.gz"))
    demux_r2        <- file.path(output_folder, paste0("demultiplexed_completematch_", r2_name, ".fastq.gz"))
  }
  reference       <- reference
  cat(demux_r1)
  sc_aligning(ref = reference,
              R1 = demux_r1,
              R2 = demux_r2,
              nthreads  = nthreads,
              output_folder = output_folder)

  if (!is.null(barcode_fastq)) {
    bam_to_tag  <- file.path(output_folder, paste0("demux_", r1_name, "_aligned.bam"))
  } else {
    bam_to_tag  <- file.path(output_folder, paste0("demultiplexed_completematch_", r1_name, "_aligned.bam"))
  }

  sc_atac_bam_tagging(inbam         = bam_to_tag,
                      output_folder = output_folder,
                      bam_tags      = list(bc="CB", mb="OX"),
                      nthreads      =  nthreads)
  
  if (!is.null(barcode_fastq)) {
    sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted.bam"))
  } else {
    sorted_tagged_bam <- file.path(output_folder, paste0("demultiplexed_completematch_", r1_name, "_aligned_tagged_sorted.bam"))
  }

  if (isTRUE(remove_duplicates)) {
    removed <- sc_atac_remove_duplicates(inbam = sorted_tagged_bam,
                              samtools_path = samtools_path,
                              output_folder = output_folder)
    removed <- TRUE
    if (!removed) return()
    if (!is.null(barcode_fastq)) {
      sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted.bam"))
    } else {
      sorted_tagged_bam <- file.path(output_folder, paste0("demultiplexed_completematch_", r1_name, "_aligned_tagged_sorted_markdup.bam"))
    }
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
                            genome_size   = genome_size,
                            promoters_file = promoters_file,
                            tss_file       = tss_file,
                            enhs_file      = enhs_file,
                            gene_anno_file = gene_anno_file,
                            bin_size      = bin_size,
                            yieldsize     = yieldsize,
                            mapq          = mapq,
                            exclude_regions = exclude_regions,
                            excluded_regions_filename = excluded_regions_filename,
                            output_folder = output_folder,
                            fix_chr       = fix_chr,
                            lower         = lower,
                            min_uniq_frags = min_uniq_frags,
                            max_uniq_frags = max_uniq_frags,
                            min_frac_peak = min_frac_peak,
                            min_frac_tss = min_frac_tss,
                            min_frac_enhancer = min_frac_enhancer,
                            min_frac_promoter = min_frac_promoter,
                            max_frac_mito = max_frac_mito)
  
  sce <- sc_atac_create_sce(input_folder = output_folder,
                            organism     = organism,
                            feature_type = feature_type,
                            pheno_data   = NULL,
                            report       = report)
  return(sce)
}

#' @name sc_atac_pipeline_quick_test
#' @title A function that tests the pipeline on a small test sample (without duplicate removal)
#'
sc_atac_pipeline_quick_test <- function() {
  data.folder <- system.file("extdata", package = "scPipe", mustWork = TRUE)
  output_folder <- file.path(getwd(), "scPipe-atac-output")
  out <- tryCatch(
    {
      sce <- sc_atac_pipeline(r1 = file.path(data.folder, "small_chr21_R1.fastq.gz"),
                              r2 = file.path(data.folder, "small_chr21_R3.fastq.gz"),
                              barcode_fastq = file.path(data.folder, "small_chr21_R2.fastq.gz"),
                              organism = "hg38",
                              reference = file.path(data.folder, "small_chr21.fa"),
                              feature_type = "peak",
                              remove_duplicates = FALSE,
                              min_uniq_frags = 0,
                              min_frac_peak = 0,
                              min_frac_promoter = 0,
                              output_folder = output_folder)
      cat("Successfully ran pipeline.\n")
    },
    error = function(e) {
      message(e)
    },
    finally = {
      # system2("rm", c("-rf", output_folder))
    }
  )
}
