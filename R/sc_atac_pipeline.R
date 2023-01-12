#' @name sc_atac_pipeline
#' 
#' @title A convenient function for running the entire pipeline
#' 
#' @param r1 The first read fastq file
#' @param r2 The second read fastq file
#' @param bc_file the barcode information, can be either in a \code{fastq} format (e.g. from 10x-ATAC) or
#' from a \code{.csv} file (here the barcode is expected to be on the second column). 
#' Currently, for the fastq approach, this can be a list of barcode files.
#' @param valid_barcode_file optional file path of the valid (expected) barcode sequences to be found in the bc_file (.txt, can be txt.gz). 
#' Must contain one barcode per line on the second column separated by a comma (default ="").
#' If given, each barcode from bc_file is matched against the barcode of
#' best fit (allowing a hamming distance of 1). If a FASTQ \code{bc_file} is provided, barcodes with a higher mapping quality, as given by
#' the fastq reads quality score are prioritised.
#' @param id1_st barcode start position (0-indexed) for read 1, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param id2_st barcode start position (0-indexed) for read 2, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param id1_len barcode length for read 1, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param id2_len barcode length for read 2, which is an extra parameter that is needed if the
#' \code{bc_file} is in a \code{.csv} format.
#' @param rmN ogical, whether to remove reads that contains N in UMI or cell barcode.
#' @param rmlow logical, whether to remove reads that have low quality barcode sequences.
#' @param organism The name of the organism e.g. hg38
#' @param reference The reference genome file
#' @param feature_type The feature type (either `genome_bin` or `peak`)
#' @param remove_duplicates Whether or not to remove duplicates (samtools is required) 
#' @param samtools_path A custom path of samtools to use for duplicate removal
#' @param bin_size The size of the bins for feature counting with the `genome_bin` feature type
#' @param yieldsize The number of reads to read in for feature counting
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
#' @returns None (invisible `NULL`)
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
#' \dontrun{
#' sc_atac_pipeline(
#'   r1 = r1,
#'   r2 = r2,
#'   bc_file = barcode_fastq
#' )
#' }
#'
#' @export
#' 
sc_atac_pipeline <- function(r1,
                            r2,
                            bc_file,
                            valid_barcode_file = "",
                            id1_st = -0,
                            id1_len = 16,
                            id2_st = 0,
                            id2_len = 16,
                            rmN           = TRUE,
                            rmlow         = TRUE,
                            organism = NULL,
                            reference = NULL,
                            feature_type = NULL,
                            remove_duplicates = FALSE,
                            samtools_path = NULL,
                            genome_size   = NULL,
                            bin_size      = NULL,
                            yieldsize     = 1000000,
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
        return(paste(vec[seq_len(name.size)], collapse = "."))
    }

    r1_name <- get_filename_without_extension(r1, extension_length = 2)
    r2_name <- get_filename_without_extension(r2, extension_length = 2)

    sc_atac_trim_barcode (r1            = r1,
                        r2            = r2,
                        bc_file = bc_file,
                        valid_barcode_file = valid_barcode_file,
                        id1_st = id1_st,
                        id1_len = id1_len,
                        id2_st = id2_st,
                        id2_len = id2_len,
                        rmN           = TRUE,
                        rmlow         = TRUE,
                        output_folder = output_folder)
    
    
    demux_r1        <- file.path(output_folder, paste0("demux_completematch_", r1_name, ".fastq.gz"))
    demux_r2        <- file.path(output_folder, paste0("demux_completematch_", r2_name, ".fastq.gz"))
    
    reference       <- reference
    bam_to_tag <- sc_aligning(ref = reference,
                    tech = "atac",
                    R1 = demux_r1,
                    R2 = demux_r2,
                    nthreads  = nthreads,
                    output_folder = output_folder)

    sorted_tagged_bam <- sc_atac_bam_tagging(inbam = bam_to_tag,
                        output_folder = output_folder,
                        bam_tags      = list(bc="CB", mb="OX"),
                        nthreads      =  nthreads)
    
    if (isTRUE(remove_duplicates)) {
        removed <- sc_atac_remove_duplicates(inbam = sorted_tagged_bam,
                                samtools_path = samtools_path,
                                output_folder = output_folder)
        if (!isFALSE(removed))
        sorted_tagged_bam <- removed
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

    sc_atac_feature_counting (fragment_file = file.path(output_folder, "fragments.bed"),
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
    
    # sce <- sc_atac_create_sce(input_folder = output_folder,
    #                           organism     = organism,
    #                           feature_type = feature_type,
    #                           pheno_data   = NULL,
    #                           report       = report)
    # return(sce)
}

#' @name sc_atac_pipeline_quick_test
#' @title A function that tests the pipeline on a small test sample (without duplicate removal)
#' @returns None (invisible `NULL`)
sc_atac_pipeline_quick_test <- function() {
    data.folder <- system.file("extdata", package = "scPipe", mustWork = TRUE)
    output_folder <- file.path(getwd(), "scPipe-atac-output")
    out <- tryCatch({
        sce <- sc_atac_pipeline(r1 = file.path(data.folder, "small_chr21_R1.fastq.gz"),
                                r2 = file.path(data.folder, "small_chr21_R3.fastq.gz"),
                                bc_file = file.path(data.folder, "small_chr21_R2.fastq.gz"),
                                organism = "hg38",
                                reference = file.path(data.folder, "small_chr21.fa"),
                                feature_type = "peak",
                                remove_duplicates = FALSE,
                                min_uniq_frags = 0,
                                min_frac_peak = 0,
                                min_frac_promoter = 0,
                                output_folder = output_folder)
        message("Successfully ran pipeline.")
        },
        error = function(e) {
        message(e)
        },
        finally = {
        # system2("rm", c("-rf", output_folder))
        }
    )
}
