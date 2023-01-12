## ---- echo=FALSE, warning=FALSE, message=FALSE, out.width="100%"--------------
library(knitr)

knitr::include_graphics(system.file("scATAC-seq_workflow.png", package = "scPipe"))

## ---- message=FALSE-----------------------------------------------------------
library(scPipe)
data.folder <- system.file("extdata", package = "scPipe", mustWork = TRUE)


## -----------------------------------------------------------------------------
output_folder <- tempdir()

## -----------------------------------------------------------------------------
# data fastq files
r1      <- file.path(data.folder, "small_chr21_R1.fastq.gz") 
r2      <- file.path(data.folder, "small_chr21_R3.fastq.gz") 

# barcodes in fastq format
barcode_fastq      <- file.path(data.folder, "small_chr21_R2.fastq.gz") 

# barcodes in .csv format
barcode_1000       <- file.path(data.folder, "chr21_modified_barcode_1000.csv")


## -----------------------------------------------------------------------------
sc_atac_trim_barcode (r1            = r1, 
                      r2            = r2, 
                      bc_file       = barcode_fastq, 
                      rmN           = TRUE,
                      rmlow         = TRUE,
                      output_folder = output_folder)

## ----eval=FALSE---------------------------------------------------------------
#  sc_atac_trim_barcode (r1            = r1,
#                        r2            = r2,
#                        bc_file       = barcode_1000,
#                        id1_st        = -1,
#                        id1_len       = -1,
#                        id2_st        = -1,
#                        id2_len       = -10,
#                        output_folder = output_folder,
#                        rmN           = TRUE)

## -----------------------------------------------------------------------------
demux_r1        <- file.path(output_folder, "demux_completematch_small_chr21_R1.fastq.gz")
demux_r2        <- file.path(output_folder, "demux_completematch_small_chr21_R3.fastq.gz")
reference = file.path(data.folder, "small_chr21.fa")


aligned_bam <- sc_aligning(ref = reference, 
                R1 = demux_r1, 
                R2 = demux_r2, 
                nthreads  = 12,
                output_folder = output_folder)

## -----------------------------------------------------------------------------

sorted_tagged_bam <- sc_atac_bam_tagging (inbam = aligned_bam, 
                       output_folder =  output_folder, 
                       bam_tags      = list(bc="CB", mb="OX"), 
                       nthreads      =  12)


## ---- eval=FALSE--------------------------------------------------------------
#  
#  removed <- sc_atac_remove_duplicates(sorted_tagged_bam,
#                            output_folder = output_folder)
#  if (!isFALSE(removed))
#    sorted_tagged_bam <- removed
#  

## -----------------------------------------------------------------------------
sc_atac_create_fragments(inbam = sorted_tagged_bam,
                         output_folder = output_folder)


## -----------------------------------------------------------------------------
sc_atac_peak_calling(inbam = sorted_tagged_bam, 
                     ref = reference,
                     genome_size = NULL,
                     output_folder = output_folder)


## ----eval=FALSE---------------------------------------------------------------
#  features          <- file.path(output_folder, "NA_peaks.narrowPeak")
#  
#  sc_atac_feature_counting (fragment_file = file.path(output_folder, "fragments.bed"),
#                            feature_input = features,
#                            bam_tags      = list(bc="CB", mb="OX"),
#                            feature_type  = "peak",
#                            organism      = "hg38",
#                            yieldsize     = 1000000,
#                            mapq          = 20,
#                            exclude_regions = TRUE,
#                            output_folder = output_folder,
#                            fix_chr       = "none"
#                            )

## ----eval=FALSE---------------------------------------------------------------
#  reference       <- file.path(data.folder, "small_chr21.fa")
#  sc_atac_feature_counting (fragment_file = file.path(output_folder, "fragments.bed"),
#                            feature_input = reference,
#                            bam_tags      = list(bc="CB", mb="OX"),
#                            feature_type  = "genome_bin",
#                            organism      = "hg38",
#                            cell_calling  = FALSE,
#                            yieldsize     = 1000000,
#                            exclude_regions = TRUE,
#                            output_folder = output_folder,
#                            fix_chr       = "none"
#                            )
#  

## ----eval = FALSE-------------------------------------------------------------
#  features          <- file.path(output_folder, "NA_peaks.narrowPeak")
#  
#  sc_atac_feature_counting (fragment_file = file.path(output_folder, "fragments.bed"),
#                            feature_input = features,
#                            bam_tags      = list(bc="CB", mb="OX"),
#                            feature_type  = "peak",
#                            organism      = "hg38",
#                            cell_calling  = "filter",
#                            min_uniq_frags = 0,
#                            min_frac_peak = 0,
#                            min_frac_promoter = 0,
#                            yieldsize     = 1000000,
#                            exclude_regions = TRUE,
#                            output_folder = output_folder,
#                            fix_chr       = "none"
#                            )

