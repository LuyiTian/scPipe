## ---- message=FALSE-----------------------------------------------------------
devtools::load_all()
library(scPipe)
data.folder <- system.file("data", package = "scPipe", mustWork = TRUE)


## -----------------------------------------------------------------------------
output_folder <- file.path(data.folder, "scPipe-atac-output")



## ----eval=TRUE----------------------------------------------------------------
# data fastq files
r1      <- file.path(data.folder, "testfastq_S1_L001_R1_001.fastq.gz") 
r2      <- file.path(data.folder, "testfastq_S1_L001_R3_001.fastq.gz") 

get_filename_without_extension <- function(path, extension_length = 1) {
  vec <- strsplit(basename(path), "\\.")[[1]]
  name.size <- length(vec) - extension_length
  return(paste(vec[1:name.size], collapse = "."))
}

r1_name <- get_filename_without_extension(r1, extension_length = 2)
r2_name <- get_filename_without_extension(r2, extension_length = 2)

reference       <- file.path(data.folder, "genome.fa")

# barcodes in fastq format
barcode_fastq      <- file.path(data.folder, "testfastq_S1_L001_R2_001.fastq.gz") 

# barcodes in .csv format
barcode_1000       <- file.path(data.folder, "testfastq_modified_barcode_1000.csv")



## ----eval=TRUE----------------------------------------------------------------
sc_atac_trim_barcode (r1            = r1, 
                      r2            = r2, 
                      bc_file       = barcode_fastq, 
                      rmN           = TRUE,
                      rmlow         = TRUE,
                      output_folder = output_folder)


## ----eval=FALSE---------------------------------------------------------------
## sc_atac_trim_barcode (r1            = r1,
##                       r2            = r2,
##                       bc_file       = barcode_1000,
##                       id1_st        = -1,
##                       id1_len       = -1,
##                       id2_st        = -1,
##                       id2_len       = -10,
##                       output_folder = output_folder,
##                       rmN           = TRUE)


## ----eval=TRUE----------------------------------------------------------------

demux_r1        <- file.path(output_folder, paste0("demux_", r1_name, ".fastq.gz"))
demux_r2        <- file.path(output_folder, paste0("demux_", r2_name, ".fastq.gz"))



sc_atac_aligning(ref       = reference, 
                 readFile1 = demux_r1, 
                 readFile2 = demux_r2, 
                 nthreads  = 12,
                 output_folder = output_folder)


## ---- eval=TRUE---------------------------------------------------------------
bam_to_tag  <- file.path(output_folder, paste0("demux_", r1_name, "_aligned.bam"))

sc_atac_bam_tagging (inbam         = bam_to_tag, 
                     output_folder =  output_folder, 
                     bam_tags      = list(bc="CB", mb="OX"), 
                     nthreads      =  12)


## ---- eval=TRUE---------------------------------------------------------------
sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted.bam"))

sc_atac_remove_duplicates(sorted_tagged_bam,
                          output_folder = output_folder)

sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted_markdup.bam"))

if (!file.exists(sorted_tagged_bam))
  sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted.bam"))



## -----------------------------------------------------------------------------
sc_atac_create_fragments(inbam = sorted_tagged_bam,
                         output_folder = output_folder)



## -----------------------------------------------------------------------------
sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted.bam"))
sc_atac_peak_calling(inbam = sorted_tagged_bam, 
                     ref = reference,
                     genome_size = NULL,
                     output_folder = output_folder)



## ----eval=TRUE----------------------------------------------------------------
sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted.bam"))
features          <- file.path(output_folder, "NA_peaks.narrowPeak")

sc_atac_feature_counting (insortedbam   = sorted_tagged_bam,
                          feature_input = features, 
                          bam_tags      = list(bc="CB", mb="OX"), 
                          feature_type  = "peak", 
                          organism      = "hg38",
                          cell_calling  = FALSE,
                          genome_size   = NULL,
                          bin_size      = NULL, 
                          yieldsize     = 1000000,
                          mapq          = 20,
                          exclude_regions = TRUE,
                          output_folder = output_folder,
                          fix_chr       = "none"
                          )


## -----------------------------------------------------------------------------
reference <- file.path(data.folder, "genome.fa")

sc_atac_feature_counting (insortedbam   = sorted_tagged_bam,
                          feature_input = reference,
                          bam_tags      = list(bc="CB", mb="OX"),
                          feature_type  = "genome_bin",
                          organism      = "hg38",
                          cell_calling  = FALSE,
                          genome_size   = NULL,
                          bin_size      = NULL,
                          yieldsize     = 1000000,
                          mapq          = 20,
                          exclude_regions = TRUE,
                          output_folder = output_folder,
                          fix_chr       = "none"
                          )



## -----------------------------------------------------------------------------
sorted_tagged_bam <- file.path(output_folder, paste0("demux_", r1_name, "_aligned_tagged_sorted.bam"))
features          <- file.path(output_folder, "NA_peaks.narrowPeak")

sc_atac_feature_counting (insortedbam   = sorted_tagged_bam,
                          feature_input = features, 
                          bam_tags      = list(bc="CB", mb="OX"), 
                          feature_type  = "peak",
                          organism      = "hg38",
                          cell_calling  = "cellranger",
                          genome_size   = NULL,
                          bin_size      = NULL, 
                          yieldsize     = 1000000,
                          mapq          = 30,
                          exclude_regions = TRUE,
                          output_folder = output_folder,
                          fix_chr       = "none"
                          )


## ----eval=FALSE---------------------------------------------------------------
## feature_matrix <- readRDS(file.path(output_folder, "feature_matrix.rds"))
## dplyr::giimpse(feature_matrix)
## 
## sparseM <- readMM(file.path(output_folder, "sparse_matrix.mtx"))
## dplyr::glimpse(sparseM)


## ---- eval=FALSE--------------------------------------------------------------
## sce <- sc_atac_create_sce(input_folder = output_folder,
##                    organism     = "hg38",
##                    feature_type = "peak",
##                    pheno_data   = NULL,
##                    report       = TRUE)


## -----------------------------------------------------------------------------
sessionInfo()

