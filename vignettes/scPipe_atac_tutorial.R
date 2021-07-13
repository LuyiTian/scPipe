## ---- message=FALSE-----------------------------------------------------------
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
#  sc_atac_trim_barcode (r1            = r1,
#                        r2            = r2,
#                        bc_file       = barcode_1000,
#                        id1_st        = -1,
#                        id1_len       = -1,
#                        id2_st        = -1,
#                        id2_len       = -10,
#                        output_folder = output_folder,
#                        rmN           = TRUE)

