## ---- echo=FALSE, warning=FALSE, message=FALSE, out.width="100%"--------------
knitr::include_graphics(system.file("scATAC-seq_workflow.png", package = "scPipe"))

## ---- message=FALSE-----------------------------------------------------------
library(scPipe)
data.folder <- system.file("extdata", package = "scPipe", mustWork = TRUE)


## -----------------------------------------------------------------------------
#output_folder <- tempdir()
output_folder <- "D:/WEHI/output/vignette"

## -----------------------------------------------------------------------------
# data fastq files
r1      <- file.path(data.folder, "small_chr21_R1.fastq.gz") 
r2      <- file.path(data.folder, "small_chr21_R3.fastq.gz") 

# barcodes in fastq format
barcode_fastq      <- file.path(data.folder, "small_chr21_R2.fastq.gz") 

# barcodes in .csv format
barcode_1000       <- file.path(data.folder, "chr21_modified_barcode_1000.csv")


