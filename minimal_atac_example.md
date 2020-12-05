# Minimal example -----------------

Open RStudio

Load the project
```r
#devtools::install()
#devtools::load_all()
library(here)
library(scPipe)
devtools::load_all(here())
```
Install the packages from sc_atac_setup.R

-- However this package installation needs to be automatically handled within the package later

# Demultiplexing -----------------

## define inputs

```r
# original read files
r1      <- here(
  "data","testfastq_S1_L001_R1_001.fastq.gz") 
r2      <- here(
  "data","testfastq_S1_L001_R3_001.fastq.gz") 
  
# read files with barcode starting from 1st position
r1_barcode <- here(
  "data","testfastq_S1_L001_R1_001_withbarcode.fastq.gz") 
r2_barcode <- here(
  "data","testfastq_S1_L001_R3_001_withbarcode.fastq.gz") 

# barcodes in fastq format
barcode_fastq      <- here(
  "data","testfastq_S1_L001_R2_001.fastq.gz") 

# barcodes in .csv format
barcode_csv_sample <- here("data","barcode.csv")
barcode_csv        <- here("data", "testfastq_modified_barcode.csv") 
barcode_csv_small  <- here("data", "testfastq_modified_barcode_small.csv") 
barcode_csv_incorr <- here("data", "testfastq_modified_barcode_corrupt.csv") 

```
### if using a barcode fastq file ----------------

```r
sc_atac_trim_barcode (r1 = r1, r2 =  r2, bc_file =  barcode_fastq, output_folder = "")
```
### if using a barcode .csv file -----------------

program expects the barcode to be on the second column of the .csv file

```r
sc_atac_trim_barcode(r1 = r1_barcode, r2 = r2_barcode, bc_file = barcode_csv_small, bc_start = 3, bc_length = 16, output_folder = "", rmN = FALSE)
```
# Aligning to reference -----------------

## define inputs

```r
reference       <- here("data", "genome.fa")
demux_r1        <- here("scPipe-atac-output", "demux_testfastq_S1_L001_R1_001.fastq.gz")
demux_r2        <- here("scPipe-atac-output", "demux_testfastq_S1_L001_R3_001.fastq.gz")
```
## run function

```r
sc_atac_aligning(ref = reference, readFile1 = demux_r1, readFile2 = demux_r2)
```
# Tagging the aligned BAM file -----------------

```r
bam_to_tag  <- here("scPipe-atac-output", "demux_testfastq_S1_L001_R1_001_aligned_sorted.bam")

sc_atac_bam_tagging(inbam = bam_to_tag, outfolder="", bam_tags = list(bc="CB", mb="OX"), nthreads=6)
```
# Feature counting (currently being actively developed, so beware when trying)-----------------

```r
sorted_tagged_bam <- here("scPipe-atac-output", "demux_testfastq_S1_L001_R1_001_aligned_sorted_tagged_sorted.bam")
features          <- here("data", "extdata", "NA_peaks.narrowPeak")

sc_atac_feature_counting (insortedbam      = sorted_tagged_bam,
                          feature_input    = features, 
                          bam_tags         = list(bc="CB", mb="OX"), 
                          feature_type     = "peak", 
                          organism         = NULL,
                          cell_calling     = "emptydrops",
                          genome_size      = NULL,
                          bin_size         = NULL, 
                          qc_per_bc_file   = NULL,
                          yieldsize        = 1000000,
                          mapq             = 20,
                          exclude_regions  = NULL,
                          output_folder    = "",
                          fix_chr          = "none" # should be either one of these: c("none", "blacklist", "feature", "both")
                          )

```
# Currently ends here...
