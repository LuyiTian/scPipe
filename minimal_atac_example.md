# Minimal example -----------------

Open RStudio

Load the project

```r
devtools::load_all()
```

# Demultiplexing 

## define inputs

```r
data.folder <- system.file("data", package = "scPipe", mustWork = TRUE)

# original read files
r1 <- file.path(data.folder,"testfastq_S1_L001_R1_001.fastq.gz") 
r2 <- file.path(data.folder,"testfastq_S1_L001_R3_001.fastq.gz") 
  
# read files with barcode starting from 1st position
r1_barcode <- file.path(data.folder,"testfastq_S1_L001_R1_001_withbarcode.fastq.gz") 
r2_barcode <- file.path(data.folder,"testfastq_S1_L001_R3_001_withbarcode.fastq.gz") 

# barcodes in fastq format
barcode_fastq      <- file.path(data.folder, "testfastq_S1_L001_R2_001.fastq.gz") 

# barcodes in .csv format
barcode_csv_sample <- file.path(data.folder,"barcode.csv")
barcode_csv        <- file.path(data.folder, "testfastq_modified_barcode.csv") 
barcode_csv_small  <- file.path(data.folder, "testfastq_modified_barcode_small.csv") 
barcode_csv_incorr <- file.path(data.folder, "testfastq_modified_barcode_corrupt.csv") 

output_folder <- NULL

```
### if using a barcode fastq file 

```r
sc_atac_trim_barcode (r1            = r1, 
                      r2            = r2, 
                      bc_file       = barcode_fastq, 
                      output_folder = output_folder)
```
### if using a barcode .csv file 

program expects the barcode to be on the second column of the .csv file

```r
sc_atac_trim_barcode(r1 = r1_barcode, 
                     r2 = r2_barcode, 
                     bc_file = barcode_csv_small, 
                     bc_start = 3, bc_length = 16, 
                     output_folder = output_folder, 
                     rmN = FALSE)
```
# Aligning to reference 

## define inputs

```r
reference       <- file.path(data.folder, "genome.fa")
demux_r1        <- file.path(output_folder, "demux_testfastq_S1_L001_R1_001.fastq.gz")
demux_r2        <- file.path(output_folder, "demux_testfastq_S1_L001_R3_001.fastq.gz")
```
## run function

```r
sc_atac_aligning(ref = reference, 
                 R1 = demux_r1, 
                 R2 = demux_r2, 
                 output_folder = output_folder,
                 nthreads = 6)
```
# Tagging the aligned BAM file

```r
bam_to_tag  <- file.path(output_folder, "demux_testfastq_S1_L001_R1_001_aligned.bam")

sc_atac_bam_tagging (inbam         = bam_to_tag, 
                     output_folder = output_folder, 
                     bam_tags      = list(bc="CB", mb="OX"), 
                     nthreads      = 6)
```

# Remove PCR duplicates (requires samtools)
```r
sorted_tagged_bam <- file.path(output_folder, "demux_testfastq_S1_L001_R1_001_aligned_tagged_sorted.bam")

sc_atac_remove_duplicates(inbam         = sorted_tagged_bam,
                          samtools_path = NULL, # can specify custom path here
                          output_folder = output_folder)

sorted_tagged_bam <- file.path(output_folder, "demux_testfastq_S1_L001_R1_001_aligned_tagged_sorted_markdup.bam")
```

# Gemerating a fragment file
```r
sc_atac_create_fragments(inbam = sorted_tagged_bam,
                         output_folder = output_folder)
```

# Peak calling
```r
sc_atac_peak_calling(inbam = sorted_tagged_bam, 
                     ref = reference,
                     genome_size = NULL,
                     output_folder = output_folder)
```

# Feature counting

```r
sorted_tagged_bam <- file.path(output_folder, "demux_testfastq_S1_L001_R1_001_aligned_tagged_sorted.bam")
features          <- file.path(data.folder, "extdata", "NA_peaks.narrowPeak")

sc_atac_feature_counting (insortedbam   = sorted_tagged_bam,
                          feature_input = features, 
                          bam_tags      = list(bc="CB", mb="OX"), 
                          feature_type  = "peak",
                          organism      = "hg38",
                          cell_calling  = "filter",
                          genome_size   = NULL,
                          bin_size      = NULL, 
                          yieldsize     = 1000000,
                          mapq          = 30,
                          exclude_regions = TRUE,
                          output_folder = output_folder,
                          fix_chr       = "none"
                          )
```

# Generating the *Single-cell Experiment (SCE)* object

```r
sce <- sc_atac_create_sce(input_folder = output_folder,
                          organism     = "hg38",
                          feature_type = "peak",
                          pheno_data   = NULL,
                          report       = TRUE
                          )
```
