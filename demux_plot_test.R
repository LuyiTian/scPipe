output_folder <- "/stornext/Home/data/allstaff/y/yang.p/scPipe_testing/scPipe_atac_output_csv"

cat("Found barcode non-matches demultiplexed FASTQ files. Proceeding to align them.\n")
R1 <- file.path(output_folder, "demux_partial_nomatch_R1.fastq.gz")
R3 <- file.path(output_folder, input_folder_files[grep("demultiplexed_nomatch.+R3.+fastq", input_folder_files)])
fileNameWithoutExtension <- paste0(output_folder, "/", strsplit(basename(R1), "\\.")[[1]][1])
outbam <- paste0(fileNameWithoutExtension, "_aligned.bam")
Rsubread::align(
  index = file.path(output_folder, "genome_index"),
  readfile1 = R1,
  readfile2 = R3,
  sortReadsByCoordinates = TRUE,
  type = "DNA",
  nthreads = 12,
  output_file = outbam)

# Extract columns


aligned_completematch_bam <- file.path(output_folder, "demultiplexed_completematch_testfastq_S1_L001_R1_001_withbarcode_aligned.bam")

bam_tags = list(bc="CB", mb="OX")
param <- Rsamtools::ScanBamParam(tag = as.character(bam_tags),  mapqFilter=20)
bamfl <- open(Rsamtools::BamFile(aligned_completematch_bam))
params <- Rsamtools::ScanBamParam(what=c("flag"), tag=c("CB"))
bam0 <- Rsamtools::scanBam(bamfl, param = params)

# Need to convert sam flag into alignment status

# Then, count number of reads with each alignment status

# Then add more row to the DF which represent the number of reads that do have a corresponding barcode