# --------- Trim the barcodes ---------
input_folder_test <- "/stornext/Home/data/allstaff/y/yang.p/scPipe_testing/test_data"
output_folder <- "/stornext/Home/data/allstaff/y/yang.p/scPipe_testing/scPipe_atac_output_csv"
get_filename_without_extension <- function(path, extension_length = 1) {
  vec <- strsplit(basename(path), "\\.")[[1]]
  name.size <- length(vec) - extension_length
  return(paste(vec[1:name.size], collapse = "."))
}
r1 = file.path(input_folder_test, "testfastq_S1_L001_R1_001_withbarcode.fastq.gz")
r2 = file.path(input_folder_test, "testfastq_S1_L001_R3_001_withbarcode.fastq.gz")
barcode_csv = file.path(input_folder_test, "testfastq_modified_barcode_1000.csv")

r1_nobc = file.path(input_folder_test, "testfastq_S1_L001_R1_001.fastq.gz")
r3_nobc = file.path(input_folder_test, "testfastq_S1_L001_R3_001.fastq.gz")

barcode_fastq = file.path(input_folder_test, "testfastq_S1_L001_R2_001.fastq.gz")

r1_name <- get_filename_without_extension(r1, extension_length = 2)
r2_name <- get_filename_without_extension(r2, extension_length = 2)

atac_bc_white_path <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/software/cellranger/cellranger-arc-2.0.0/lib/python/atac/barcodes/737K-arc-v1.txt.gz"

sc_atac_trim_barcode (r1            = r1_nobc,
                      r2            = r3_nobc,
                      bc_file       = barcode_fastq,
                      id1_st = 0,
                      id1_len = 16,
                      id2_st = 0,
                      id2_len = 16,
                      rmN           = TRUE,
                      rmlow         = TRUE,
                      output_folder = output_folder)

sc_atac_trim_barcode (r1            = r1_nobc,
                      r2            = r3_nobc,
                      bc_file       = barcode_fastq,
                      valid_barcode_file = atac_bc_white_path,
                      rmN           = TRUE,
                      rmlow         = TRUE,
                      output_folder = output_folder)

demux_r1 <- file.path(output_folder, "demultiplexed_completematch_testfastq_S1_L001_R1_001_withbarcode.fastq.gz")
demux_r3 <- file.path(output_folder, "demultiplexed_completematch_testfastq_S1_L001_R3_001_withbarcode.fastq.gz")


sc_aligning(R1 = demux_r1,
            R2 = demux_r1,
            tech = "atac",
            ref = "/stornext/General/data/user_managed/grpu_mritchie_1/PhilYang/genome_data/hg38.fa",
            nthreads = 8,
            output_folder = output_folder)


data.folder <- "/stornext/Home/data/allstaff/y/yang.p/scPipe_testing/test_data"
output_folder <- "/stornext/Home/data/allstaff/y/yang.p/scPipe_testing/scPipe_atac_output_csv"
sce <- sc_atac_pipeline(r1 = file.path(data.folder, "testfastq_S1_L001_R1_001_withbarcode.fastq.gz"),
                        r2 = file.path(data.folder, "testfastq_S1_L001_R3_001_withbarcode.fastq.gz"),
                        barcode_csv = file.path(data.folder, "testfastq_modified_barcode_1000.csv"),
                        organism = "hg38",
                        reference = file.path(data.folder, "chr21.fa"),
                        feature_type = "peak",
                        remove_duplicates = TRUE,
                        min_uniq_frags = 0,
                        min_frac_enhancer = 0,
                        min_frac_peak = 0,
                        min_frac_promoter = 0,
                        samtools_path = "/stornext/System/data/apps/samtools/samtools-1.14/bin/samtools",
                        output_folder = output_folder)






# ------------------- plot demux stats
data <- read.csv(file.path(output_folder, "scPipe_atac_stats", "demultiplexing_stats.csv"))

DT::datatable(data, options=list(paging=FALSE, searching=FALSE))

data$prop <- data$count/sum(data$count) 

ggplot(data, aes_string(x="status", y="prop", fill="status")) + 
  scale_fill_brewer(palette="Set1") +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(y = prop, label = percent(prop)), vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x=element_text(angle = 50, hjust = 1),
        axis.ticks.x=element_blank(),
        panel.border = element_blank()) +  
  xlab("status") +
  ylab("percentage") +
  expand_limits(y = 1) +
  ggtitle(paste0("Overall alignment mapping statistics of demultiplexed reads"))





# -----Test valid

output_folder <- "/stornext/Home/data/allstaff/y/yang.p/scPipe_testing/scPipe_atac_output_csv"
demux_r1 <- file.path(output_folder, "demultiplexed_completematch_testfastq_S1_L001_R1_001_withbarcode.fastq.gz")
nomatch_r1 <- file.path(output_folder, "demultiplexed_nomatch_testfastq_S1_L001_R1_001_withbarcode.fastq.gz")
partial_match_r1 <- file.path(output_folder, "demultiplexed_partialmatch_testfastq_S1_L001_R1_001_withbarcode.fastq.gz")
length(readLines(nomatch_r1))
