## ---- message=FALSE-----------------------------------------------------------
library(scPipe)
library(SingleCellExperiment)
data_dir = tempdir()

## -----------------------------------------------------------------------------
# file path:
ERCCfa_fn = system.file("extdata", "ERCC92.fa", package = "scPipe")
ERCCanno_fn = system.file("extdata", "ERCC92_anno.gff3", package = "scPipe")
barcode_annotation_fn = system.file("extdata", "barcode_anno.csv", package = "scPipe")

## ----eval=TRUE----------------------------------------------------------------
fq_R1 = system.file("extdata", "simu_R1.fastq.gz", package = "scPipe")
fq_R2 = system.file("extdata", "simu_R2.fastq.gz", package = "scPipe")

## ----eval=TRUE----------------------------------------------------------------
sc_trim_barcode(file.path(data_dir, "combined.fastq.gz"),
                fq_R1,
                fq_R2,
                read_structure = list(bs1=-1, bl1=0, bs2=6, bl2=8, us=0, ul=6))

## ----eval=TRUE----------------------------------------------------------------
if(.Platform$OS.type != "windows"){
  Rsubread::buildindex(basename=file.path(data_dir, "ERCC_index"), reference=ERCCfa_fn)

  Rsubread::align(index=file.path(data_dir, "ERCC_index"),
      readfile1=file.path(data_dir, "combined.fastq.gz"),
      output_file=file.path(data_dir, "out.aln.bam"), phredOffset=64)
}

## ----eval=TRUE----------------------------------------------------------------
if(.Platform$OS.type != "windows"){
  sc_exon_mapping(file.path(data_dir, "out.aln.bam"),
                file.path(data_dir, "out.map.bam"),
                ERCCanno_fn)
}

## ----eval=TRUE----------------------------------------------------------------
if(.Platform$OS.type != "windows"){
  sc_demultiplex(file.path(data_dir, "out.map.bam"), data_dir, barcode_annotation_fn, has_UMI=FALSE)

  sc_gene_counting(data_dir, barcode_annotation_fn)
}

## -----------------------------------------------------------------------------
if(.Platform$OS.type != "windows"){
sce = create_sce_by_dir(data_dir)
dim(sce)
}

## -----------------------------------------------------------------------------
data("sc_sample_data")
data("sc_sample_qc")
sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data))) # generate new sce with gene count matrix
QC_metrics(sce) = sc_sample_qc
demultiplex_info(sce) = cell_barcode_matching
UMI_dup_info(sce) = UMI_duplication

## ---- fig.height=7, fig.width=7-----------------------------------------------
plot_demultiplex(sce)

## ---- fig.height=7, fig.width=7-----------------------------------------------
plot_UMI_dup(sce)

## ---- warning=FALSE, message=FALSE--------------------------------------------
sce = calculate_QC_metrics(sce)
sce = detect_outlier(sce)

## ---- fig.height=7, fig.width=7-----------------------------------------------
plot_mapping(sce, percentage = TRUE, dataname = "sc_sample_data")

## ---- warning=FALSE, message=FALSE, fig.height=7, fig.width=7-----------------
plot_QC_pairs(sce)

## -----------------------------------------------------------------------------
sce = remove_outliers(sce)
dim(sce)

## -----------------------------------------------------------------------------
sessionInfo()

