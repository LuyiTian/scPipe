

#' sc_trim_barcode
#'
#' @description This function will reformat the fastq file and move the barcode 
#' and UMI sequence to read names.
#' 
#' @details The default read structure in this function represents CEL-seq 
#' paired-ended reads, with one cell barcode in read two start from 6bp and 
#' UMI sequence in read two start from the first bp. So the read structure 
#' will be : `list(bs1=-1, bl1=0, bs2=6, bl2=8, us=0, ul=6)`. `bs1=-1, bl1=0` 
#' means we dont have index in read one so we set a negative value to start 
#' position and zero to the length. `bs2=6, bl2=8` means we have index in 
#' read two which start at 6bp with 8bp in its length. `us=0, ul=6` means 
#' we have UMI from the start of read two and the length in 6bp. NOTE: we 
#' use the zero based index so the index of the sequence start from zero. 
#' For a typical Drop-seq experiment the setting will be 
#' `list(bs1=-1, bl1=0, bs2=0, bl2=12, us=12, ul=8)`, which means the 
#' read one only contains transcript, the first 12bp in read two are index, 
#' followed by 8bp UMIs.
#'
#' @name sc_trim_barcode
#' @param outfq the output fastq file, which reformat the barcode and UMI 
#' into the read name.
#' @param r1 read one for pair-end reads. This read should contains the transcript.
#' @param r2 read two for pair-end reads, default to be `NULL` for single reads.
#' @param read_structure a list contains read structure configuration:\itemize{
#'  \item{bs1} is the starting position of barcode in read one, if there 
#'  is no reads in read one set it to -1.
#'  \item{bl1} is the length of barcode in read one, if there is no 
#'  barcode in read one this number is used 
#'  for trimming the beginning of read one.
#'  \item{bs2} is the starting position of barcode in read two
#'  \item{bl2} is the length of barcode in read two
#'  \item{us} is the starting position of UMI
#'  \item{ul} is the length of UMI
#'  }
#' @param filter_settings A list contains read filter settings:\itemize{
#'  \item{rmlow} whether to remove the low quality reads.
#'  \item{rmN} whether to remove reads that contains N in UMI or cell barcode.
#'  \item{minq} the minimum base pair quality that we allowed
#'  \item{numbq} the maximum number of base pair that have quality 
#'  below \code{numbq}
#'  }
#' @export
#' @return no return
#'
#' @examples
#' data_dir="celseq2_demo"
#' \dontrun{
#' # for the complete workflow, refer to the vignettes
#' ...
#' sc_trim_barcode(file.path(data_dir, "combined.fastq"),
#'    file.path(data_dir, "simu_R1.fastq"),
#'    file.path(data_dir, "simu_R2.fastq"))
#' ...
#' }
sc_trim_barcode = function(outfq, r1, r2=NULL,
                           read_structure = list(
                             bs1=-1, bl1=0, bs2=6, bl2=8, us=0, ul=6),
                           filter_settings = list(
                             rmlow=TRUE, rmN=TRUE, minq=20, numbq=2)) {
  if (filter_settings$rmlow) {
    i_rmlow = 1
  }
  else {
    i_rmlow = 0
  }
  if (filter_settings$rmN) {
    i_rmN = 1
  }
  else {
    i_rmN = 0
  }

  if (!is.null(r2)) {
    if (!file.exists(r1)) {stop("read1 fastq file does not exists.")}
    if (!file.exists(r2)) {stop("read2 fastq file does not exists.")}
    rcpp_sc_trim_barcode_paired(outfq, r1, r2,
                                read_structure$bs1,
                                read_structure$bl1,
                                read_structure$bs2,
                                read_structure$bl2,
                                read_structure$us,
                                read_structure$ul,
                                i_rmlow,
                                i_rmN,
                                filter_settings$minq,
                                filter_settings$numbq)
  }
  else {
    stop("not implemented.")
  }
}


#' sc_exon_mapping
#' 
#' @description This function will take the alinged read and map them to exons.
#' The result will be written into optional fields in bam file with different 
#' tags. Following this link for more information regarding to bam file format:
#' http://samtools.github.io/hts-specs
#' 
#' @name sc_exon_mapping
#' @param inbam input bam file. This is the output of alignment program.
#' @param outbam output bam file with gene and barcode tag
#' @param annofn genome annotation gff file. It can have multiple files 
#' names in a vector.
#' @param am mapping status tag (default: YE)
#' @param ge gene id tag (default: GE)
#' @param bc cell barcode tag (default: BC)
#' @param mb molecular barcode tag (default: OX)
#' @param bc_len total barcode length
#' @param UMI_len UMI length
#' @param stnd perform strand specific mapping or not
#' @param fix_chr add `chr` to chromosome names, fix inconsistant names.
#'
#' @export
#' @return no return
#' @examples
#' data_dir="celseq2_demo"
#' ERCCanno_fn = system.file("extdata", "ERCC92_anno.gff3", 
#'     package = "scPipe")
#' \dontrun{
#' # for the complete workflow, refer to the vignettes
#' ...
#' sc_exon_mapping(file.path(data_dir, "out.aln.bam"),
#'                 file.path(data_dir, "out.map.bam"),
#'                 ERCCanno_fn)
#' ...
#' }
#' 
sc_exon_mapping = function(inbam, outbam, annofn,
                            am="YE", ge="GE", bc="BC", mb="OX",
                            bc_len=8, UMI_len=6, stnd=TRUE, fix_chr=FALSE) {
  if (stnd) {
    i_stnd = 1
  }
  else {
    i_stnd = 0
  }
  if (fix_chr) {
    i_fix_chr = 1
  }
  else {
    i_fix_chr = 0
  }

  if (!file.exists(inbam)) {stop("input bam file does not exists.")}
  if (!file.exists(annofn)) {stop("genome annotation file does not exists.")}

  rcpp_sc_exon_mapping(inbam, outbam, annofn, am, ge, bc, mb, bc_len, 
                       UMI_len, stnd, fix_chr)
}


#' sc_demultiplex
#' 
#' @description This function will process bam file by cell barcode, 
#' output to outdir/count/[cell_id].csv.
#' the output contains information for all reads that can be mapped to exons.
#' including the gene id, UMI of that read and the distance to transcript end 
#' position.
#'
#' @param inbam input bam file. This should be the output of 
#' \code{sc_exon_mapping}
#' @param outdir output folder
#' @param bc_anno barcode annotation, first column is cell id, second column 
#' is cell barcode sequence
#' @param max_mis maximum mismatch allowed in barcode. Default to be 1
#' @param am mapping status tag (default: YE)
#' @param ge gene id tag (default: GE)
#' @param bc cell barcode tag (default: BC)
#' @param mb molecular barcode tag (default: OX)
#' @param mito mitochondrial chromosome name. 
#' This should be consistant with the chromosome names in the bam file.
#' @param has_UMI whether the protocol contains UMI (default: TRUE)
#'
#' @export
#' @return no return
#' @examples
#' data_dir="celseq2_demo"
#' barcode_annotation_fn = system.file("extdata", "barcode_anno.csv", 
#'     package = "scPipe")
#' \dontrun{
#' # refer to the vignettes for the complete workflow
#' ...
#' sc_demultiplex(file.path(data_dir, "out.map.bam"), 
#'     data_dir, 
#'     barcode_annotation_fn,has_UMI=FALSE)
#' ...
#' }
#' 
sc_demultiplex = function(inbam, outdir, bc_anno,
                          max_mis=1,
                          am="YE", ge="GE",
                          bc="BC",
                          mb="OX",
                          mito="MT",
                          has_UMI=TRUE) {
  dir.create(file.path(outdir, "count"), showWarnings = FALSE)
  dir.create(file.path(outdir, "stat"), showWarnings = FALSE)
  if (!file.exists(inbam)) {stop("input bam file does not exists.")}
  if (!file.exists(bc_anno)) {stop("barcode annotation file does not exists.")}
  rcpp_sc_demultiplex(inbam, outdir, bc_anno, max_mis, am, ge, bc, mb, 
                      mito, has_UMI)
}


#' sc_gene_counting
#' 
#' @description This function will merge UMI and generate gene count matrix
#'
#' @param outdir output folder that contains \code{sc_demultiplex} output
#' @param bc_anno barcode annotation, first column is cell id, second column 
#' is cell barcode sequence
#' @param UMI_cor correct UMI sequence error: 0 means no correction, 1 means 
#' simple correction and merge UMI with distance 1.
#' @param gene_fl whether to remove low abundant gene count. Low abundant is 
#' defined as only one copy of one UMI for this gene
#'
#' @export
#' @return no return
#' @examples
#' data_dir="celseq2_demo"
#' barcode_annotation_fn = system.file("extdata", "barcode_anno.csv", 
#' package = "scPipe")
#' \dontrun{
#' # refer to the vignettes for the complete workflow
#' ...
#' sc_gene_counting(data_dir, barcode_annotation_fn)
#' ...
#' }
#' 
sc_gene_counting = function(outdir, bc_anno, UMI_cor=1, gene_fl=FALSE) {
  dir.create(file.path(outdir, "count"), showWarnings = FALSE)
  dir.create(file.path(outdir, "stat"), showWarnings = FALSE)
  if (gene_fl) {
    i_gene_fl = 1
  }
  else {
    i_gene_fl = 0
  }
  if (!file.exists(bc_anno)) {stop("barcode annotation file does not exists.")}
  rcpp_sc_gene_counting(outdir, bc_anno, UMI_cor, i_gene_fl)
}



#' sc_detect_bc
#' 
#' @description Detect cell barcode and generate the barcode annotation
#'
#' @param infq input fastq file, shoule be the output file of 
#' \code{sc_trim_barcode}
#' @param outcsv output barcode annotation
#' @param suffix the suffix of cell name, default to be `CELL_`. The cell name
#' will be CELL_001, CELL_002 accordingly.
#' @param bc_len the length of cell barcode, should be consistent with bl1+bl2
#' in \code{sc_trim_barcode}
#' @param max_reads the maximum of reads processed, default is 1000,000, set to
#' "all" to process all reads (may spend more time)
#' @param min_count minimum counts to keep, barcode will be discarded if 
#' it has lower count. Default value is 10. This should be set according 
#' to \code{max_reads}.
#' @param max_mismatch the maximum mismatch allowed. Barcodes within this 
#' number will be considered as sequence error and merged, default to be 1.
#' @export
#' @return no return
#' @examples
#' \dontrun{
#' # `sc_detect_bc`` should run before `sc_demultiplex` for 
#' # Drop-seq or 10X protocols
#' sc_detect_bc("input.fastq","output.cell_index.csv",bc_len=8)
#' sc_demultiplex(...,"output.cell_index.csv")
#' }
#' 
#' 
sc_detect_bc = function(infq, outcsv, suffix="CELL_", bc_len, 
                        max_reads=1000000, min_count=10, max_mismatch=1) {
  if (!file.exists(infq)) {stop("input fastq file does not exists.")}
  if (max_reads=="all") {m_r = -1}
  else {m_r=max_reads}
  rcpp_sc_detect_bc(infq, outcsv, suffix, bc_len, m_r, min_count, max_mismatch)
}
