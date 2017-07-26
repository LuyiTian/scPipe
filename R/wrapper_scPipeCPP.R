

#' sc_trim_barcode
#'
#' reformat the fastq file and move the barcode and UMI to read names.
#'
#' @name sc_trim_barcode
#' @param outfq the output fastq file, which reformat the barcode and UMI into the read name
#' @param r1 read one for pair-end reads. This read should contains the transcript
#' @param r2 read two for pair-end reads. default to be `NULL` for single reads
#' @param read_structure a list contains read structure configuration:\itemize{
#'  \item{bs1} is the starting position of barcode in read one, if there is no reads in read one set it to -1.
#'  \item{bl1} is the length of barcode in read one, if there is no barcode in read one this number is used for trimming the beginning of read one.
#'  \item{bs2} is the starting position of barcode in read two
#'  \item{bl2} is the length of barcode in read two
#'  \item{us} is the starting position of UMI
#'  \item{ul} is the length of UMI
#'  }
#' @param filter_settings a list contains read filter settings:
#' @export
#'
#' @examples
#' #TODO
sc_trim_barcode <- function(outfq, r1, r2=NULL,
                           read_structure = list(bs1=-1,bl1=2, bs2=6, bl2=8, us=0, ul=6),
                           filter_settings = list(rmlow=TRUE, rmN=TRUE, minq=20, numbq=2)){
  if (filter_settings$rmlow)
  {
    i_rmlow = 1
  }
  else
  {
    i_rmlow = 0
  }
  if (filter_settings$rmN)
  {
    i_rmN = 1
  }
  else
  {
    i_rmN = 0
  }

  if (!is.null(r2)){
    if(!file.exists(r1)){stop("read1 fastq file does not exists.")}
    if(!file.exists(r2)){stop("read2 fastq file does not exists.")}
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
  else
  {
    stop("not implemented.")
  }
}


#' sc_exon_mapping
#' 
#' take the alinged read and map them to exons. the result will be written into optional fields
#' in bam file with different tags. for more information regarding to bam file format:
#' http://samtools.github.io/hts-specs/SAMv1.pdf
#' @name sc_exon_mapping
#' @param inbam input bam file. this is the output of alignment program
#' @param outbam output bam file with gene and barcode tag
#' @param annofn genome annotation, gff file. can have multiple files
#' @param am mapping status tag (default: YE)
#' @param ge gene id tag (default: GE)
#' @param bc cell barcode tag (default: YC)
#' @param mb molecular barcode tag (default: YM)
#' @param bc_len total barcode length
#' @param UMI_len UMI length
#' @param stnd perform strand specific mapping or not
#' @param fix_chr add `chr` to chromosome names, fix inconsistant names.
#'
#' @export
#'
#' @examples
#' #TODO
sc_exon_mapping <- function(inbam, outbam, annofn,
                            am="YE", ge="GE", bc="YC", mb="YM",
                            bc_len=8, UMI_len=6, stnd=TRUE, fix_chr=FALSE){
  if (stnd)
  {
    i_stnd = 1
  }
  else
  {
    i_stnd = 0
  }
  if (fix_chr)
  {
    i_fix_chr = 1
  }
  else
  {
    i_fix_chr = 0
  }

  if(!file.exists(inbam)){stop("input bam file does not exists.")}
  if(!file.exists(annofn)){stop("genome annotation file does not exists.")}

  rcpp_sc_exon_mapping(inbam, outbam, annofn, am, ge, bc, mb, bc_len, UMI_len, stnd, fix_chr)
}



#' sc_demultiplex
#' 
#' separate bam file by cell, output to outdir/count/[cell_id].csv.
#' the output contains information for all reads that can be mapped to exons.
#' including the gene id, UMI of that read and the distance to transcript end position.
#'
#' @param inbam input bam file. this should be the output of \code{sc_exon_mapping}
#' @param outdir output folder
#' @param bc_anno barcode annotation, first column is cell id, second column is cell barcode sequence
#' @param max_mis maximum mismatch allowed in barcode. default to be 1
#' @param am mapping status tag (default: YE)
#' @param ge gene id tag (default: GE)
#' @param bc cell barcode tag (default: YC)
#' @param mb molecular barcode tag (default: YM)
#' @param mito mitochondrial chromosome name, should be consistant with the chromosome name in bam file
#'
#' @export
#'
#' @examples
#' #TODO
sc_demultiplex <- function(inbam, outdir, bc_anno,
                          max_mis=1,
                          am="YE", ge="GE",
                          bc="YC",
                          mb="YM",
                          mito="MT"){
  dir.create(file.path(outdir, "count"), showWarnings = FALSE)
  dir.create(file.path(outdir, "stat"), showWarnings = FALSE)
  if(!file.exists(inbam)){stop("input bam file does not exists.")}
  if(!file.exists(bc_anno)){stop("barcode annotation file does not exists.")}
  rcpp_sc_demultiplex(inbam, outdir, bc_anno, max_mis, am, ge, bc, mb, mito)
}



#' sc_gene_counting
#' 
#' merge UMI and generate gene counting matrix
#'
#' @param outdir output folder that contains \code{sc_demultiplex} output
#' @param bc_anno barcode annotation, first column is cell id, second column is cell barcode sequence
#' @param UMI_cor correct UMI sequence error: 0 means no correction, 1 means simple correction and merge UMI with distance 1.
#' @param gene_fl whether to remove low abundant gene count. low abundant is defined as only one copy of one UMI for this gene
#'
#' @export
#'
#' @examples
#' #TODO
sc_gene_counting <- function(outdir, bc_anno, UMI_cor=1, gene_fl=FALSE){
  dir.create(file.path(outdir, "count"), showWarnings = FALSE)
  dir.create(file.path(outdir, "stat"), showWarnings = FALSE)
  if (gene_fl)
  {
    i_gene_fl = 1
  }
  else
  {
    i_gene_fl = 0
  }
  if(!file.exists(bc_anno)){stop("barcode annotation file does not exists.")}
  rcpp_sc_gene_counting(outdir, bc_anno, UMI_cor, i_gene_fl)

}



#' sc_celseq2_simulator
#'
#' @param r1fn simulated read1 fastq file
#' @param r2fn simulated read2 fastq file
#' @param annofn gff3 genome annotation
#' @param bc_anno barcode annotation
#' @param fafn genome fasta file
#' @param UMI_len UMI length
#' @param r_len read length
#' @param frag_mean the average fragment length, for CEL-seq2 it is usually 200~500
#' @param dup_mean the average PCR duplicate number
#' @param ran_dist the distribution used to generate random gene count
#' @param param a list of parameter of random distribution
#' @param seed random seed, used for reproducing data
#'
#' @export
#'
#' @examples
#' #TODO
sc_celseq2_simulator <- function(r1fn, r2fn, annofn, bc_anno, fafn,
                                UMI_len=6,
                                r_len=75,
                                frag_mean=300,
                                dup_mean=5,
                                ran_dist="gamma",
                                param=list(shape=0.2, scale=50),
                                seed=NA){

  if (ran_dist=="gamma"){
    r_dist = "gamma_random"
    pam = c(param$shape, param$scale)
  }
  else
  {
    stop("not implemented.") #TODO
  }
  if (is.na(seed)){
    rcpp_generate_celseq2_data(r1fn, r2fn, annofn, bc_anno, fafn, UMI_len, r_len, frag_mean, dup_mean, r_dist, pam, 0)
  }
  else{
    rcpp_generate_celseq2_data(r1fn, r2fn, annofn, bc_anno, fafn, UMI_len, r_len, frag_mean, dup_mean, r_dist, pam, seed)
  }

}



#' sc_detect_bc
#' 
#' detect cell barcode and generate the barcode annotation
#'
#' @param infq input fastq file, shoule be the output file of \code{sc_trim_barcode}
#' @param outcsv output barcode annotation
#' @param surfix the surfix of cell name, default to be `CELL_`, the cell name will be CELL_001, CELL_002 accordingly
#' @param bc_len the length of cell barcode, should be consistent with bl1+bl2 in \code{sc_trim_barcode}
#' @param max_reads the maximum of reads processed, default is 1000,000, set to "all" to process all reads (may spend more time)
#' @param min_count minimum counts to keep, barcode will be discarded if it has lower count. default is 10. this should be set according to \code{max_reads}
#' @param max_mismatch the maximum mismatch allowed, barcodes within this number will be considered as sequence error and merged, default is 1.
#' @export
#'
#' @examples
#' #TODO
sc_detect_bc <- function(infq, outcsv, surfix="CELL_", bc_len, max_reads=1000000, min_count=10, max_mismatch=1){

  if(!file.exists(infq)){stop("input fastq file does not exists.")}
  if(max_reads=="all"){m_r = -1}
  else{m_r=max_reads}
  rcpp_sc_detect_bc(infq, outcsv, surfix, bc_len, m_r, min_count, max_mismatch)
  
}

