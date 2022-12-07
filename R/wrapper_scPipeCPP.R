#' sc_trim_barcode
#'
#' @description Reformat fastq files so barcode and UMI sequences are moved from
#'   the sequence into the read name.
#'
#'@details Positions used in this function are 0-indexed, so they start from 0
#'  rather than 1. The default read structure in this function represents
#'  CEL-seq paired-ended reads. This contains a transcript in the first read, a
#'  UMI in the first 6bp of the second read followed by a 8bp barcode. So the
#'  read structure will be : \code{list(bs1=-1, bl1=0, bs2=6, bl2=8, us=0,
#'  ul=6)}. \code{bs1=-1, bl1=0} indicates negative start position and zero
#'  length for the barcode on read one, this is used to denote "no barcode" on
#'  read one. \code{bs2=6, bl2=8} indicates there is a barcode in read two that
#'  starts at the 7th base with length 8bp. \code{us=0, ul=6} indicates a UMI
#'  from first base of read two and the length in 6bp.
#'
#'  For a typical Drop-seq experiment the read structure will be
#'  \code{list(bs1=-1, bl1=0, bs2=0, bl2=12, us=12, ul=8)}, which means the read
#'  one only contains transcript, the first 12bp in read two are cell barcode, followed
#'  by a 8bp UMI.
#'
#' @name sc_trim_barcode
#' @param outfq the output fastq file, which reformat the barcode and UMI into
#'   the read name. Files ending in \code{.gz} will be automatically compressed.
#' @param r1 read one for pair-end reads. This read should contain
#'   the transcript.
#' @param r2 read two for pair-end reads, NULL if single read.
#'   (default: NULL)
#' @param read_structure a list containing the read structure configuration:
#'   \itemize{
#'     \item{bs1}: starting position of barcode in read one. -1 if no barcode in
#'       read one.
#'     \item{bl1}: length of barcode in read one, if there is no
#'       barcode in read one this number is used for trimming beginning of read
#'       one.
#'     \item{bs2}: starting position of barcode in read two
#'     \item{bl2}: length of barcode in read two
#'     \item{us}: starting position of UMI
#'     \item{ul}: length of UMI
#'   }
#' @param filter_settings A list contains read filter settings:\itemize{
#'  \item{rmlow} whether to remove the low quality reads.
#'  \item{rmN} whether to remove reads that contains N in UMI or cell barcode.
#'  \item{minq} the minimum base pair quality that we allowed
#'  \item{numbq} the maximum number of base pair that have quality
#'  below \code{numbq}
#'  }
#' @export
#' @return generates a trimmed fastq file named \code{outfq}
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
sc_trim_barcode <- function(outfq, r1, r2=NULL,
                            read_structure = list(
                                bs1=-1, bl1=0, bs2=6, bl2=8, us=0, ul=6),
                            filter_settings = list(
                                rmlow=TRUE, rmN=TRUE, minq=20, numbq=2)) {

    outdir <- regmatches(outfq, regexpr(".*/", outfq))
    if (outdir != character(0) && !dir.exists(outdir))
        dir.create(outdir, recursive = TRUE)

    if (filter_settings$rmlow) {
        i_rmlow <- 1
    }
    else {
        i_rmlow <- 0
    }
    if (filter_settings$rmN) {
        i_rmN <- 1
    }
    else {
        i_rmN <- 0
    }

    if (substr(outfq, nchar(outfq) - 2, nchar(outfq)) == ".gz") {
        write_gz <- TRUE
    }
    else {
        write_gz <- FALSE
    }

    if (!is.null(r2)) {
        if (!file.exists(r1)) {stop("read1 fastq file does not exists.")}
        if (!file.exists(r2)) {stop("read2 fastq file does not exists.")}

        # expand tilde to home path for downstream gzopen() call
        r1 <- path.expand(r1)
        r2 <- path.expand(r2)


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
                                    filter_settings$numbq,
                                    write_gz)
    }
    else {
        stop("not implemented.")
    }
}


#' sc_exon_mapping
#'
#' @description Map aligned reads to exon annotation.
#' The result will be written into optional fields in bam file with different
#' tags. Following this link for more information regarding to bam file format:
#' http://samtools.github.io/hts-specs
#'
#' The function can accept multiple bam file as input, if multiple bam file is
#' provided and the `bc_len` is zero, then the function will use the barcode in
#' the `barcode_vector` to insert into the `bc` bam tag. So the length of
#' `barcode_vector` and the length of `inbam` should be the same
#' If this is the case then the `max_mis` argument in `sc_demultiplex`
#' should be zero. If `be_len` is larger than zero, then the function will
#' still seek for barcode in fastq headers with given length. In this case
#' each bam file is not treated as from a single cell.
#'
#'
#' @name sc_exon_mapping
#' @param inbam input aligned bam file. can have multiple files as input
#' @param outbam output bam filename
#' @param annofn single string or vector of gff3 annotation filenames,
#'   data.frame in SAF format or GRanges object containing complete gene_id
#'   metadata column.
#' @param bam_tags list defining BAM tags where mapping information is
#'   stored.
#'   \itemize{
#'     \item "am": mapping status tag
#'     \item "ge": gene id
#'     \item "bc": cell barcode tag
#'     \item "mb": molecular barcode tag
#'   }
#' @param bc_len total barcode length
#' @param barcode_vector a list of barcode if each individual bam is a single
#' cell. (default: NULL). The barcode should be of the same length for each
#' cell.
#' @param UMI_len UMI length
#' @param stnd TRUE to perform strand specific mapping. (default: TRUE)
#' @param fix_chr TRUE to add `chr` to chromosome names, MT to chrM. (default: FALSE)
#' @param nthreads number of threads to use. (default: 1)
#'
#' @export
#' @return generates a bam file with exons assigned
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
sc_exon_mapping <- function(inbam, outbam, annofn,
                                bam_tags = list(am="YE", ge="GE", bc="BC", mb="OX"),
                                bc_len=8, barcode_vector="", UMI_len=6, stnd=TRUE, fix_chr=FALSE,
                                nthreads=1) {
    if (stnd) {
        i_stnd <- 1
    }
    else {
        i_stnd <- 0
    }
    if (fix_chr) {
        i_fix_chr <- 1
    }
    else {
        i_fix_chr <- 0
    }

    if (any(!file.exists(inbam))) {
        stop("At least one input bam file does not exist")
    } else {
        inbam <- path.expand(inbam)
    }

    outbam <- path.expand(outbam)

    # if (length(inbam) > 1) {
    #   stop("Only one bam file can be used as input")
    # }

    if (!is(annofn, "character") &&
        !is(annofn, "GRanges") &&
        !is(annofn, "data.frame")) {
        stop("'annofn' must be either character vector, GRanges, or data.frame object")
    }

    if (is(annofn, "character")) {
        if (any(!file.exists(annofn))) {
            stop("At least one genome annotation file does not exist")
        } else {
            annofn <- path.expand(annofn)
        }
        rcpp_sc_exon_mapping_df_anno(inbam, outbam, anno_import(annofn), bam_tags$am, bam_tags$ge, bam_tags$bc, bam_tags$mb, bc_len,
                                    barcode_vector, UMI_len, stnd, fix_chr, nthreads)
    } else if (is(annofn, "GRanges")) {
        rcpp_sc_exon_mapping_df_anno(inbam, outbam, anno_to_saf(annofn), bam_tags$am, bam_tags$ge, bam_tags$bc, bam_tags$mb, bc_len,
                                    barcode_vector, UMI_len, stnd, fix_chr, nthreads)
    } else if (is(annofn, "data.frame")) {
        validate_saf(annofn)
        rcpp_sc_exon_mapping_df_anno(inbam, outbam, annofn, bam_tags$am, bam_tags$ge, bam_tags$bc, bam_tags$mb, bc_len,
                                    barcode_vector, UMI_len, stnd, fix_chr, nthreads)
    }
}


#' sc_demultiplex
#'
#' @description Process bam file by cell barcode,
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
#' @param max_mis maximum mismatch allowed in barcode. (default: 1)
#' @param bam_tags list defining BAM tags where mapping information is
#'   stored.
#'   \itemize{
#'     \item "am": mapping status tag
#'     \item "ge": gene id
#'     \item "bc": cell barcode tag
#'     \item "mb": molecular barcode tag
#'   }
#' @param mito mitochondrial chromosome name.
#' This should be consistent with the chromosome names in the bam file.
#' @param has_UMI whether the protocol contains UMI (default: TRUE)
#' @param nthreads number of threads to use. (default: 1)
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
sc_demultiplex <- function(inbam, outdir, bc_anno,
                            max_mis=1,
                            bam_tags = list(am="YE", ge="GE", bc="BC", mb="OX"),
                            mito="MT",
                            has_UMI=TRUE,
                            nthreads = 1) {
    dir.create(file.path(outdir, "count"), showWarnings = FALSE)
    dir.create(file.path(outdir, "stat"), showWarnings = FALSE)

    if (!dir.exists(outdir))
        dir.create(outdir, recursive = TRUE)

    outdir <- path.expand(outdir)

    if (!file.exists(inbam)) {
        stop("input bam file does not exists.")
    } else {
        inbam <- path.expand(inbam)
    }
    if (!file.exists(bc_anno)) {
        stop("barcode annotation file does not exists.")
    } else {
        bc_anno <- path.expand(bc_anno)
    }
    rcpp_sc_demultiplex(inbam, outdir, bc_anno, max_mis,
                        bam_tags$am, bam_tags$ge, bam_tags$bc, bam_tags$mb,
                        mito, has_UMI, nthreads)
}


#' sc_correct_bam_bc
#'
#' @description update the cell barcode tag in bam file with corrected barcode
#' output to a new bam file. the function will be useful for methods
#' that use the cell barcode information from bam file, such as `Demuxlet`
#'
#' @param inbam input bam file. This should be the output of
#' \code{sc_exon_mapping}
#' @param outbam output bam file with updated cell barcode
#' @param bc_anno barcode annotation, first column is cell id, second column
#' is cell barcode sequence
#' @param max_mis maximum mismatch allowed in barcode. (default: 1)
#' @param bam_tags list defining BAM tags where mapping information is
#'   stored.
#'   \itemize{
#'     \item "am": mapping status tag
#'     \item "ge": gene id
#'     \item "bc": cell barcode tag
#'     \item "mb": molecular barcode tag
#'   }
#' @param mito mitochondrial chromosome name.
#' This should be consistent with the chromosome names in the bam file.
#' @param nthreads number of threads to use. (default: 1)
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
#' sc_correct_bam_bc(file.path(data_dir, "out.map.bam"),
#'     file.path(data_dir, "out.map.clean.bam"),
#'     barcode_annotation_fn)
#' ...
#' }
#'
sc_correct_bam_bc <- function(inbam, outbam, bc_anno,
                            max_mis=1,
                            bam_tags = list(am="YE", ge="GE", bc="BC", mb="OX"),
                            mito="MT",
                            nthreads = 1) {

    if (!file.exists(inbam)) {
        stop("input bam file does not exists.")
    } else {
        inbam <- path.expand(inbam)
    }

    outbam <- path.expand(outbam)

    if (!file.exists(bc_anno)) {
        stop("barcode annotation file does not exists.")
    } else {
        bc_anno <- path.expand(bc_anno)
    }
    rcpp_sc_clean_bam(inbam, outbam, bc_anno, max_mis,
                        bam_tags$am, bam_tags$ge, bam_tags$bc, bam_tags$mb,
                        mito, nthreads)
}


#' sc_gene_counting
#'
#' @description Generate gene counts matrix with UMI deduplication
#'
#' @param outdir output folder containing \code{sc_demultiplex} output
#' @param bc_anno barcode annotation comma-separated-values, first column is
#'   cell id, second column is cell barcode sequence
#' @param UMI_cor correct UMI sequencing error: 0 means no correction, 1 means
#'   simple correction and merge UMI with distance 1. 2 means merge on both UMI
#'   alignment position match.
#' @param gene_fl whether to remove low abundance genes. A gene is considered to
#'   have low abundance if only one copy of one UMI is associated with it.
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
sc_gene_counting <- function(outdir, bc_anno, UMI_cor=2, gene_fl=FALSE) {
    if (!dir.exists(outdir))
        dir.create(outdir, recursive = TRUE)

    outdir <- path.expand(outdir)

    dir.create(file.path(outdir, "count"), showWarnings = FALSE)
    dir.create(file.path(outdir, "stat"), showWarnings = FALSE)

    if (gene_fl) {
        i_gene_fl <- 1
    } else {
        i_gene_fl <- 0
    }

    if (!file.exists(bc_anno)) {
        stop("barcode annotation file does not exists.")
    } else {
        bc_anno <- path.expand(bc_anno)
    }

    rcpp_sc_gene_counting(outdir, bc_anno, UMI_cor, i_gene_fl)
}



#' sc_detect_bc
#'
#' @description Detect cell barcode and generate the barcode annotation
#'
#' @param infq input fastq file, should be the output file of
#' \code{sc_trim_barcode}
#' @param outcsv output barcode annotation
#' @param prefix the prefix of cell name (default: `CELL_`)
#' @param bc_len the length of cell barcode, should be consistent with bl1+bl2
#' in \code{sc_trim_barcode}
#' @param max_reads the maximum of reads processed (default: 1,000,000)
#' @param min_count minimum counts to keep, barcode will be discarded if
#' it has lower count. Default value is 10. This should be set according
#' to \code{max_reads}.
#' @param number_of_cells number of cells kept in result. (default: 10000)
#' @param max_mismatch the maximum mismatch allowed. Barcodes within this
#' number will be considered as sequence error and merged. (default: 1)
#' @param white_list_file a file that list all the possible barcodes
#' each row is a barcode sequence. the list for 10x can be found at:
#' https://community.10xgenomics.com/t5/Data-Sharing/List-of-valid-cellular-barcodes/td-p/527
#' (default: NULL)
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
sc_detect_bc <- function(infq, outcsv, prefix="CELL_", bc_len,
                            max_reads=1000000, min_count=10, number_of_cells=10000,
                            max_mismatch=1,white_list_file=NULL) {
    if (!file.exists(infq)) {
        stop("input fastq file does not exists.")
    } else {
        infq <- path.expand(infq)
    }
    outcsv <- path.expand(outcsv)
    if (!is.null(white_list_file)){
        if(!file.exists(white_list_file)){
            stop("input whitelist file does not exists.")
        }
    } else {
        white_list_file <- ""
    }
    if (max_reads=="all") {m_r <- -1}
    else {m_r<-max_reads}
    rcpp_sc_detect_bc(infq, outcsv, prefix, bc_len, m_r, number_of_cells,
        min_count, max_mismatch, white_list_file)
}

#' sc_demultiplex_and_count
#'
#' @description Wrapper to run \code{\link{sc_demultiplex}} and
#'   \code{\link{sc_gene_counting}} with a single command
#'
#' @inheritParams sc_demultiplex
#' @inheritParams sc_gene_counting
#'
#' @return no return
#'
#' @examples
#' \dontrun{
#' refer to the vignettes for the complete workflow, replace demultiplex and
#' count with single command:
#' ...
#' sc_demultiplex_and_count(
#'    file.path(data_dir, "out.map.bam"),
#'    data_dir,
#'    barcode_annotation_fn,
#'    has_UMI = FALSE
#' )
#' ...
#' }
#'
#' @export
#'
sc_demultiplex_and_count <- function(
    inbam, outdir, bc_anno,
    max_mis = 1,
    bam_tags = list(am="YE", ge="GE", bc="BC", mb="OX"),
    mito = "MT", has_UMI = TRUE, UMI_cor = 1, gene_fl = FALSE,
    nthreads = 1) {
    sc_demultiplex(
        inbam = inbam,
        outdir = outdir,
        bc_anno = bc_anno,
        max_mis = max_mis,
        bam_tags = bam_tags,
        mito = mito,
        has_UMI = has_UMI,
        nthreads = nthreads
    )

    sc_gene_counting(
        outdir = outdir,
        bc_anno = bc_anno,
        UMI_cor = UMI_cor,
        gene_fl = gene_fl
    )
}

#' sc_count_aligned_bam
#'
#' @description Wrapper to run \code{\link{sc_exon_mapping}},
#'   \code{\link{sc_demultiplex}} and \code{\link{sc_gene_counting}} with a
#'   single command
#'
#' @inheritParams sc_exon_mapping
#' @inheritParams sc_demultiplex
#' @inheritParams sc_gene_counting
#' @param keep_mapped_bam TRUE if feature mapped bam file should be retained.
#'
#' @return no return
#'
#' @examples
#' \dontrun{
#' sc_count_aligned_bam(
#'   inbam = "aligned.bam",
#'   outbam = "mapped.bam",
#'   annofn = c("MusMusculus-GRCm38p4-UCSC.gff3", "ERCC92_anno.gff3"),
#'   outdir = "output",
#'   bc_anno = "barcodes.csv"
#' )
#' }
#'
#' @export
#'
sc_count_aligned_bam <- function(
    inbam, outbam, annofn,
    bam_tags = list(am="YE", ge="GE", bc="BC", mb="OX"),
    bc_len = 8, UMI_len = 6,
    stnd = TRUE, fix_chr = FALSE,
    outdir,
    bc_anno,
    max_mis = 1,
    mito = "MT",
    has_UMI = TRUE, UMI_cor = 1, gene_fl = FALSE,
    keep_mapped_bam = TRUE,
    nthreads = 1
    ) {
    sc_exon_mapping(
        inbam = inbam,
        outbam = outbam,
        annofn = annofn,
        bam_tags = bam_tags,
        bc_len = bc_len,
        UMI_len = UMI_len,
        stnd = stnd,
        fix_chr = fix_chr,
        nthreads = nthreads
    )

    sc_demultiplex_and_count(
        inbam = outbam,
        outdir = outdir,
        bc_anno = bc_anno,
        max_mis = max_mis,
        bam_tags = bam_tags,
        mito = mito,
        has_UMI = has_UMI,
        UMI_cor = UMI_cor,
        gene_fl = gene_fl,
        nthreads = nthreads
    )

    if (!keep_mapped_bam) {
        unlink(outbam)
    }
}
