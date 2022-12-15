#####################################################################
# Generating Fragments for Aligned and Demultiplexed scATAC-Seq Reads
#####################################################################
#' @name sc_atac_create_fragments
#' @title Generating the popular fragments for scATAC-Seq data
#' @description Takes in a tagged and sorted BAM file and outputs the associated fragments in a .bed file
#'
#' @param inbam The tagged, sorted and duplicate-free input BAM file
#' @param output_folder The path of the output folder
#' @param min_mapq : int
#'    Minimum MAPQ to retain fragment
#' @param nproc : int, optional
#'    Number of processors to use. Default is 1.
#' @param cellbarcode : str
#'   Tag used for cell barcode. Default is CB (used by cellranger)
#' @param chromosomes : str, optional
#'    Regular expression used to match chromosome names to include in the
#'    output file. Default is "(?i)^chr" (starts with "chr", case-insensitive).
#'    If None, use all chromosomes in the BAM file.
#' @param readname_barcode : str, optional
#'    Regular expression used to match cell barocde stored in read name.
#'    If None (default), use read tags instead. Use "[^:]*" to match all characters
#'    before the first colon (":").
#' @param cells : str
#'    File containing list of cell barcodes to retain. If None (default), use all cell barcodes
#'    found in the BAM file.
#' @param max_distance : int, optional
#'    Maximum distance between integration sites for the fragment to be retained.
#'    Allows filtering of implausible fragments that likely result from incorrect
#'    mapping positions. Default is 5000 bp.
#' @param min_distance : int, optional
#'    Minimum distance between integration sites for the fragment to be retained.
#'    Allows filtering implausible fragments that likely result from incorrect
#'    mapping positions. Default is 10 bp.
#' @param chunksize : int
#'    Number of BAM entries to read through before collapsing and writing
#'    fragments to disk. Higher chunksize will use more memory but will be
#'    faster.
#'
#' @returns returns NULL
#'
#' @export
sc_atac_create_fragments <- function(
        inbam,
        output_folder="",
        min_mapq=30,
        nproc=1,
        cellbarcode="CB",
        chromosomes="^chr",
        readname_barcode=NULL,
        cells=NULL,
        max_distance=5000,
        min_distance=10,
        chunksize=500000) {
    if(output_folder == ''){
        output_folder <- file.path(getwd(), "scPipe-atac-output")
    }

    if (!dir.exists(output_folder)){
        dir.create(output_folder,recursive=TRUE)
        message("Output Directory Does Not Exist. Created Directory: ", output_folder)
    }

    log_and_stats_folder <- file.path(output_folder, "scPipe_atac_stats")
    dir.create(log_and_stats_folder, showWarnings = FALSE)
    log_file <- file.path(log_and_stats_folder, "log_file.txt")
    if(!file.exists(log_file)) file.create(log_file)
    write(
    c(
        "sc_atac_create_fragments starts at ",
        as.character(Sys.time()),
        "\n"
    ),
    file = log_file, append = TRUE)

    # output = paste0(output_folder, "/fragments.bed")

    chrom <- get_chromosomes(inbam, keep_contigs=chromosomes) # List with names as contigs and elements as lengths
    cells <- read_cells(cells) # character vector / StringVector

    sc_atac_create_fragments_cpp(inbam,
                                output_folder,
                                names(chrom),
                                as.integer(chrom),
                                min_mapq,
                                nproc,
                                cellbarcode,
                                chromosomes,
                                readname_barcode,
                                cells,
                                max_distance,
                                min_distance,
                                chunksize)

    invisible()
}

#' Get Chromosomes
#'
#' Gets a list of NamedList of chromosomes and the reference length
#' acquired through the bam index file.
#'
#' @param bam file path to the bam file to get data from
#' @param keep_contigs regular expression used with grepl to filter reference names
#'
#' @return a named list where element names are chromosomes reference names and elements are integer lengths
get_chromosomes <- function(bam, keep_contigs="^chr") {
    if (is.null(keep_contigs)) {
        keep_contigs <- "."
    }
    if (!file.exists(paste0(bam, ".bai"))) {
        Rsamtools::indexBam(bam)
    }

    idxstats <- Rsamtools::idxstatsBam(bam)
    # regex = (?i)^chr
    # above matches chr at start of string ignoring case ((?i) means ignore case)
    # R version: grepl("^chr", ignore.case=TRUE)
    contigs <- grepl(keep_contigs, idxstats$seqnames, ignore.case=TRUE)
    mapped <- idxstats$mapped > 0
    keep_contigs <- idxstats[contigs & mapped, ]
    conlen <- stats::setNames(as.list(keep_contigs$seqlength), keep_contigs$seqnames)

    return(conlen)
}

#' Read Cell barcode file
#'
#' @param cells the file path to the barcode file. Assumes
#' one barcode per line or barcode csv.
#' Or, cells can be a comma delimited string of barcodes
#' @return a character vector of the provided barcodes
#' @importFrom data.table fread
read_cells <- function(cells) {
    if (is.null(cells)) return(NULL)

    if (file.exists(cells)) {
        barcodes <- data.table::fread(cells, header=FALSE, data.table=FALSE)[[1]]
    } else {
        barcodes <- strsplit(cells, ",")[[1]]
    }

    return (barcodes)

}
