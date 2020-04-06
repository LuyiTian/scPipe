#' Import gene annotation
#'
#' Imports and GFF3 or GTF gene annotation file and transforms it into a SAF
#' formatted data.frame. SAF described at
#' http://bioinf.wehi.edu.au/featureCounts/. SAF contains positions for exons,
#' strand and the GeneID they are associated with.
#'
#'
#' @param filename The name of the annotation gff3 or gtf file. File can be
#'   gzipped.
#'
#' @description Because of the variations in data format depending on annotation
#'   source, this function has only been tested with human annotation from
#'   ENSEMBL, RefSeq and Gencode. If it behaves unexpectedly with any annotation
#'   please submit an issue at www.github.com/LuyiTian/scPipe with details.
#'
#' @return data.frame containing exon information in SAF format
#'
#' @importFrom tools file_ext
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' ens_chrY <- anno_import(system.file("extdata", "ensembl_hg38_chrY.gtf.gz", package = "scPipe"))
#'
anno_import <- function(filename) {
    accepted_formats <- c("gff", "gff3", "gtf")
    file_format <- tools::file_ext(stringr::str_remove(filename, ".gz$"))
    all_valid_format <- any(file_format %in% accepted_formats)

    if (!all_valid_format) {
        stop("only files the following annotation formats are accepted: ", paste(accepted_formats, collapse = ", "), " and their gzipped variants")
    }

    anno <- lapply(filename, rtracklayer::import)
    anno <- lapply(anno, anno_to_saf)

    # gene_id column is present and contains necessary information
    # return SAF converted data.frame
    return(do.call(rbind, anno))
}

#' Convert annotation from GenomicRanges to Simple Annotation Format (SAF)
#'
#' Convert a GRanges object containing type and gene_id information into a SAF
#' format data.frame. SAF described at http://bioinf.wehi.edu.au/featureCounts/.
#' SAF contains positions for exons, strand and the GeneID they are associated
#' with.
#'
#' @param anno The GRanges object containing exon information
#'
#' @description This function converts a GRanges object into a data.frame of the
#'   SAF format for scPipe's consumption. The GRanges object should contain a
#'   "type" column where at least some features are annotated as "exon", in
#'   addition there should be a gene_id column specifying the gene to which the
#'   exon belongs. In the SAF only the gene ID, chromosome, start, end and
#'   strand are recorded, this is a gene-exon centric format, with all entries
#'   containing the same gene ID treated as exons of that gene. It is possible
#'   to count alternative features by setting the gene_id column to an arbitrary
#'   feature name and having alternative features in the SAF table, the main
#'   caveat is that the features are still treated as exons, and the mapping
#'   statistics for exon and intron will not reflect biological exons and
#'   introns but rather the annotation features.
#'
#' @return data.frame containing exon information in SAF format
#'
#' @importFrom GenomicRanges mcols
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' anno <- system.file("extdata", "ensembl_hg38_chrY.gtf.gz", package = "scPipe")
#' saf_chrY <- anno_to_saf(rtracklayer::import(anno))
#' }
#'
anno_to_saf <- function(anno) {
    stopifnot(is(anno, "GRanges"))

    meta_cols <- colnames(GenomicRanges::mcols(anno))
    if (!"type" %in% meta_cols) {
        stop("'type' column missing from GRanges metadata")
    }

    anno <- infer_gene_ids(anno)
    meta_cols <- colnames(GenomicRanges::mcols(anno))

    if (!"gene_id" %in% meta_cols) {
        stop("'gene_id' column missing from GRanges metadata and could not be inferred")
    }

    anno_df <- as.data.frame(anno)

    n_exons <- anno_df %>%
        dplyr::filter(.data$type == "exon") %>%
        nrow()

    if (n_exons == 0) {
        stop("no exons found in annotation, must be labelled 'exon' under 'type' column")
    }

    saf <- anno_df %>%
        dplyr::filter(.data$type == "exon") %>%
        dplyr::select("gene_id", "seqnames", "start", "end", "strand") %>%
        dplyr::rename(
            GeneID = .data$gene_id,
            Chr = .data$seqnames,
            Start = .data$start,
            End = .data$end,
            Strand = .data$strand
        )

    if (anyNA(saf$GeneID)) {
        orig_rows <- nrow(saf)
        saf <- saf %>%
            dplyr::filter(!is.na(GeneID))
        filt_rows <- nrow(saf)
        message(glue::glue("NA found in GeneID of {orig_rows - filt_rows} of {orig_rows} entries, automatically removing these entries"))
    }

    saf %>% dplyr::select(.data$GeneID, dplyr::everything())
}

infer_gene_ids <- function(anno) {
    # check if gene_id column is present, this is missing from RefSeq annotations
    no_gene_ids <- is.null(anno$gene_id)
    has_dbx <- !is.null(anno$Dbxref)

    if (!no_gene_ids) {
        exons_missing_gene_ids <- anno %>%
            as.data.frame() %>%
            dplyr::filter(.data$type == "exon") %>%
            dplyr::pull("gene_id") %>%
            is.na() %>%
            any()

        if (!exons_missing_gene_ids) {
            return(anno)
        }
    }
    if (!no_gene_ids && all(!is.na(anno$gene_id))) {
        return(anno)
    }

    if (no_gene_ids && has_dbx) {
        anno <- infer_gene_id_from_dbx(anno)
        return(anno)
    }

    # check if every entry has a gene_id value, this is not true for ENSEMBL gff3
    incomplete_gene_ids <- anyNA(anno$gene_id[anno$type == "exon"])
    has_parent <- !is.null(anno$Parent)

    if ((no_gene_ids || incomplete_gene_ids) && has_parent) {
        anno <- infer_gene_id_from_parent(anno)
        return(anno)
    }

    if (incomplete_gene_ids) {
        stop("some exons have missing gene_id, and gene_id could not be inferred")
    }

    return(anno)
}

# shorthand for stringr string interpolation
fmt_str <- stringr::str_interp

infer_gene_id_from_dbx <- function(anno) {
    extract_gene_id <- function(anno) {
        sapply(anno, function(x) x[stringr::str_detect(x, "GeneID")][1]) %>%
            stringr::str_extract("GeneID:[^,]+") %>%
            stringr::str_remove("GeneID:")
    }

    gene_ids <- extract_gene_id(anno$Dbxref)

    anno$gene_id <- gene_ids

    anno
}

infer_gene_id_from_parent <- function(anno) {
    # create hash map from transcript id to parent gene id
    transcript_hash <- local({
        transcripts <- anno %>%
            as.data.frame() %>%
            dplyr::filter(!is.na(.data$transcript_id)) %>%
            dplyr::select("transcript_id", "Parent") %>%
            dplyr::mutate(
                transcript_id = paste0("transcript:", .data$transcript_id),
                Parent = stringr::str_remove(.data$Parent, "gene:")
            )

        if (any(is.na(transcripts$Parent))) {
            warning("there are entries in annotation with transcript_id but no parent")
        }

        hash::hash(transcripts$transcript_id, transcripts$Parent)
    })

    # get parent ids of all exons
    parents <- anno %>%
        as.data.frame() %>%
        dplyr::filter(type == "exon") %>%
        dplyr::pull("Parent") %>%
        unlist()
    
    gene_ids <- character(length(parents))
    
    for (i in seq_along(parents)) {
        gene_ids[i] <- transcript_hash[[parents[i]]]
    }

    anno$gene_id[anno$type == "exon"] <- gene_ids

    anno
}

validate_saf <- function(saf_df) {
    if (!is(saf_df, "data.frame")) {
        stop("annotation must object of class data.frame")
    }

    expected_colnames <- c("GeneID", "Chr", "Start", "End", "Strand")
    if (!identical(colnames(saf_df), expected_colnames)) {
        stop("columns of SAF data.frame must be: ", paste(expected_colnames, collapse = ", "))
    }

    if (anyNA(saf_df)) {
        stop("SAF data.frame must not contain any NA")
    }
}
