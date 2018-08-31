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
#' @export
#'
#' @examples
#' ens_chrY <- anno_import(system.file("extdata", "ensembl_hg38_chrY.gtf.gz", package = "scPipe"))
#'
anno_import <- function(filename) {
    anno <- rtracklayer::import(filename)

    # check if gene_id column is present, this is missing from RefSeq annotataions
    no_gene_ids <- is.null(anno$gene_id)
    has_dbx <- !is.null(anno$Dbxref)

    if (no_gene_ids && has_dbx) {
        anno <- infer_gene_id_from_dbx(anno)
        return(anno_to_saf(anno))
    }

    # check if every entry has a gene_id value, this is not true for ENSEMBL gff3
    incomplete_gene_ids <- anyNA(anno$gene_id[anno$type == "exon"])

    if (incomplete_gene_ids) {
        anno <- infer_gene_id_from_parent(anno)
        return(anno_to_saf(anno))
    }

    # gene_id column is present and contains necessary information
    # return SAF converted data.frame
    return(anno_to_saf(anno))
}

# convert annotation from GenomicRanges to Simple Annotation Format(SAF)
# http://bioinf.wehi.edu.au/featureCounts/ for details on SAF
anno_to_saf <- function(full_anno) {
    full_anno %>%
        as.data.frame() %>%
        dplyr::filter(type == "exon") %>%
        dplyr::select(gene_id, seqnames, start, end, strand) %>%
        dplyr::rename(
            GeneID = gene_id,
            Chr = seqnames,
            Start = start,
            End = end,
            Strand = strand
        ) %>%
        dplyr::select(GeneID, dplyr::everything())
}

# shorthand for stringr string interpolation
fmt_str <- stringr::str_interp

infer_gene_id_from_dbx <- function(anno) {
    extract_gene_id <- function(x) {
        sapply(x, function(x) x[1]) %>%
            stringr::str_extract("GeneID:[^,]+") %>%
            stringr::str_remove("GeneID:")
    }

    gene_ids <- extract_gene_id(anno$Dbxref)

    anno$gene_id <- gene_ids

    anno
}

infer_gene_id_from_parent <- function(anno) {
    transcript_hash <- local({
        transcripts <- anno %>%
            as.data.frame() %>%
            filter(!is.na(transcript_id)) %>%
            select(transcript_id, Parent) %>%
            mutate(
                transcript_id = paste0("transcript:", transcript_id),
                Parent = stringr::str_remove(Parent, "gene:")
            )

        if (any(is.na(transcripts$Parent))) {
            warning("there are entries in annotation with transcript_id but no parent")
        }

        with(transcripts, hashmap::hashmap(transcript_id, Parent))
    })

    gene_ids <- anno %>%
        as.data.frame() %>%
        dplyr::filter(type == "exon") %>%
        dplyr::select(Parent) %>%
        dplyr::mutate(gene_id = transcript_hash[[unlist(Parent)]]) %>%
        dplyr::pull(gene_id)

    anno$gene_id[anno$type == "exon"] <- gene_ids

    anno
}
