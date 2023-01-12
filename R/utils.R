#' Get ERCC annotation table
#'
#' Helper function to retrieve ERCC annotation as a dataframe in SAF format
#'
#' @return data.frame containing ERCC annotation
#'
#' @examples
#' ercc_anno <- get_ercc_anno()
#'
#' @export
#'
get_ercc_anno <- function() {
    anno_import(system.file("extdata", "ERCC92_anno.gff3", package = "scPipe"))
}
