#' Get read structure for particular scRNA-seq protocol
#'
#' The supported protocols are:
#' \itemize{
#'  \item CelSeq
#'  \item CelSeq2
#'  \item DropSeq
#'  \item 10x (also called ChromiumV1)
#'}
#' If you know the structure of a specific protocol and would like it supported,
#' please leave a issue post at www.github.com/luyitian/scPipe.
#'
#' @param protocol name of the protocol
#'
#' @return list of UMI and Barcode locations for use in other scPipe functions
#' @export
#'
#' @examples
#' get_read_str("celseq")
get_read_str <- function(protocol) {
    original_protocol <- protocol
    error_msg <- paste("No read structure for", original_protocol, "found")

    # strip non-alphanumeric characters and convert to lower case
    protocol <- gsub("[^[:alnum:]_]", "", protocol)
    protocol <- tolower(protocol)

    switch(
        protocol,
        "celseq"         = list(bs1 = -1, bl1 = 0, bs2 = 6, bl2 = 8,  us = 0,  ul = 6),
        "celseq2"        = list(bs1 = -1, bl1 = 0, bs2 = 6, bl2 = 8,  us = 0,  ul = 6),
        "dropseq"        = list(bs1 = -1, bl1 = 0, bs2 = 0, bl2 = 12, us = 12, ul = 8),
        "10x"            = list(bs1 = -1, bl1 = 0, bs2 = 0, bl2 = 16, us = 16, ul = 10),
        "chromiumv1"     = list(bs1 = -1, bl1 = 0, bs2 = 0, bl2 = 16, us = 16, ul = 10),
        stop(error_msg)
    )
}
