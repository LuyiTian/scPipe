#' Check Valid Barcode Start Position
#'
#' @description Checks to see if the given barcode start position (\code{bstart})
#' is valid for the fastq file. If the found barcode percentage is less than
#' the given \code{threshold}, a new barcode start position is searched for by
#' checking every position from the start of each read to 10 bases after the \code{bstart}
#' @param fastq file containing reads
#' @param barcode_file csv file
#' @param barcode_file_realname the real name of the csv file
#' @param bstart the start position for barcodes in the given reads
#' @param blength length of each barcode
#' @param search_lines the number of fastq lines to use for the check
#' @param threshold the minimum percentage of found barcodes to accept for the
#' program to continue
#'
#' @return Boolean; TRUE if program can continue execution, FALSE otherwise.
#'
check_barcode_start_position <- function(fastq, barcode_file, barcode_file_realname, bstart, blength, search_lines, threshold) {
    # check if the given bstart param is valid
    continue <- check_barcode_reads(fastq, barcode_file, barcode_file_realname, bstart, blength, search_lines, threshold)

    return (continue)
}
