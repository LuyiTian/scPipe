#' Check Valid Barcode Start Position
#'
#' @description Checks to see if the given barcode start position (\code{bstart})
#' is valid for the fastq file. If the found barcode percentage is less than
#' the given \code{threshold}, a new barcode start position is searched for by
#' checking every position from the start of each read to 10 bases after the \code{bstart}
#' @param fastq file containing reads
#' @param barcodes csv file
#' @param bstart the start position for barcodes in the given reads
#' @param blength length of each barcode
#' @param search_lines the number of fastq lines to use for the check
#' @param threshold the minimum percentage of found barcodes to accept for the
#' program to continue
#'
#' @return Boolean; TRUE if program can continue execution, FALSE otherwise.
#'
check_barcode_start_position <- function(fastq, barcodes, bstart, blength, search_lines, threshold) {
  # read in barcodes to desired format and write to file.
  # we need to be able to handle multiple barcode files
  read_barcodes = NULL
  for (bc_file in barcodes) {
    if (is.null(read_barcodes)) read_barcodes <- read.csv(bc_file, header=FALSE, strip.white=TRUE)
    else read_barcodes <- c(read_barcodes, read.csv(bc_file, header=FALSE, strip.white=TRUE))
  }
  tmp_barcode_file = paste0(getwd(), "/tmp_barcodes.txt")
  on.exit({if(file.exists(tmp_barcode_file)) file.remove(tmp_barcode_file)}, add=TRUE)
  # combine and write all barcode files to the one file
  write(read_barcodes$V2, file=tmp_barcode_file, ncolumns=1, sep="\n")

  # check if the given bstart param is valid
  continue = check_barcode_reads(fastq, tmp_barcode_file, bstart, blength, search_lines, threshold)

  return (continue)
}
