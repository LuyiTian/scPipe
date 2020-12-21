# @export
#test_check <- function(fastq, barcodes, bstart, blength, search_lines) {
  # valid code for running our C++ program
#  barcodes <- read.csv(barcodes, header=FALSE, strip.white=TRUE)
#  tmp_barcode_file = paste0(getwd(), "/tmp_barcodes.txt")
#  on.exit({if(file.exists(tmp_barcode_file)) file.remove(tmp_barcode_file)}, add=TRUE)

#  write(barcodes$V2, file=tmp_barcode_file, ncolumns=1, sep="\n")


#  contin = check_barcode_reads(fastq, tmp_barcode_file, bstart, blength, search_lines)

#  return (contin)
#}

#' Check Valid bstart
#'
#' @param fastq file containing reads
#' @param barcodes csv file
#' @param bstart the start position for barcodes in the given reads
#' @param blength length of each barcode
#'
#' @return Boolean; TRUE if program can continue execution, FALSE otherwise.
#'
check_barcode_start_position <- function(fastq, barcodes, bstart, blength, search_lines) {
  # read in barcodes to desired format and write to file.
  # we need to be able to handle multiple barcode files
  read_barcodes = NULL
  for (bc_file in barcodes) {
    if (is.null(read_barcodes)) read_barcodes <- read.csv(bc_file, header=FALSE, strip.white=TRUE)
    else read_barcodes <- c(read_barcodes, read.csv(bc_file, header=FALSE, strip.white=TRUE))
  }
  #barcodes <- read.csv(barcodes, header=FALSE, strip.white=TRUE)
  tmp_barcode_file = paste0(getwd(), "/tmp_barcodes.txt")
  on.exit({if(file.exists(tmp_barcode_file)) file.remove(tmp_barcode_file)}, add=TRUE)
  # combine and write all barcode files to the one file
  write(read_barcodes$V2, file=tmp_barcode_file, ncolumns=1, sep="\n")

  # check if the given bstart param is valid
  continue = check_barcode_reads(fastq, tmp_barcode_file, bstart, blength, search_lines)

  return (continue)
}
