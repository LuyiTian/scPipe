#' @export
test_check <- function(fastq, barcodes, bstart, blength, tmp_text) {
  # this function is a test at the moment, not real representation of how it will
  # work. This is just to test the c code.
  # write the barcodes as a text file
  barcodes <- read.csv(barcodes, header=FALSE)
  write(barcodes$V2, file=tmp_text, ncolumns=1, sep="\n")


  check_barcode_reads(fastq, tmp_text, bstart, blength)

  invisible()
}
