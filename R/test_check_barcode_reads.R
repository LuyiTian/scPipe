#' @export
test_check <- function(fastq, barcodes, bstart, blength, search_lines) {
  # valid code for running our C++ program
  barcodes <- read.csv(barcodes, header=FALSE, strip.white=TRUE)
  tmp_barcode_file = paste0(getwd(), "/tmp_barcodes.txt")
  on.exit({if(file.exists(tmp_barcode_file)) file.remove(tmp_barcode_file)}, add=TRUE)

  write(barcodes$V2, file=tmp_barcode_file, ncolumns=1, sep="\n")


  contin = check_barcode_reads(fastq, tmp_barcode_file, bstart, blength, search_lines)

  return (contin)
}
