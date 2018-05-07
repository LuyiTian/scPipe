#' scPipe - single cell RNA-seq pipeline
#'
#' @description The scPipe will do cell barcode demultiplexing, UMI deduplication and quality control on
#' fastq data generated from all protocols
#'
#' @author Luyi Tian <tian.l@wehi.edu.au>; Shian Su <su.s@wehi.edu.au>
#' @docType package
#' @name scPipe
#' @import Rhtslib
#' @import SingleCellExperiment
#' @import Rcpp
#' @useDynLib scPipe, .registration = TRUE
#' @aliases scPipe scPipe-package
NULL
