### all classes defined for the scPipe package


#' The "Single Cell Data" (SCData)  class
#'
#' S4 class and the main class used by scPipe to hold single cell expression
#' data, other phenotype data and qualiy control information.
#' SCData extends the basic Bioconductor ExpressionSet class.
#'
#' This class is initialized from the output dir
#'
#' Methods that operate on SCData objects constitute the basic scPipe workflow.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{logged}:}{Scalar of class \code{"logical"}, indicating whether
#'    or not the expression data in the `exprs` slot have been log2-transformed
#'    or not.}
#'    \item{\code{gene_id_type}}{the type of gene id. could be `ensembl_gene_id` or
#'    `external_gene_name` or other type which can be assessed by \code{listAttributes} from \code{biomaRt} package}
#'    \item{\code{logExprsOffset}:}{Scalar of class \code{"numeric"}, providing an offset
#'    applied to expression data in the `exprs` slot when undergoing log2-transformation
#'    to avoid trying to take logs of zero.}
#'    \item{\code{FACSData}:}{contains the index sorted FACSData, rows are cell and columns
#'    are each cell surface markers}
#'    \item{\code{reducedExprDimension}:}{Matrix of class \code{"numeric"}, containing
#'    reduced-dimension coordinates for gene expression of single cell, where each row
#'    represent a cell and each column represent a dimension}
#'    \item{\code{reducedFACSDimension}:}{Matrix of class \code{"numeric"}, containing
#'    reduced-dimension coordinates for surface marker phenotypes of single cell.}
#'    \item{\code{onesense}:}{Matrix of class \code{"numeric"}, containing
#'    onesense tSNE coordinates for both phenotype and transcriptome.}
#'    \item{\code{QualityControlInfo}:}{Data frame of class
#'    \code{"AnnotatedDataFrame"} that can contain QC metrics used for in QC.}
#'    \item{\code{useForExprs}:}{Character string (one of 'exprs','tpm','counts' or 'fpkm') indicating
#'    which expression representation both internal methods and external packages should use.
#'    Defaults to 'exprs'.}
#'}
#' @name SCData
#' @rdname SCData
#' @inheritParams Biobase ExpressionSet
#' @aliases SCData-class
#' @references  the SCData class is adapted from SCESet from scater
#' (github.com/davismcc/scater/). Thank Davis for creating such a wonderful package
#' @exportClass SCData
#' @examples
#' ## list all dataset. the names can be used in organism attribute
#' library(biomaRt)
#' ensembl=useMart("ensembl")
#' listDatasets(ensembl)
#'
setClass("SCData",
         contains = "ExpressionSet",
         slots = c(logged = "logical",
                   gene_id_type = "character",
                   logExprsOffset = "numeric",
                   FACSData = "AnnotatedDataFrame",
                   reducedExprDimension = "matrix",
                   reducedFACSDimension = "matrix",
                   organism="character",
                   onesense = "matrix",
                   QualityControlInfo = "AnnotatedDataFrame",
                   useForExprs = "character"),
         prototype = prototype(new("VersionedBiobase",
                                   versions = c(classVersion("ExpressionSet"),
                                                 SCData = "0.99"))))


