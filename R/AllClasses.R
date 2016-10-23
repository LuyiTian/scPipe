### all classes defined for the scPipe package

################################################################################
### defining the SCData class

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
#'    \item{\code{logExprsOffset}:}{Scalar of class \code{"numeric"}, providing an offset 
#'    applied to expression data in the `exprs` slot when undergoing log2-transformation
#'    to avoid trying to take logs of zero.}
#'    \item{\code{reducedExprDimension}:}{Matrix of class \code{"numeric"}, containing
#'    reduced-dimension coordinates for gene expression of single cell.}
#'    \item{\code{reducedPhenDimension}:}{Matrix of class \code{"numeric"}, containing
#'    reduced-dimension coordinates for phenotype of single cell.}
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
#' @import Biobase
#' @aliases SCData-class
#' @references  Thanks to the scater package 
#' (github.com/davismcc/scater/) for their SCESet class, 
#' which provided the inspiration and template for SCData.
#' @exportClass SCData
setClass("SCData",
         contains = "ExpressionSet",
         slots = c(logged = "logical",
                   logExprsOffset = "numeric",
                   reducedExprDimension = "matrix",
                   reducedPhenDimension = "matrix",
                   onesense = "matrix",
                   QualityControlInfo = "AnnotatedDataFrame",
                   useForExprs = "character"),
         prototype = prototype(new("VersionedBiobase",
                                   versions = c(classVersion("ExpressionSet"),
                                                SCData = "0.99")))
)