#' @name DimReducedExpr
#' @export
#' @docType methods
#' @rdname DimReducedExpr
setGeneric("DimReducedExpr", function(object) {
  standardGeneric("DimReducedExpr")
})

#' @name DimReducedExpr<-
#' @export
#' @docType methods
#' @rdname DimReducedExpr
setGeneric("DimReducedExpr<-", function(object, value) {
  standardGeneric("DimReducedExpr<-")
})

#' @name QCMetrics
#' @export
#' @docType methods
#' @rdname QCMetrics
setGeneric("QCMetrics", function(object) {
  standardGeneric("QCMetrics")
})

#' @name QCMetrics<-
#' @export
#' @docType methods
#' @rdname QCMetrics
setGeneric("QCMetrics<-", function(object, value) {
  standardGeneric("QCMetrics<-")
})

#' @name organism
#' @aliases organism
#' @export
#' @docType methods
#' @rdname organism
setGeneric("organism", function(object) {
  standardGeneric("organism")
})

#' @name organism<-
#' @aliases organism
#' @export
#' @docType methods
#' @rdname organism
setGeneric("organism<-", function(object, value) {
  standardGeneric("organism<-")
})

#' @name FACSData
#' @aliases FACSData
#' @export
#' @docType methods
#' @rdname FACSData
setGeneric("FACSData", function(object) {
  standardGeneric("FACSData")
})

#' @name FACSData<-
#' @aliases FACSData
#' @export
#' @docType methods
#' @rdname FACSData
setGeneric("FACSData<-", function(object, value) {
  standardGeneric("FACSData<-")
})

#' @name gene_id_type
#' @aliases gene_id_type
#' @export
#' @docType methods
#' @rdname gene_id_type
setGeneric("gene_id_type", function(object) {
  standardGeneric("gene_id_type")
})

#' @name gene_id_type<-
#' @aliases gene_id_type
#' @export
#' @docType methods
#' @rdname gene_id_type
setGeneric("gene_id_type<-", function(object, value) {
  standardGeneric("gene_id_type<-")
})

#' @name tpm
#' @aliases tpm
#' @export
#' @docType methods
#' @return a matrix of transcripts-per-million data
#' @rdname tpm
setGeneric("tpm", function(object) {standardGeneric("tpm")})

#' @name tpm<-
#' @aliases tpm
#' @export
#' @docType methods
#' @rdname tpm
setGeneric("tpm<-", function(object, value) {standardGeneric("tpm<-")})

#' @name cpm
#' @aliases cpm
#' @export
#' @docType methods
#' @return a matrix of counts-per-million values
#' @rdname cpm
setGeneric("cpm", function(object) {standardGeneric("cpm")})

#' @name cpm<-
#' @aliases cpm
#' @export
#' @docType methods
#' @rdname cpm
setGeneric("cpm<-", function(object, value) {standardGeneric("cpm<-")})

#' @name fpkm
#' @aliases fpkm
#' @export
#' @docType methods
#' @return a matrix of FPKM values
#' @rdname fpkm
setGeneric("fpkm", function(object) {standardGeneric("fpkm")})

#' @name fpkm<-
#' @aliases fpkm
#' @export
#' @docType methods
#' @rdname fpkm
setGeneric("fpkm<-", function(object, value) {standardGeneric("fpkm<-")})
