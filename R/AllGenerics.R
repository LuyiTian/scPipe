#' @name DimReduceExpr
#' @export
#' @docType methods
#' @rdname DimReduceExpr
setGeneric("DimReduceExpr", function(object) {
  standardGeneric("DimReduceExpr")
})

#' @name DimReduceExpr<-
#' @export
#' @docType methods
#' @rdname DimReduceExpr
setGeneric("DimReduceExpr<-", function(object, value) {
  standardGeneric("DimReduceExpr<-")
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
#' @export
#' @docType methods
#' @rdname organism
setGeneric("organism", function(object) {
  standardGeneric("organism")
})

#' @name organism<-
#' @export
#' @docType methods
#' @rdname organism
setGeneric("organism<-", function(object, value) {
  standardGeneric("organism<-")
})

#' @name FACSData
#' @export
#' @docType methods
#' @rdname FACSData
setGeneric("FACSData", function(object) {
  standardGeneric("FACSData")
})

#' @name FACSData<-
#' @export
#' @docType methods
#' @rdname FACSData
setGeneric("FACSData<-", function(object, value) {
  standardGeneric("FACSData<-")
})

#' @name gene_id_type
#' @export
#' @docType methods
#' @rdname gene_id_type
setGeneric("gene_id_type", function(object) {
  standardGeneric("gene_id_type")
})

#' @name geneIdTdye<-
#' @export
#' @docType methods
#' @rdname gene_id_type
setGeneric("gene_id_type<-", function(object, value) {
  standardGeneric("gene_id_type<-")
})

#' @name tpm
#' @export
#' @docType methods
#' @return a matrix of transcripts-per-million data
#' @rdname tpm
setGeneric("tpm", function(object) {standardGeneric("tpm")})

#' @name tpm<-
#' @export
#' @docType methods
#' @rdname tpm
setGeneric("tpm<-", function(object, value) {standardGeneric("tpm<-")})

#' @name cpm
#' @export
#' @docType methods
#' @return a matrix of counts-per-million values
#' @rdname cpm
setGeneric("cpm", function(object) {standardGeneric("cpm")})

#' @name cpm<-
#' @export
#' @docType methods
#' @rdname cpm
setGeneric("cpm<-", function(object, value) {standardGeneric("cpm<-")})

#' @name fpkm
#' @export
#' @docType methods
#' @return a matrix of FPKM values
#' @rdname fpkm
setGeneric("fpkm", function(object) {standardGeneric("fpkm")})

#' @name fpkm<-
#' @export
#' @docType methods
#' @rdname fpkm
setGeneric("fpkm<-", function(object, value) {standardGeneric("fpkm<-")})

