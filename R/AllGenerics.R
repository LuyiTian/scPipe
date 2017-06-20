

#' @name DimRd_expr
#' @export
#' @docType methods
#' @rdname DimRd_expr
setGeneric("DimRd_expr", function(object) {
  standardGeneric("DimRd_expr")
})

#' @name DimRd_expr<-
#' @export
#' @docType methods
#' @rdname DimRd_expr
setGeneric("DimRd_expr<-", function(object, value) {
  standardGeneric("DimRd_expr<-")
})

#' @name QC_metrics
#' @export
#' @docType methods
#' @rdname QC_metrics
setGeneric("QC_metrics", function(object) {
  standardGeneric("QC_metrics")
})

#' @name QC_metrics<-
#' @export
#' @docType methods
#' @rdname QC_metrics
setGeneric("QC_metrics<-", function(object, value) {
  standardGeneric("QC_metrics<-")
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

#' @name gene_id_type<-
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

