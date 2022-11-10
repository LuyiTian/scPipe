
#' @name QC_metrics
#' @aliases QC_metrics
#' @export
#' @docType methods
#' @return a dataframe of quality control matrics
#' @rdname QC_metrics
setGeneric("QC_metrics", function(object) standardGeneric("QC_metrics"))

#' @name QC_metrics
#' @aliases QC_metrics<-
#' @export
#' @docType methods
#' @rdname QC_metrics
setGeneric("QC_metrics<-", function(object, value) standardGeneric("QC_metrics<-"))


#' @name demultiplex_info
#' @aliases demultiplex_info
#' @export
#' @docType methods
#' @return a dataframe of cell barcode demultiplex information
#' @rdname demultiplex_info
setGeneric("demultiplex_info", function(object) standardGeneric("demultiplex_info"))

#' @name demultiplex_info
#' @aliases demultiplex_info<-
#' @export
#' @docType methods
#' @rdname demultiplex_info
setGeneric("demultiplex_info<-", function(object, value) standardGeneric("demultiplex_info<-"))



#' @name UMI_dup_info
#' @aliases UMI_dup_info
#' @export
#' @docType methods
#' @return a dataframe of cell UMI duplication information
#' @rdname UMI_dup_info
setGeneric("UMI_dup_info", function(object) standardGeneric("UMI_dup_info"))

#' @name UMI_dup_info
#' @aliases UMI_dup_info<-
#' @export
#' @docType methods
#' @rdname UMI_dup_info
setGeneric("UMI_dup_info<-", function(object, value) standardGeneric("UMI_dup_info<-"))

#' @name gene_id_type
#' @aliases gene_id_type
#' @export
#' @docType methods
#' @return the gene id type used by Biomart
#' @rdname gene_id_type
setGeneric("gene_id_type", function(object) standardGeneric("gene_id_type"))

#' @name gene_id_type
#' @aliases gene_id_type<-
#' @export
#' @docType methods
#' @rdname gene_id_type
setGeneric("gene_id_type<-", function(object, value) standardGeneric("gene_id_type<-"))
