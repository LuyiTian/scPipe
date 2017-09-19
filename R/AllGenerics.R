
#' @name accessors
#' @aliases QC_metrics
#' @export
#' @docType methods
#' @return a dataframe of quality control matrics
#' @rdname accessors
setGeneric("QC_metrics", function(object) standardGeneric("QC_metrics"))

#' @name accessors
#' @aliases QC_metrics<-
#' @export
#' @docType methods
#' @rdname accessors
setGeneric("QC_metrics<-", function(object, value) standardGeneric("QC_metrics<-"))


#' @name accessors
#' @aliases demultiplex_info
#' @export
#' @docType methods
#' @return a dataframe of cell barcode demultiplex information
#' @rdname accessors
setGeneric("demultiplex_info", function(object) standardGeneric("demultiplex_info"))

#' @name accessors
#' @aliases demultiplex_info<-
#' @export
#' @docType methods
#' @rdname accessors
setGeneric("demultiplex_info<-", function(object, value) standardGeneric("demultiplex_info<-"))



#' @name accessors
#' @aliases UMI_dup_info
#' @export
#' @docType methods
#' @return a dataframe of cell UMI duplication information
#' @rdname accessors
setGeneric("UMI_dup_info", function(object) standardGeneric("UMI_dup_info"))

#' @name accessors
#' @aliases UMI_dup_info<-
#' @export
#' @docType methods
#' @rdname accessors
setGeneric("UMI_dup_info<-", function(object, value) standardGeneric("UMI_dup_info<-"))
