# guess the organism and species from input data
.guess_attr = function(expr_mat) {
  hsp_ensembl = length(grep("^ENSG", rownames(expr_mat)))
  mm_ensembl = length(grep("^ENSMUSG", rownames(expr_mat)))
  if ((hsp_ensembl>0) & (hsp_ensembl>mm_ensembl)) {
    return(list(organism="hsapiens_gene_ensembl", gene_id_type="ensembl_gene_id"))
  }
  else if ((mm_ensembl>0) & (mm_ensembl>hsp_ensembl)) {
    return(list(organism="mmusculus_gene_ensembl", gene_id_type="ensembl_gene_id"))
  }
  else {
    return(list(organism=NA, gene_id_type=NA))
  }
}



#' Get or set quality control metrics in a SingleCellExperiment object
#' @rdname QCMetrics
#' @param object An \code{SingleCellExperiment} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return A DataFrame of quality control metrics.
#' @author Luyi Tian
#'
#' @export
#'
#' @examples
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' QualityControlInfo = new("AnnotatedDataFrame", data = as.data.frame(sc_sample_qc))
#' scd = newSCData(countData = as.matrix(sc_sample_data),
#'                QualityControlInfo = QualityControlInfo,
#'                useForExprs = "counts",
#'                organism = "mmusculus_gene_ensembl",
#'                gene_id_type = "external_gene_name")
#' QCMetrics(scd)
#'
QCMetrics.sce <- function(object) {
  if(!("scPipe" %in% names(object@int_metadata))){
    stop("`scPipe` not in `int_metadata`. cannot identify quality control columns")
  }else if(!("QC_cols" %in% names(object@int_metadata$scPipe))){
    stop("The int_metadata$scPipe does not have `QC_cols`. 
      cannot identifyquality control  columns")
  }
  return(object@int_colData[, object@int_metadata$scPipe$QC_cols])
}

#' @rdname QCMetrics
#' @aliases QCMetrics
#' @export
#'
setMethod("QCMetrics", signature(object = "SingleCellExperiment"),
          QCMetrics.sce)

#' @rdname QCMetrics
#' @aliases QCMetrics
#' @export
setReplaceMethod("QCMetrics",
                 signature="SingleCellExperiment",
                 function(object, value) {
                   if(!("scPipe" %in% names(object@int_metadata))){
                     object@int_metadata[["scPipe"]] = list(QC_cols=colnames(value))
                   }else{
                     object@int_metadata$scPipe$QC_cols = colnames(value)
                   }
                   object@int_colData = DataFrame(value)
                   #validObject(object) # could add other checks
                   return(object)
                 })


#' Get or set cell barcode demultiplx results in a SingleCellExperiment object
#' @rdname demultiplx_info
#' @param object An \code{SingleCellExperiment} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return A DataFrame of cell barcode demultiplx results.
#' @author Luyi Tian
#'
#' @export
#'
#' @examples
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' QualityControlInfo = new("AnnotatedDataFrame", data = as.data.frame(sc_sample_qc))
#' scd = newSCData(countData = as.matrix(sc_sample_data),
#'                QualityControlInfo = QualityControlInfo,
#'                useForExprs = "counts",
#'                organism = "mmusculus_gene_ensembl",
#'                gene_id_type = "external_gene_name")
#' QCMetrics(scd)
#'
demultiplx_info.sce <- function(object) {
  if(!("scPipe" %in% names(object@int_metadata))){
    stop("`scPipe` not in `int_metadata`. cannot find columns for cell barcode demultiplx results")
  }else if(!("demultiplx_info" %in% names(object@int_metadata$scPipe))){
    stop("The int_metadata$scPipe does not have `demultiplx_info`.")
  }
  return(object@int_metadata$scPipe$demultiplx_info)
  }

#' @rdname demultiplx_info
#' @aliases demultiplx_info
#' @export
#'
setMethod("demultiplx_info", signature(object = "SingleCellExperiment"),
          demultiplx_info.sce)

#' @rdname demultiplx_info
#' @aliases demultiplx_info
#' @export
setReplaceMethod("demultiplx_info",
                 signature="SingleCellExperiment",
                 function(object, value) {
                   if(!("scPipe" %in% names(object@int_metadata))){
                     object@int_metadata[["scPipe"]] = list(demultiplx_info=value)
                   }else{
                     object@int_metadata$scPipe$demultiplx_info = value
                   }
                   #validObject(object) # could add other checks
                   return(object)
                 })




#' Get or set UMI duplication results in a SingleCellExperiment object
#' @rdname UMI_dup_info
#' @param object An \code{SingleCellExperiment} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return A DataFrame of UMI duplication results.
#' @author Luyi Tian
#'
#' @export
#'
#' @examples
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' QualityControlInfo = new("AnnotatedDataFrame", data = as.data.frame(sc_sample_qc))
#' scd = newSCData(countData = as.matrix(sc_sample_data),
#'                QualityControlInfo = QualityControlInfo,
#'                useForExprs = "counts",
#'                organism = "mmusculus_gene_ensembl",
#'                gene_id_type = "external_gene_name")
#' QCMetrics(scd)
#'
UMI_dup_info.sce <- function(object) {
  if(!("scPipe" %in% names(object@int_metadata))){
    stop("`scPipe` not in `int_metadata`. cannot find columns for cell barcode demultiplx results")
  }else if(!("UMI_dup_info" %in% names(object@int_metadata$scPipe))){
    stop("The int_metadata$scPipe does not have `UMI_dup_info`.")
  }
  return(object@int_metadata$scPipe$UMI_dup_info)
}

#' @rdname UMI_dup_info
#' @aliases UMI_dup_info
#' @export
#'
setMethod("UMI_dup_info", signature(object = "SingleCellExperiment"),
          UMI_dup_info.sce)

#' @rdname UMI_dup_info
#' @aliases UMI_dup_info
#' @export
setReplaceMethod("UMI_dup_info",
                 signature="SingleCellExperiment",
                 function(object, value) {
                   if(!("scPipe" %in% names(object@int_metadata))){
                     object@int_metadata[["scPipe"]] = list(UMI_dup_info=value)
                   }else{
                     object@int_metadata$scPipe$UMI_dup_info = value
                   }
                   #validObject(object) # could add other checks
                   return(object)
                 })



#' Get or set \code{organism} from an SCData object
#' @rdname organism
#' @param object A \code{SingleCellExperiment} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return organism string
#' @author Luyi Tian
#' @export
#' @examples
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' QualityControlInfo = new("AnnotatedDataFrame", data = as.data.frame(sc_sample_qc))
#' scd = newSCData(countData = as.matrix(sc_sample_data),
#'                QualityControlInfo = QualityControlInfo,
#'                useForExprs = "counts",
#'                organism = "mmusculus_gene_ensembl",
#'                gene_id_type = "external_gene_name")
#' organism(scd)
#'
organism.sce <- function(object) {
  return(object@int_metadata$Biomart$organism)
}


#' @aliases organism
#' @rdname organism
#' @export
setMethod("organism", signature(object="SingleCellExperiment"),
          organism.sce)


#' @aliases organism
#' @rdname organism
#' @export
#' @export
setReplaceMethod("organism",
           signature="SingleCellExperiment",
           function(object, value) {
               object@int_metadata$Biomart$organism = value
               validObject(object) # could add other checks
               return(object)
             })



#' Get or set \code{gene_id_type} from an SCData object
#' @rdname gene_id_type
#' @param object An \code{\link{SCData}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return gene id type string
#' @author Luyi Tian
#'
#' @export
#'
#' @examples
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' QualityControlInfo = new("AnnotatedDataFrame", data = as.data.frame(sc_sample_qc))
#' scd = newSCData(countData = as.matrix(sc_sample_data),
#'                QualityControlInfo = QualityControlInfo,
#'                useForExprs = "counts",
#'                organism = "mmusculus_gene_ensembl",
#'                gene_id_type = "external_gene_name")
#' gene_id_type(scd)
#'
gene_id_type.sce <- function(object) {
  return(object@int_metadata$Biomart$gene_id_type)
}


#' @rdname gene_id_type
#' @aliases gene_id_type
#' @export
setMethod("gene_id_type", signature(object = "SingleCellExperiment"),
          gene_id_type.sce)


#' @aliases gene_id_type
#' @rdname gene_id_type
#' @export
setReplaceMethod("gene_id_type",
                 signature="SingleCellExperiment",
                 function(object, value) {
                   object@int_metadata$Biomart$gene_id_type = value
                   validObject(object) # could add other checks
                   return(object)
                 })


