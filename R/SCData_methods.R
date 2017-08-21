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



#' Get or set quality control metrics from an SingleCellExperiment object
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
QCMetrics.SCData <- function(object) {
  return(object@QualityControlInfo)
}

#' @rdname QCMetrics
#' @aliases QCMetrics
#' @export
#'
setMethod("QCMetrics", signature(object = "SCData"),
          QCMetrics.SCData)

#' @rdname QCMetrics
#' @aliases QCMetrics
#' @export
setReplaceMethod("QCMetrics",
                 signature="SCData",
                 function(object, value) {
                   object@QualityControlInfo = new("AnnotatedDataFrame", data = as.data.frame(value))
                   validObject(object) # could add other checks
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
                   return(object)
                 })


