# guess the organism and species from input data
.guess_attr = function(row_names) {
  hsp_ensembl = length(grep("^ENSG",row_names))
  mm_ensembl = length(grep("^ENSMUSG", row_names))
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



validObject = function(object){
  if (!is(object, "SingleCellExperiment")) {
    stop("object must be an `SingleCellExperiment` object.")
  }
  
  if(!("scPipe" %in% names(object@int_metadata))){
    object@int_metadata$scPipe$version = packageVersion("scPipe")  # set version information
  }
  if(!("QC_cols" %in% names(object@int_metadata$scPipe))){
    QC_metrics(object) = DataFrame(row.names = colnames(object)) # create a empty QC metrics if not exists
  }
  
  if(min(dim(object)) == 0){
    stop("The dimension of sce should be larger than zero.")
  }else if(is.null(rownames(object)) | is.null(colnames(object))){
    stop("rowname/colname does not exists for sce.")
  }else if(!all(rownames(QC_metrics(object)) == colnames(object))){
    stop("The rownames of QC metrics is not consistant with column names of the object.")
  }
  
  if(is.null(gene_id_type(object))){ # check the gene id type
    gene_id_type(object) = NA
  }else if (!is.na(gene_id_type(object))){
    if(gene_id_type(object) == "NA"){
        gene_id_type(object) = NA
    }
  }
  
  if(is.null(organism(object))){ # check the organism
    organism(object) = NA
  }else if (!is.na(organism(object))){
    if(organism(object) == "NA"){
      organism(object) = NA
    }
  }
  
  if(is.na(organism(object)) | is.na(gene_id_type(object))){
    tmp_res = .guess_attr(rownames(object))
    if((!is.na(tmp_res$organism)) & (!is.na(tmp_res$gene_id_type))){
      gene_id_type(object) = tmp_res$gene_id_type
      organism(object) = tmp_res$organism
      print(paste("organism/gene_id_type not provided. make a guess:", 
                  tmp_res$organism,
                  "/",
                  tmp_res$gene_id_type))
    } 
  }
  return(object)
}


#' Get or set quality control metrics in a SingleCellExperiment object
#' @rdname QC_metrics
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
QC_metrics.sce <- function(object) {
  if(!("scPipe" %in% names(object@int_metadata))){
    warning("`scPipe` not in `int_metadata`. cannot identify quality control columns")
    return(NULL)
  }else if(!("QC_cols" %in% names(object@int_metadata$scPipe))){
    warning("The int_metadata$scPipe does not have `QC_cols`. 
      cannot identifyquality control  columns")
    return(NULL)
  }
  return(object@int_colData[, object@int_metadata$scPipe$QC_cols])
}

#' @rdname QC_metrics
#' @aliases QC_metrics
#' @export
#'
setMethod("QC_metrics", signature(object = "SingleCellExperiment"),
          QC_metrics.sce)

#' @rdname QC_metrics
#' @aliases QC_metrics
#' @export
setReplaceMethod("QC_metrics",
                 signature="SingleCellExperiment",
                 function(object, value) {
                   if(!("scPipe" %in% names(object@int_metadata))){
                     object@int_metadata[["scPipe"]] = list(QC_cols=colnames(value))
                   }else{
                     object@int_metadata$scPipe$QC_cols = colnames(value)
                   }
                   object@int_colData = DataFrame(value)
                   object = validObject(object) # could add other checks
                   return(object)
                 })


#' @title demultiplex_info
#' 
#' @description Get or set cell barcode demultiplx results in a SingleCellExperiment object
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
demultiplex_info.sce <- function(object) {
  if(!("scPipe" %in% names(object@int_metadata))){
    warning("`scPipe` not in `int_metadata`. cannot find columns for cell barcode demultiplex results")
    return(NULL)
  }else if(!("demultiplex_info" %in% names(object@int_metadata$scPipe))){
    warning("The int_metadata$scPipe does not have `demultiplex_info`.")
    return(NULL)
  }
  return(object@int_metadata$scPipe$demultiplex_info)
  }

#' @rdname demultiplex_info
#' @aliases demultiplex_info
#' @export
#'
setMethod("demultiplex_info", signature(object = "SingleCellExperiment"),
          demultiplex_info.sce)

#' @rdname demultiplex_info
#' @aliases demultiplex_info
#' @export
setReplaceMethod("demultiplex_info",
                 signature="SingleCellExperiment",
                 function(object, value) {
                   if(!("scPipe" %in% names(object@int_metadata))){
                     object@int_metadata[["scPipe"]] = list(demultiplex_info=value)
                   }else{
                     object@int_metadata$scPipe$demultiplex_info = value
                   }
                   object = validObject(object) # could add other checks
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
    warning("`scPipe` not in `int_metadata`. cannot find columns for cell barcode demultiplex results")
    return(NULL)
  }else if(!("UMI_dup_info" %in% names(object@int_metadata$scPipe))){
    warning("The int_metadata$scPipe does not have `UMI_dup_info`.")
    return(NULL)
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
                   object = validObject(object) # could add other checks
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
               object = validObject(object) # could add other checks
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
                   object = validObject(object) # could add other checks
                   return(object)
                 })


