# guess the organism and species from input data
.guess_attr <- function(row_names) {
    hsp_ensembl <- length(grep("^ENSG",row_names))
    mm_ensembl <- length(grep("^ENSMUSG", row_names))
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


# check the object, fix empty slot with default.
validObject <- function(object){
    if (!is(object, "SingleCellExperiment")) {
        stop("object must be an `SingleCellExperiment` object.")
    }

    if(!("scPipe" %in% names(object@metadata))){
        object@metadata$scPipe$version <- packageVersion("scPipe")  # set version information
    }

    if(min(dim(object)) == 0){
        stop("The dimension of sce should be larger than zero.")
    }else if(is.null(rownames(object)) | is.null(colnames(object))){
        stop("rowname/colname does not exists for sce.")
    }else if(!all(rownames(QC_metrics(object)) == colnames(object))){
        stop("The rownames of QC metrics is not consistent with column names of the object.")
    }


    if(any(is.null(organism(object)) || is.na(organism(object)))){
        tmp_res <- .guess_attr(rownames(object))
        if((!is.na(tmp_res$organism)) & (!is.na(tmp_res$gene_id_type))){
            gene_id_type(object) <- tmp_res$gene_id_type
            organism(object) <- tmp_res$organism
            message("organism/gene_id_type not provided. Make a guess:",
                    tmp_res$organism,
                    "/",
                    tmp_res$gene_id_type)
        } else {
            gene_id_type(object) <- "NA"
            organism(object) <- "NA"
        }
    }
    return(object)
}


#' Get or set quality control metrics in a SingleCellExperiment object
#' @rdname QC_metrics
#' @param object A \code{\link{SingleCellExperiment}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return A DataFrame of quality control metrics.
#' @author Luyi Tian
#'
#' @importFrom S4Vectors DataFrame SimpleList
#'
#' @export
#'
#' @examples
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' QC_metrics(sce) = sc_sample_qc
#'
#' head(QC_metrics(sce))
#'
QC_metrics.sce <- function(object) {
    if(!("scPipe" %in% names(object@metadata))){
        warning("`scPipe` not in `metadata`. Cannot identify quality control columns")
        return(NULL)
    } else if (!("QC_cols" %in% names(object@metadata$scPipe))) {
        warning("metadata$scPipe is missing `QC_cols`. Cannot identify quality control columns")
        return(NULL)
    }
    return(colData(object)[, object@metadata$scPipe$QC_cols])
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
setReplaceMethod(
    "QC_metrics",
    signature = "SingleCellExperiment",
    function(object, value) {
        if (!("scPipe" %in% names(object@metadata))) {
            object@metadata[["scPipe"]] <- list(QC_cols=colnames(value))
        } else {
            object@metadata$scPipe$QC_cols <- colnames(value)
        }

        colData(object)[, colnames(value)] <- DataFrame(value)

        return(object)
    }
)


#' @title demultiplex_info
#'
#' @description Get or set cell barcode demultiplex results in a SingleCellExperiment object
#' @rdname demultiplex_info
#' @param object A \code{\link{SingleCellExperiment}} object.
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
#' sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' demultiplex_info(sce) = cell_barcode_matching
#' UMI_dup_info(sce) = UMI_duplication
#'
#' demultiplex_info(sce)
#'
demultiplex_info.sce <- function(object) {
    if(!("scPipe" %in% names(object@metadata))){
        warning("`scPipe` not in `metadata`. Cannot find columns for cell barcode demultiplex results")
        return(NULL)
    }else if(!("demultiplex_info" %in% names(object@metadata$scPipe))){
        warning("The metadata$scPipe does not have `demultiplex_info`.")
        return(NULL)
    }
    return(object@metadata$scPipe$demultiplex_info)
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
                    if(!("scPipe" %in% names(object@metadata))){
                        object@metadata[["scPipe"]] <- list(demultiplex_info=value)
                    }else{
                        object@metadata$scPipe$demultiplex_info <- value
                    }
                    object <- validObject(object) # could add other checks
                    return(object)
})




#' Get or set UMI duplication results in a SingleCellExperiment object
#' @rdname UMI_dup_info
#' @param object A \code{\link{SingleCellExperiment}} object.
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
#' sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' demultiplex_info(sce) = cell_barcode_matching
#' UMI_dup_info(sce) = UMI_duplication
#'
#' head(UMI_dup_info(sce))
#'
UMI_dup_info.sce <- function(object) {
    if(!("scPipe" %in% names(object@metadata))){
        warning("`scPipe` not in `metadata`. Cannot find columns for cell barcode demultiplex results")
        return(NULL)
    }else if(!("UMI_dup_info" %in% names(object@metadata$scPipe))){
        warning("The metadata$scPipe does not have `UMI_dup_info`.")
        return(NULL)
    }
    return(object@metadata$scPipe$UMI_dup_info)
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
                    if(!("scPipe" %in% names(object@metadata))){
                        object@metadata[["scPipe"]] <- list(UMI_dup_info=value)
                    }else{
                        object@metadata$scPipe$UMI_dup_info <- value
                    }
                    object <- validObject(object) # could add other checks
                    return(object)
})



#' Get or set \code{organism} from a SingleCellExperiment object
#' @rdname organism
#' @param object A \code{\link{SingleCellExperiment}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @importFrom BiocGenerics organism organism<-
#' @return organism string
#' @author Luyi Tian
#' @export
#' @examples
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' demultiplex_info(sce) = cell_barcode_matching
#' UMI_dup_info(sce) = UMI_duplication
#'
#' organism(sce)
#'
organism.sce <- function(object) {
    return(object@metadata$Biomart$organism)
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
setReplaceMethod("organism",signature="SingleCellExperiment",
                function(object, value) {
                    if(is.null(value)){
                        object@metadata$Biomart$organism <- NA
                        }else if(value == "NA"){
                        object@metadata$Biomart$organism <- NA
                        }else{
                            object@metadata$Biomart$organism <- value
                            }
                    return(object)
                })



#' Get or set \code{gene_id_type} from a SingleCellExperiment object
#' @rdname gene_id_type
#' @param object A \code{\link{SingleCellExperiment}} object.
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
#' sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' demultiplex_info(sce) = cell_barcode_matching
#' UMI_dup_info(sce) = UMI_duplication
#'
#' gene_id_type(sce)
#'
gene_id_type.sce <- function(object) {
    return(object@metadata$Biomart$gene_id_type)
}


#' @rdname gene_id_type
#' @aliases gene_id_type
#' @export
setMethod("gene_id_type", signature(object = "SingleCellExperiment"),
            gene_id_type.sce)


#' @aliases gene_id_type
#' @rdname gene_id_type
#' @export
setReplaceMethod("gene_id_type",signature="SingleCellExperiment",
                function(object, value) {
                    if(is.null(value)){
                        object@metadata$Biomart$gene_id_type <- NA
                    }else if(value == "NA"){
                        object@metadata$Biomart$gene_id_type <- NA
                    }else{
                        object@metadata$Biomart$gene_id_type <- value
                    }
                    return(object)
                })

# -------------------------------------------------------- scATAC-seq -----------------------------------------------------------
#' Get or set \code{feature_info} from a SingleCellExperiment object
#' @rdname feature_info
#' @param object A \code{\link{SingleCellExperiment}} object.
#' @param value Value to be assigned to corresponding object.
#' @author Shani Amarasinghe
#' @return A DataFrame of feature information
#' @export
#'
feature_info.sce <- function(object) {
    if(!("scPipe" %in% names(object@metadata))){
        warning("`scPipe` not in `metadata`.")
        return(NULL)
    }else if(!("feature_cols" %in% names(object@metadata$scPipe))){
        warning("The metadata$scPipe does not have `feature_cols`.")
        return(NULL)
    }
    return(rowData(object)[, object@metadata$scPipe$feature_cols])
}


#' @rdname feature_info
#' @aliases feature_info
#' @export
setMethod("feature_info", signature(object = "SingleCellExperiment"),
        feature_info.sce)


#' @aliases feature_info
#' @rdname feature_info
#' @export
setReplaceMethod("feature_info",signature="SingleCellExperiment",
                function(object, value) {
                    feature <- NULL
                    value <- subset(value, select = -c(feature))
                    if (!("scPipe" %in% names(object@metadata))) {
                        object@metadata[["scPipe"]] <= list(feature_cols=colnames(value))
                    } else {
                        object@metadata$scPipe$feature_cols <- colnames(value)
                    }
                    rowData(object)[, colnames(value)] <- DataFrame(value)
                    return(object)
                })


#' Get or set \code{feature_type} from a SingleCellExperiment object
#' @rdname feature_type
#' @param object A \code{\link{SingleCellExperiment}} object.
#' @param value Value to be assigned to corresponding object.
#' @author Shani Amarasinghe 
#' @return A string representing the feature type
#'
#' @export
#'
feature_type.sce <- function(object) {
    return(object@metadata$scPipe$feature_type)
}


#' @rdname feature_type
#' @aliases feature_type
#' @export
setMethod("feature_type", signature(object = "SingleCellExperiment"),
            feature_type.sce)


#' @aliases feature_type
#' @rdname feature_type
#' @export
setReplaceMethod("feature_type",signature="SingleCellExperiment",
                function(object, value) {
                    if(is.null(value) || value == "NA"){
                        object@metadata$scPipe$feature_type <- NA
                    }else{
                        object@metadata$scPipe$feature_type <- value
                    }
                    return(object)
                })
