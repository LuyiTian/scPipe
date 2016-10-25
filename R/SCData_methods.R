


#' Create a new SCData object.
#'
#' SCData extends Bioconductor's ExpressionSet class, and the same basic interface is
#' supported. \code{newSCData()} expects a matrix of expression values as its first
#' argument, with rows as features (usually genes) and columns as cells.
#' Per-feature and per-cell metadata can be supplied with the featureData and
#' phenoData arguments, respectively. Use of these optional arguments is
#' strongly encouraged. The SCData also includes a slot 'counts' to store an
#' object containing raw count data.
#'
#' @return a new SCData object
#' @author Luyi Tian
#'
#' @details for more details see \code{SCData}
#' @importFrom Biobase annotatedDataFrameFrom
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase assayDataNew
#' @importFrom Biobase assayDataElement
#' @import methods
#' @export
#' @examples TODO
#'
newSCData <- function(exprsData = NULL,
                      countData = NULL,
                      tpmData = NULL,
                      fpkmData = NULL,
                      cpmData = NULL,
                      phenoData = NULL,
                      FACSData = NULL,
                      featureData = NULL,
                      experimentData = NULL,
                      logged = FALSE,
                      gene_id_type = NULL,
                      organism = NULL,
                      logExprsOffset = 1,
                      reducedExprDimension = NULL,
                      reducedFACSDimension = NULL,
                      onesense = NULL,
                      QualityControlInfo = NULL,
                      useForExprs = c("exprs","tpm","cpm","counts","fpkm")){

  # Check that at least we have the expression data we wanted
  if (missing(useForExprs)){
    stop("Need to specify which expresssion matrix to use by `useForExprs`.")
  }
  exprs_mat <- switch(useForExprs,
                      exprs = exprsData,
                      tpm = tpmData,
                      cpm = cpmData,
                      fpkm = fpkmData,
                      counts = countData)
  if (is.null(exprs_mat)){
    stop("the data related to `useForExprs` is null. \n(i.e if you choose `count` in `useForExprs` then you should set countData)")
  }

  # Generate valid phenoData, FACSData, featureData and QualityControlInfo if not provided
  if (is.null(phenoData)){
    phenoData <- annotatedDataFrameFrom(exprs_mat, byrow = FALSE)
  }
  if (is.null(FACSData)){
    FACSData <- annotatedDataFrameFrom(exprs_mat, byrow = FALSE)
  }
  if (is.null(featureData)){
    featureData <- annotatedDataFrameFrom(exprs_mat, byrow = TRUE)
  }
  if (is.null(QualityControlInfo)){
    QualityControlInfo <- annotatedDataFrameFrom(exprs_mat, byrow = FALSE)
  }

  # set ERCC spikein information
  spikein = grepl("^ERCC", rownames(exprs_mat))
  if(any(spikein)){
    featureData$isSpike = spikein
  }
  else{
    featureData$isSpike = spikein
    print("cannot detect ERCC Spikeins from data.")
  }

  # check organism names or gene_id_type is set correctly
  if(missing(organism)){
    stop("organism cannot be NULL. \n List of possible names can be \nretrieved using the function `listDatasets`from `biomaRt` package. \n(i.e `mmusculus_gene_ensembl` or `hsapiens_gene_ensembl`)")
  }
  if(missing(gene_id_type)){
    stop("gene_id_type cannot be NULL. \n A possible list of ids can be retrieved using the function `listAttributes` from `biomaRt` package. \nthe commonly used id types are `external_gene_name`, `ensembl_gene_id` or `entrezgene`.")
  }


  # Generate valid reducedExprDimension, onesense and reducedFACSDimension if not provided
  if (is.null(reducedExprDimension)){
    reducedExprDimension = matrix(0, nrow = ncol(exprs_mat), ncol = 0)
    rownames(reducedExprDimension) = colnames(exprs_mat)
  }
  if (is.null(reducedFACSDimension)){
    reducedFACSDimension = matrix(0, nrow = ncol(exprs_mat), ncol = 0)
    rownames(reducedFACSDimension) = colnames(exprs_mat)
  }
  if (is.null(onesense)){
    onesense = matrix(0, nrow = ncol(exprs_mat), ncol = 0)
    rownames(onesense) = colnames(exprs_mat)
  }

  # Check experimentData
  expData_null <- new("MIAME",
                      name = "<your name here>",
                      lab = "<your lab here>",
                      contact = "<email address>",
                      title = "<title for this dataset>",
                      abstract = "<abstract for this dataset>",
                      url = "<your website here>",
                      other = list(
                        notes = "This dataset created from ...",
                        coauthors = c("")
                      ))
  if ( !is.null( experimentData ) ) {
    if ( is(experimentData, "MIAME") )
      expData <- experimentData
    else {
      expData <- expData_null
      warning("experimentData supplied is not an 'MIAME' object. Thus, experimentData is being set to an empty MIAME object.\n Please supply a valid 'MIAME' class object containing experiment data to experimentData(object).")
    }
  } else {
    expData <- expData_null
  }

  # Generate new SCData object
  assaydata <- assayDataNew("lockedEnvironment", exprs = exprs_mat)
  scd <- new( "SCData",
                 assayData = assaydata,
                 phenoData = phenoData,
                 featureData = featureData,
                 FACSData = FACSData,
                 experimentData = expData,
                 logExprsOffset = logExprsOffset,
                 logged = logged,
                 gene_id_type = gene_id_type,
                 organism = organism,
                 reducedExprDimension = reducedExprDimension,
                 reducedFACSDimension = reducedFACSDimension,
                 onesense = onesense,
                 QualityControlInfo = QualityControlInfo,
                 useForExprs = useForExprs)

  # Add non-null slots to assayData for SCData object, omitting null slots
  if ( !is.null(tpmData) )
    tpm(scd) <- tpmData
  if ( !is.null(fpkmData) )
    fpkm(scd) <- fpkmData
  if ( !is.null(countData) )
    counts(scd) <- countData
  if ( !is.null(cpmData) )
    cpm(scd) <- cpmData


  # Check validity of object
  validObject(scd)
  return(scd)
}


# check validity for SCData class object

setValidity("SCData", function(object) {
  msg <- NULL
  valid <- TRUE

  # Check that the dimensions of the reducedExprDimension and
  # reducedFACSDimension slot are sensible
  if ( (nrow(object@reducedExprDimension) != 0) &&
       (nrow(object@reducedExprDimension) != ncol(object)) ) {
    valid <- FALSE
    msg <- c(msg, "Number of cells in reducedExprDimension doesn't match number of cells in SCData.")
  }
  if ( (nrow(object@reducedFACSDimension) != 0) &&
       (nrow(object@reducedFACSDimension) != ncol(object)) ) {
    valid <- FALSE
    msg <- c(msg, "Number of cells in reducedFACSDimension doesn't match number of cells in SCData.")
  }

  if (valid) TRUE else msg
})



#' Get or set quality control metrics from an SCData object
#' @name QC_metrics
#' @rdname QC_metrics
#' @param object An \code{\link{SCData}} object.
#'
#' @return A \code{AnnotatedDataFrame} of quality control metrics.
#' @author Luyi Tian
#'
#' @export
#' @examples
#' TODO
#'
#'
QC_metrics.SCData <- function(object) {
  return(object@QualityControlInfo)
}

#' @name QC_metrics
#' @rdname QC_metrics
#' @aliases QC_metrics
#' @export
#'
setMethod("QC_metrics", signature(object = "SCData"),
          QC_metrics.SCData)

#' @name QC_metrics<-
#' @aliases QC_metrics
#' @rdname QC_metrics
#' @exportMethod "QC_metrics<-"
setReplaceMethod("QC_metrics",
                 signature="SCData",
                 function(object, value) {
                   object@QualityControlInfo = new("AnnotatedDataFrame", data = as.data.frame(value))
                   validObject(object) # could add other checks
                   return(object)
                 })



#' Get or set \code{reducedExprDimension} from an SCData object
#' @name DimRd_expr
#' @rdname DimRd_expr
#' @param object An \code{\link{SCData}} object.
#'
#' @return A matrix of reduced-dimension coordinates for gene
#' expression of single cell.
#' @author Luyi Tian
#'
#' @export
#' @examples
#' TODO
#'
DimRd_expr.SCData <- function(object) {
  return(object@reducedExprDimension)
}

#' @name DimRd_expr
#' @rdname DimRd_expr
#' @aliases DimRd_expr
#' @export
#'
setMethod("DimRd_expr", signature(object = "SCData"),
          DimRd_expr.SCData)


#' @name DimRd_expr<-
#' @aliases DimRd_expr
#' @rdname DimRd_expr
#' @exportMethod "DimRd_expr<-"
setReplaceMethod("DimRd_expr",
                 signature="SCData",
                 function(object, value) {
                   object@reducedExprDimension = value
                   validObject(object) # could add other checks
                   return(object)
                 })



#' Get or set \code{organism} from an SCData object
#' @name organism
#' @rdname organism
#' @param object An \code{\link{SCData}} object.
#'
#' @return organism string
#' @author Luyi Tian
#'
#' @export
#' @examples
#' TODO
#'
organism.SCData <- function(object) {
  return(object@organism)
}


#' @rdname organism
#' @name organism
#' @aliases organism
#' @export
setMethod("organism", signature(object = "SCData"),
          organism.SCData)


#' @name organism<-
#' @aliases organism
#' @rdname organism
#' @exportMethod "organism<-"
setReplaceMethod("organism",
                 signature="SCData",
                 function(object, value) {
                   object@organism = value
                   validObject(object) # could add other checks
                   return(object)
                 })



#' Get or set \code{gene_id_type} from an SCData object
#' @name gene_id_type
#' @rdname gene_id_type
#' @param object An \code{\link{SCData}} object.
#'
#' @return gene id type string
#' @author Luyi Tian
#'
#' @export
#' @examples
#' TODO
#'
gene_id_type.SCData <- function(object) {
  return(object@gene_id_type)
}


#' @name gene_id_type
#' @rdname gene_id_type
#' @aliases gene_id_type
#' @export
setMethod("gene_id_type", signature(object = "SCData"),
          gene_id_type.SCData)


#' @name gene_id_type<-
#' @rdname gene_id_type
#' @aliases gene_id_type
#' @rdname gene_id_type
#' @exportMethod "gene_id_type<-"
setReplaceMethod("gene_id_type",
                 signature="SCData",
                 function(object, value) {
                   object@gene_id_type = value
                   validObject(object) # could add other checks
                   return(object)
                 })





#' Accessors for the 'counts' element of an SCData object.
#'
#' The counts element holds the count data as a matrix of non-negative integer
#' count values, one row for each feature (gene, exon, region, etc), and one
#' column for each cell. It is an element of the assayData slot of the SCData
#' object.
#'
#' @usage
#' \S4method{counts}{SCData}(object)
#'
#' \S4method{counts}{SCData,matrix}(object)<-value
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @importFrom BiocGenerics counts
#' @aliases counts counts,SCData-method counts<-,SCData,matrix-method
#'
#' @param object a \code{SCData} object.
#' @param value an integer matrix
#' @author Luyi Tian
#' @export
#' @examples
#'
counts.SCData <- function(object) {
  object@assayData$counts
}

#' @rdname counts
#' @export
setMethod("counts", signature(object = "SCData"), counts.SCData)

#' @name counts
#' @rdname counts
#' @importFrom BiocGenerics counts<-
#' @exportMethod "counts<-"
setReplaceMethod("counts", signature(object = "SCData", value = "matrix"),
                 function(object, value) {
                   Biobase::assayDataElement(object, "counts") <- value
                   validObject(object)
                   object
                 })



#' Accessors for the 'tpm' (transcripts per million) element of an SCData object.
#'
#' The \code{tpm} element of the arrayData slot in an SCData object holds
#' a matrix containing transcripts-per-million values. It has the same dimensions
#' as the 'exprs' and 'counts' elements, which hold the transformed expression
#' data and count data, respectively.
#'
#' @usage
#' \S4method{tpm}{SCData}(object)
#'
#' \S4method{tpm}{SCData,matrix}(object)<-value
#'
#' @docType methods
#' @name tpm
#' @rdname tpm
#' @aliases tpm tpm,SCData-method tpm<-,SCData,matrix-method
#'
#' @param object a \code{SCData} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Luyi Tian
#' @export
#' @aliases tpm tpm,SCData-method tpm<-,SCData,matrix-method
#'
#' @examples
#'
tpm.SCData <- function(object) {
  object@assayData$tpm
}

#' @name tpm
#' @rdname tpm
#' @export
#' @aliases tpm,SCData-method
setMethod("tpm", signature(object = "SCData"), tpm.SCData)

#' @name tpm<-
#' @rdname tpm
#' @exportMethod "tpm<-"
#' @aliases tpm<-,SCData,matrix-method
setReplaceMethod("tpm", signature(object = "SCData", value = "matrix"),
                 function(object, value) {
                   Biobase::assayDataElement(object, "tpm") <- value
                   validObject(object)
                   object
                 })



#' Accessors for the 'cpm' (counts per million) element of an SCData object.
#'
#' The \code{cpm} element of the arrayData slot in an SCData object holds
#' a matrix containing counts-per-million values. It has the same dimensions
#' as the 'exprs' and 'counts' elements, which hold the transformed expression
#' data and count data, respectively.
#'
#' @usage
#' \S4method{cpm}{SCData}(object)
#'
#' \S4method{cpm}{SCData,matrix}(object)<-value
#'
#' @docType methods
#' @name cpm
#' @rdname cpm
#' @aliases cpm cpm,SCData-method cpm<-,SCData,matrix-method
#'
#' @param object a \code{SCData} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Luyi Tian
#' @export
#' @aliases cpm cpm,SCData-method cpm<-,SCData,matrix-method
#'
#' @examples
#'
cpmSCData <- function(object) {
  object@assayData$cpm
}

#' @name cpm
#' @rdname cpm
#' @export
#' @aliases cpm,SCData-method
setMethod("cpm", signature(object = "SCData"), cpmSCData)

#' @name cpm<-
#' @rdname cpm
#' @exportMethod "cpm<-"
#' @aliases cpm<-,SCData,matrix-method
setReplaceMethod("cpm", signature(object = "SCData", value = "matrix"),
                 function(object, value) {
                   Biobase::assayDataElement(object, "cpm") <- value
                   validObject(object)
                   object
                 })



#' Accessors for the 'fpkm' (fragments per kilobase of exon per million reads mapped) element of an SCData object.
#'
#' The \code{fpkm} element of the arrayData slot in an SCData object holds
#' a matrix containing fragments per kilobase of exon per million reads mapped
#' (FPKM) values. It has the same dimensions as the 'exprs' and 'counts'
#' elements, which hold the transformed expression data and count data,
#' respectively.
#'
#' @usage
#' \S4method{fpkm}{SCData}(object)
#'
#' \S4method{fpkm}{SCData,matrix}(object)<-value
#'
#' @docType methods
#' @name fpkm
#' @rdname fpkm
#' @aliases fpkm fpkm,SCData-method fpkm<-,SCData,matrix-method
#'
#' @param object a \code{SCData} object.
#' @param value a matrix of class \code{"numeric"}
#'
#' @author Luyi Tian
#' @export
#'
#' @examples
#'
fpkm.SCData <- function(object) {
  object@assayData$fpkm
}

#' @name fpkm
#' @rdname fpkm
#' @export
#' @aliases fpkm,SCData-method
setMethod("fpkm", signature(object = "SCData"), fpkm.SCData)

#' @name fpkm<-
#' @rdname fpkm
#' @exportMethod "fpkm<-"
#' @aliases fpkm<-,SCData,matrix-method
setReplaceMethod("fpkm", signature(object = "SCData", value = "matrix"),
                 function(object, value) {
                   Biobase::assayDataElement(object, "fpkm") <- value
                   validObject(object)
                   object
                 })
