#' merge matrix for multiple SCData object
#'
.merge_mat <- function(scd_list, func) {
  if (all(unlist(lapply(scd_list, function(x) {!is.null(func(x))})))) {
    all_exprs = lapply(scd_list, function(x) {func(x)})
    all_gene_id = rownames(all_exprs[[1]])
    all_cell_id = colnames(all_exprs[[1]])
    for (i in 2:length(all_exprs)) {
      all_gene_id = union(all_gene_id, rownames(all_exprs[[i]]))
      all_cell_id = c(all_cell_id, colnames(all_exprs[[i]]))
    }
    merged_exprs = matrix(0,
                          nrow = length(all_gene_id),
                          ncol = length(all_cell_id),
                          dimnames = list(all_gene_id, all_cell_id))
    for (i in 1:length(all_exprs)) {
      merged_exprs[rownames(all_exprs[[i]]), colnames(all_exprs[[i]])] = all_exprs[[i]]
    }
  }
  else {
    stop("all data should contain expression matrix.")
  }
  return(merged_exprs)
}


#' guess the organism and species from input data
#' 
.guess_attr <- function(expr_mat) {
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
#' @examples
#' #TODO
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
                      useForExprs = c("exprs", "tpm", "cpm", "counts", "fpkm")) {

  # Check that at least we have the expression data we wanted
  if (missing(useForExprs)) {
    stop("Need to specify which expresssion matrix to use by `useForExprs`.")
  }
  exprs_mat <- switch(useForExprs,
                      exprs = exprsData,
                      tpm = tpmData,
                      cpm = cpmData,
                      fpkm = fpkmData,
                      counts = countData)
  if (is.null(exprs_mat)) {
    stop("the data related to `useForExprs` is null. \n(i.e if you choose `count` in `useForExprs` then you should set countData)")
  }

  # Generate valid phenoData, FACSData, featureData and QualityControlInfo if not provided
  if (is.null(phenoData)) {
    phenoData <- annotatedDataFrameFrom(exprs_mat, byrow = FALSE)
  }
  if (is.null(FACSData)) {
    FACSData <- annotatedDataFrameFrom(exprs_mat, byrow = FALSE)
  }
  if (is.null(featureData)) {
    featureData <- annotatedDataFrameFrom(exprs_mat, byrow = TRUE)
  }
  if (is.null(QualityControlInfo)) {
    QualityControlInfo <- annotatedDataFrameFrom(exprs_mat, byrow = FALSE)
  }

  # set ERCC spikein information
  spikein = grepl("^ERCC", rownames(exprs_mat))
  if (any(spikein)) {
    featureData$isSpike = spikein
  }
  else {
    featureData$isSpike = spikein
    print("cannot detect ERCC Spikeins from data.")
  }

  # Generate valid reducedExprDimension, onesense and reducedFACSDimension if not provided
  if (is.null(reducedExprDimension)) {
    reducedExprDimension = matrix(0, nrow = ncol(exprs_mat), ncol = 0)
    rownames(reducedExprDimension) = colnames(exprs_mat)
  }
  if (is.null(reducedFACSDimension)) {
    reducedFACSDimension = matrix(0, nrow = ncol(exprs_mat), ncol = 0)
    rownames(reducedFACSDimension) = colnames(exprs_mat)
  }
  if (is.null(onesense)) {
    onesense = matrix(0, nrow = ncol(exprs_mat), ncol = 0)
    rownames(onesense) = colnames(exprs_mat)
  }

  # Check experimentData
  expData_null = new("MIAME",
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
      expData = experimentData
    else {
      expData = expData_null
      warning("experimentData supplied is not an 'MIAME' object. Thus, experimentData is being set to an empty MIAME object.\n Please supply a valid 'MIAME' class object containing experiment data to experimentData(object).")
    }
  } else {
    expData = expData_null
  }

  # Generate new SCData object
  assaydata = assayDataNew("lockedEnvironment", exprs = exprs_mat)
  scd = new("SCData",
            assayData = assaydata,
            phenoData = phenoData,
            featureData = featureData,
            FACSData = FACSData,
            experimentData = expData,
            logExprsOffset = logExprsOffset,
            logged = logged,
            reducedExprDimension = reducedExprDimension,
            reducedFACSDimension = reducedFACSDimension,
            onesense = onesense,
            QualityControlInfo = QualityControlInfo,
            useForExprs = useForExprs)

  # check organism names or gene_id_type is set correctly
  tmp_res = .guess_attr(exprs_mat)
  if (is.null(organism)) {
    if (is.na(tmp_res$organism)) {
      stop("organism cannot be NULL. \n List of possible names can be \nretrieved using the function `listDatasets`from `biomaRt` package. \n(i.e `mmusculus_gene_ensembl` or `hsapiens_gene_ensembl`)")
    }
    else {
      print(paste("organism not provided. make a guess:", tmp_res$organism))
      organism(scd) = tmp_res$organism
    }
  }
  else {
    # set organism
    organism(scd) = organism
  }
  if (is.null(gene_id_type)) {
    if (is.na(tmp_res$gene_id_type)) {
      stop("gene_id_type cannot be NULL. \n A possible list of ids can be retrieved using the function `listAttributes` from `biomaRt` package. \nthe commonly used id types are `external_gene_name`, `ensembl_gene_id` or `entrezgene`.")
    }
    else {
      print(paste("gene_id_type not provided. make a guess:", tmp_res$gene_id_type))
      gene_id_type(scd) = tmp_res$gene_id_type
    }
  }
  else {
    gene_id_type(scd) = gene_id_type
  }


  # Add non-null slots to assayData for SCData object, omitting null slots
  if ( !is.null(tpmData) )
    tpm(scd) = tpmData
  if ( !is.null(fpkmData) )
    fpkm(scd) = fpkmData
  if ( !is.null(countData) )
    counts(scd) = countData
  if ( !is.null(cpmData) )
    cpm(scd) = cpmData

  # remove cells that have zero counts
  scd = scd[, colSums(exprs_mat)>0]

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
  if ( (length(object@reducedExprDimension) != 0) &&
       (nrow(object@reducedExprDimension) != ncol(object)) ) {
    valid <- FALSE
    msg <- c(msg, "Number of cells in reducedExprDimension doesn't match number of cells in SCData.")
  }
  if ( (length(object@reducedFACSDimension) != 0) &&
       (nrow(object@reducedFACSDimension) != ncol(object)) ) {
    valid <- FALSE
    msg <- c(msg, "Number of cells in reducedFACSDimension doesn't match number of cells in SCData.")
  }

  if (valid) TRUE else msg
})


#' Subsetting SCData Objects
#'
#' Subset method for SCData objects, which subsets both the expression data,
#' phenotype data, feature data and other slots in the object.
#'
#' @return an SCData object
#' @rdname SCData-subset
#' @name SCData-subset
#' @inheritParams base::Extract
#' @param i,j,... indices specifying elements to extract or replace. Indices
#' are numeric or character vectors or empty (missing) or \code{NULL}. Numeric
#' values are coerced to integer as by \code{\link[base]{as.integer}} (and hence
#' truncated towards zero). Character vectors will be matched to the names of
#' the object (or for matrices/arrays, the dimnames): see
#' \code{\link[base]{Extract}} for further details.
#'
#' For \code{[}-indexing only: \code{i, j, ...} can be logical vectors, indicating
#' elements/slices to select. Such vectors are recycled if necessary to match
#' the corresponding extent. \code{i, j, ...} can also be negative integers,
#' indicating elements/slices to leave out of the selection. When indexing
#' arrays by \code{[} a single argument i can be a matrix with as many columns
#' as there are dimensions of \code{x}; the result is then a vector with
#' elements corresponding to the sets of indices in each row of \code{i}. An
#' index value of \code{NULL} is treated as if it were \code{integer(0)}.
#'
#' @aliases [,SCData,ANY-method [,SCData,ANY,ANY-method [,SCData,ANY,ANY,ANY-method
#' @rdname SCData-subset
#' @export
#' @seealso \code{\link[base]{Extract}}
#'
setMethod('[', 'SCData', function(x, i, j, drop=FALSE) {
  if ( !missing(i) && missing(j) ) {
    # select features
    x = selectMethod('[', 'ExpressionSet')(x, i, , drop = drop)
  }
  else if ( missing(i) && !missing(j) ) {
    # select cells
    x = selectMethod('[', 'ExpressionSet')(x, , j, drop = drop)
    if ( length(x@reducedExprDimension) != 0 ) {
      x@reducedExprDimension =
        as.matrix(x@reducedExprDimension[j, , drop = drop])
    }
    if ( length(x@reducedFACSDimension) != 0 ) {
      x@reducedFACSDimension =
        as.matrix(x@reducedFACSDimension[j, , drop = drop])
    }
    if ( length(x@onesense) != 0 ) {
      x@onesense =
        as.matrix(x@onesense[j, , drop = drop])
    }
    if ( length(x@QualityControlInfo) != 0 ) {
      x@QualityControlInfo =
        x@QualityControlInfo[j, , drop = drop]
    }
    x@FACSData = x@FACSData[j, , drop = drop]
  }
  else if ( !missing(i) && !missing(j) ) {
    # selcet features (i) and cells (j)
    x <- selectMethod('[', 'ExpressionSet')(x, i, j, drop = drop)
    if ( length(x@reducedExprDimension) != 0 ) {
      x@reducedExprDimension =
        as.matrix(x@reducedExprDimension[j, , drop = drop])
    }
    if ( length(x@reducedFACSDimension) != 0 ) {
      x@reducedFACSDimension =
        as.matrix(x@reducedFACSDimension[j, , drop = drop])
    }
    if ( length(x@onesense) != 0 ) {
      x@onesense =
        as.matrix(x@onesense[j, , drop = drop])
    }
    if ( length(x@QualityControlInfo) != 0 ) {
      x@QualityControlInfo =
        x@QualityControlInfo[j, , drop = drop]
    }
    x@FACSData = x@FACSData[j, , drop = drop]
  }
  ## Check validity of object
  validObject(x)
  return(x)
})

#' merge multiple SCData object
#' @rdname mergeSCData
#' @name mergeSCData
#' @param ... multiple SCDatas. They shold have the same value for class attribute.
#' @param all only contains interset for features or union.
#' @param batch (optional) batch information
#'
#' @importFrom Biobase varLabels 
#'
#' @export
#'
mergeSCData <- function(...,
                        all = TRUE,
                        batch = NULL) {
  scd_list <- list(...)
  if (!is.null(batch)) {
    if (!(length(batch) == length(scd_list))) {
      stop("the length of batch should equal to the number of dataset")
    }
  }
  else {
    batch = 1:length(scd_list)
  }


  if (length(scd_list) < 2) {
    stop("should at least contain two SCData object.")
  }
  if (!all(unlist(lapply(scd_list, function(x) {is(x, "SCData")})))) {
    stop("all data should be SCData object")
  }
  logged = scd_list[[1]]@logged
  if (!all(unlist(lapply(scd_list, function(x) {x@logged == logged})))) {
    stop("data do not have the same value for the 'logged' slot.")
  }

  useForExprs = scd_list[[1]]@useForExprs
  if (!all(unlist(lapply(scd_list, function(x) {x@useForExprs == useForExprs})))) {
    stop("data do not have the same value for the 'useForExprs' slot.")
  }

  logExprsOffset = scd_list[[1]]@logExprsOffset
  if (!all(unlist(lapply(scd_list, function(x) {x@logExprsOffset == logExprsOffset})))) {
    stop("data do not have the same value for the 'logExprsOffset' slot.")
  }

  the_organism = organism.SCData(scd_list[[1]])
  if (!all(unlist(lapply(scd_list, function(x) {organism.SCData(x) == the_organism})))) {
    stop("data do not have the same value for the 'organism' slot.")
  }

  gene_id_type = scd_list[[1]]@gene_id_type
  if (!all(unlist(lapply(scd_list, function(x) {x@gene_id_type == gene_id_type})))) {
    stop("data do not have the same value for the 'gene_id_type' slot.")
  }

  print("merge expression matrix")
  merged_exprs = .merge_mat(scd_list, exprs)

  print("merge phenotype")
  ph_col = varLabels(scd_list[[1]])
  merged_ph = NULL
  if (length(ph_col)>0) {
    if (all(unlist(lapply(scd_list, function(x) {varLabels(x) == ph_col})))) {
      merged_ph =
        AnnotatedDataFrame(data=Reduce(rbind, lapply(scd_list, function(x) {pData(x)})))
    }
    else {
      stop("the colnames in phenoData should be the same for all data.")
    }
  }
  if ("batch" %in% colnames(merged_ph)) {
    print("already contains batch information. ignore batch argument.")
  }
  else {
    batch_num =unname(unlist(lapply(scd_list, function(x) {nrow(pData(x))})))
    merged_ph$batch = rep(batch, times=as.vector(batch_num))
  }


  print("merge quality control metrics")
  qc_col = varLabels(QCMetrics(scd_list[[1]]))
  merged_qc = NULL
  if (length(qc_col)>0) {
    if (all(unlist(lapply(scd_list, function(x) {varLabels(QCMetrics(x)) == qc_col})))) {
      merged_qc =
        AnnotatedDataFrame(data=Reduce(rbind, (lapply(scd_list, function(x) {pData(QCMetrics(x))}))))
    }
    else {
      stop("the colnames in QCMetrics should be the same for all data.")
    }
  }

  print("merge FACS data")
  fac_col = varLabels(FACSData(scd_list[[1]]))
  merged_fac = NULL
  if (length(fac_col)>0) {
    if (all(unlist(lapply(scd_list, function(x) {varLabels(FACSData(x)) == fac_col})))) {
      merged_fac =
        AnnotatedDataFrame(data=Reduce(rbind, (lapply(scd_list, function(x) {pData(FACSData(x))}))))
    }
    else {
      stop("the colnames in phenoData should be the same for all data.")
    }
  }
  print("create new SCData object")
  new_scd <- newSCData(exprsData = merged_exprs,
                       phenoData = merged_ph,
                       #featureData = NULL, #TODO
                       FACSData = merged_fac,
                       #experimentData = NULL, #TODO
                       logExprsOffset = logExprsOffset,
                       logged = logged,
                       gene_id_type = gene_id_type,
                       organism = the_organism,
                       #reducedExprDimension = NULL,
                       #reducedFACSDimension = NULL,
                       #onesense = NULL,
                       QualityControlInfo = merged_qc,
                       useForExprs = "exprs")

  if (!is.null(fpkm(scd_list[[1]]))) {
    fpkm(new_scd) = .merge_mat(scd_list, fpkm)
  }

  if (!is.null(cpm(scd_list[[1]]))) {
    cpm(new_scd) = .merge_mat(scd_list, cpm)
  }

  if (!is.null(counts(scd_list[[1]]))) {
    counts(new_scd) = .merge_mat(scd_list, counts)
  }

  if (!is.null(tpm(scd_list[[1]]))) {
    tpm(new_scd) = .merge_mat(scd_list, tpm)
  }
  return(new_scd)
}

#' Get or set quality control metrics from an SCData object
#' @name QCMetrics
#' @rdname QCMetrics
#' @param object An \code{\link{SCData}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return A \code{AnnotatedDataFrame} of quality control metrics.
#' @author Luyi Tian
#'
#' @export
#'
#'
QCMetrics.SCData <- function(object) {
  return(object@QualityControlInfo)
}

#' @name QCMetrics
#' @rdname QCMetrics
#' @aliases QCMetrics
#' @export
#'
setMethod("QCMetrics", signature(object = "SCData"),
          QCMetrics.SCData)

#' @name QCMetrics<-
#' @aliases QCMetrics
#' @rdname QCMetrics
#' @exportMethod "QCMetrics<-"
setReplaceMethod("QCMetrics",
                 signature="SCData",
                 function(object, value) {
                   object@QualityControlInfo = new("AnnotatedDataFrame", data = as.data.frame(value))
                   validObject(object) # could add other checks
                   return(object)
                 })


#' Get or set \code{organism} from an SCData object
#' @name organism
#' @rdname organism
#' @param object An \code{\link{SCData}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return organism string
#' @author Luyi Tian
#' @export
#'
organism.SCData <- function(object) {
  return(object@organism)
}


#' @name organism
#' @aliases organism
#' @rdname organism
#' @export
setMethod("organism", signature(object="SCData"),
          organism.SCData)


#' @name organism<-
#' @aliases organism
#' @rdname organism
#' @exportMethod "organism<-"
#' @export
setReplaceMethod("organism",
           signature="SCData",
           function(object, value) {
               object@organism = value
               validObject(object) # could add other checks
               return(object)
             })


#' Get or set FACs data for an SCData object
#' @name FACSData
#' @rdname FACSData
#' @param object An \code{\link{SCData}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return A \code{AnnotatedDataFrame} of FACS data.
#' @author Luyi Tian
#'
#' @export
#'
#'
FACSData.SCData <- function(object) {
  return(object@FACSData)
}

#' @name FACSData
#' @rdname FACSData
#' @aliases FACSData
#' @export
#'
setMethod("FACSData", signature(object = "SCData"),
          FACSData.SCData)

#' @name FACSData<-
#' @aliases FACSData
#' @rdname FACSData
#' @exportMethod "FACSData<-"
setReplaceMethod("FACSData",
                 signature="SCData",
                 function(object, value) {
                   object@FACSData = new("AnnotatedDataFrame", data = as.data.frame(value))
                   validObject(object) # could add other checks
                   return(object)
                 }
)

#' Get or set \code{reducedExprDimension} from an SCData object
#' @name DimReduceExpr
#' @rdname DimReduceExpr
#' @param object An \code{\link{SCData}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return A matrix of reduced-dimension coordinates for gene
#' expression of single cell.
#' @author Luyi Tian
#'
#' @export
#'
DimReduceExpr.SCData <- function(object) {
  return(object@reducedExprDimension)
}

#' @name DimReduceExpr
#' @rdname DimReduceExpr
#' @aliases DimReduceExpr
#' @export
#'
setMethod("DimReduceExpr", signature(object = "SCData"),
          DimReduceExpr.SCData)


#' @name DimReduceExpr<-
#' @aliases DimReduceExpr
#' @rdname DimReduceExpr
#' @exportMethod "DimReduceExpr<-"
setReplaceMethod("DimReduceExpr",
                 signature="SCData",
                 function(object, value) {
                   object@reducedExprDimension = value
                   validObject(object) # could add other checks
                   return(object)
                 })

#' Get or set \code{gene_id_type} from an SCData object
#' @name gene_id_type
#' @rdname gene_id_type
#' @param object An \code{\link{SCData}} object.
#' @param value Value to be assigned to corresponding object.
#'
#' @return gene id type string
#' @author Luyi Tian
#'
#' @export
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
