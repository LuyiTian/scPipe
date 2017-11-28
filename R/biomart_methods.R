# biomaRt_func.R

#' Get genes related to certain GO terms from biomart database
#'
#' @param returns the gene id which is set as return. default to be ensembl id
#' A possible list of attributes can be retrieved using the
#' function \code{listAttributes} from \code{biomaRt} package. the commonly used
#' id types are `external_gene_name`, `ensembl_gene_id` or `entrezgene`.
#' @param dataset Dataset you want to use. List of possible datasets can be
#' retrieved using the function \code{listDatasets} from \code{biomaRt} package.
#' @param go a vector of GO terms
#' @details Get genes related to certain GO terms from biomart database
#'
#' @return a vector of gene ids.
#'
#' @import biomaRt
#' @importFrom biomaRt useDataset getBM
#'
#' @export
#' @examples
#' # get all genes under GO term GO:0005739 in mouse, return ensembl gene id
#' get_genes_by_GO(returns="ensembl_gene_id",
#'     dataset="mmusculus_gene_ensembl",
#'     go=c('GO:0005739'))
#'
get_genes_by_GO <- function(returns="ensembl_gene_id",
                            dataset="mmusculus_gene_ensembl",
                            go=NULL) {
  if (is.null(go)) {
    stop("must provide GO term. (i.e go=c('GO:0005739'))")
  }
  mart = tryCatch({mart <- useDataset(dataset, useMart("ensembl")) },
           error = function(e){
             cat(paste0("cannot connect to the ensembl database. ERROR:\n", e ))
             return(c())
           })
  if(!is(mart,"Mart")){
    return(c())
  }
  G_list <- getBM(filters="go",
                  attributes=c(returns),
                  values=go,
                  mart=mart)
  return(G_list[, returns])
}


#' convert the gene ids of a SingleCellExperiment object
#'
#' @param sce a SingleCellExperiment object
#' @param returns the gene id which is set as return. default to be `external_gene_name`
#' A possible list of attributes can be retrieved using the
#' function \code{listAttributes} from \code{biomaRt} package. the commonly used
#' id types are `external_gene_name`, `ensembl_gene_id` or `entrezgene`.
#' @param all logic. for genes that cannot covert to new gene id, keep them with the old
#' id or delete them. the default is keep them.
#' @details convert the gene id of all datas in the SingleCellExperiment object
#'
#' @return sce with converted id
#'
#' @importFrom biomaRt useDataset getBM
#' @importFrom utils head
#' @importFrom methods is
#'
#' @export
#' @examples
#' # the gene id in example data are `external_gene_name`
#' # the following example will convert it to `external_gene_name`.
#' data("sc_sample_data")
#' data("sc_sample_qc")
#' sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
#' organism(sce) = "mmusculus_gene_ensembl"
#' gene_id_type(sce) = "ensembl_gene_id"
#' QC_metrics(sce) = sc_sample_qc
#' demultiplex_info(sce) = cell_barcode_matching
#' UMI_dup_info(sce) = UMI_duplication
#' head(rownames(sce))
#' sce = convert_geneid(sce, return="external_gene_name")
#' head(rownames(sce))
#'
convert_geneid <- function(sce,
                           returns="external_gene_name",
                           all=TRUE) {
  sce = validObject(sce)
  if (returns == gene_id_type(sce)) {
    stop("SingleCellExperiment already has genes in such id type. (gene_id_type(sce) == returns)")
  }

  organism = organism(sce)
  if (organism == "NA") {
    print("organism not provided.")
    return(sce)
  }
  mart = tryCatch({mart <- useDataset(organism, useMart("ensembl")) },
           error = function(e){
             cat(paste0("cannot connect to the ensembl database. ERROR:\n", e ))
    return(sce)
  })
  if(!is(mart,"Mart")){
    return(sce)
  }
  
  G_list <- getBM(filters=gene_id_type(sce), attributes=c(gene_id_type(sce), returns, "description"), values=rownames(sce), mart=mart)

  G_list <- G_list[match(rownames(sce), G_list[, gene_id_type(sce)]), ]
  na_num <- sum(is.na(G_list[, returns]))
  dup_ids <- duplicated(G_list[, returns]) | duplicated(G_list[, returns], fromLast=TRUE)
  dup_num <- (sum(dup_ids)-na_num)/2
  print(paste0("number of NA in new gene id: ", na_num, ". duplicated id: ", dup_num))
  if (dup_num>0) {
    print("first 5 duplicated:")
    print(head(G_list[dup_ids & !(is.na(G_list[, returns])), ]))
  }
  G_list[, returns][dup_ids] <- NA
  if (all | (na_num+dup_num==0)) {
    # replace NA with old id
    G_list[, returns][is.na(G_list[, returns])] = rownames(sce)[is.na(G_list[, returns])]
    if (!(gene_id_type(sce) %in% colnames((sce@int_elementMetadata)))) {
      sce@int_elementMetadata[, gene_id_type(sce)] = rownames(sce)
    }
    if (!(returns %in% colnames(sce@int_elementMetadata))) {
      sce@int_elementMetadata[, returns] = G_list[, returns]
    }
    if (!("description" %in% colnames(sce@int_elementMetadata))) {
      sce@int_elementMetadata[, "description"] <- G_list[, "description"]
    }
    rownames(sce) <- G_list[, returns]
  }
  else {
    G_list <- G_list[!is.na(G_list[, returns]), ]
    sce <- sce[!is.na(G_list[, returns]), ]
    if (!(gene_id_type(sce) %in% colnames(sce@int_elementMetadata))) {
      sce@int_elementMetadata[, gene_id_type(sce)] <- rownames(sce)
    }
    if (!(returns %in% colnames(sce@int_elementMetadata))) {
      sce@int_elementMetadata[, returns] <- G_list[, returns]
    }
    if (!("description" %in% colnames(sce@int_elementMetadata))) {
      sce@int_elementMetadata[, "description"] <- G_list[, "description"]
    }
    rownames(sce) <- G_list[, returns]
  }
  gene_id_type(sce) <- returns
  return(sce)
}
