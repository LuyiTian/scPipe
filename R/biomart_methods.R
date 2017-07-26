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
#' #TODO
#'
get_genes_by_GO <- function(returns="ensembl_gene_id",
                            dataset="mmusculus_gene_ensembl",
                            go=NULL) {
  if (is.null(go)) {
    stop("must provide GO term. (i.e go=c('GO:0005739'))")
  }
  mart <- useDataset(dataset, useMart("ensembl"))
  G_list <- getBM(filters="go",
                  attributes=c(returns),
                  values=go,
                  mart=mart)
  return(G_list[, returns])
}


#' convert the gene ids of a SCData object
#'
#' @param scd an SCData object
#' @param returns the gene id which is set as return. default to be `external_gene_name`
#' A possible list of attributes can be retrieved using the
#' function \code{listAttributes} from \code{biomaRt} package. the commonly used
#' id types are `external_gene_name`, `ensembl_gene_id` or `entrezgene`.
#' @param all logic. for genes that cannot covert to new gene id, keep them with the old
#' id or delete them. the default is keep them.
#' @details convert the gene id of all datas in SCData object
#'
#' @return scd with converted id
#'
#' @importFrom biomaRt useDataset getBM
#' @importFrom utils head
#' @importFrom Biobase ExpressionSet fData<-
#'
#' @export
#' @examples
#' #TODO
#'
convert_geneid <- function(scd,
                           returns="external_gene_name",
                           all=TRUE) {
  if (!is(scd, "SCData")) {
    stop("scd must be an SCData object.")
  }
  if (returns == gene_id_type(scd)) {
    stop("SCData already in this id type. (scd@gene_id_type == returns)")
  }

  species <- organism.SCData(scd)
  if(species == "NA"){
    print("species not provided.")
    return(scd)
  }
  mart <- useDataset(species, useMart("ensembl"))
  G_list <- getBM(filters=gene_id_type(scd), attributes=c(gene_id_type(scd), returns, "description"), values=rownames(scd), mart=mart)

  G_list <- G_list[match(rownames(scd), G_list[, gene_id_type(scd)]), ]
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
    G_list[, returns][is.na(G_list[, returns])] <- rownames(scd)[is.na(G_list[, returns])]
    if (!(gene_id_type(scd) %in% colnames(Biobase::fData(scd)))) {
      Biobase::fData(scd)[, gene_id_type(scd)] <- rownames(scd)
    }
    if (!(returns %in% colnames(Biobase::fData(scd)))) {
      Biobase::fData(scd)[, returns] <- G_list[, returns]
    }
    if (!("description" %in% colnames(Biobase::fData(scd)))) {
      Biobase::fData(scd)[, "description"] <- G_list[, "description"]
    }
    rownames(scd) <- G_list[, returns]
  }
  else{
    G_list <- G_list[!is.na(G_list[, returns]), ]
    scd <- scd[!is.na(G_list[, returns]), ]
    if (!(gene_id_type(scd) %in% colnames(Biobase::fData(scd)))) {
      Biobase::fData(scd)[, gene_id_type(scd)] <- rownames(scd)
    }
    if (!(returns %in% colnames(Biobase::fData(scd)))) {
      Biobase::fData(scd)[, returns] <- G_list[, returns]
    }
    if (!("description" %in% colnames(Biobase::fData(scd)))) {
      Biobase::fData(scd)[, "description"] <- G_list[, "description"]
    }
    rownames(scd) <- G_list[, returns]
  }
  gene_id_type(scd) <- returns
  return(scd)
}
