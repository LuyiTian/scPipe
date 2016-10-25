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
#'
#' @export
#' @examples
#' TODO
#'
get_genes_by_GO = function(returns="ensembl_gene_id",
                dataset="mmusculus_gene_ensembl",
                go=NULL){
  if(is.na(go)){
    stop("must provide GO term. (i.e go=c('GO:0005739'))")
  }
  mart <- useDataset(dataset, useMart("ensembl"))
  G_list <- getBM(filters= "go_id",
                  attributes= c(returns),
                  values=go,
                  mart=mart)
  return(G_list[,returns])

}



