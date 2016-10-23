
#' clustering and visualization based on PCA, diffusion map and tsne
#' 
#' @param scd an SCData object, contains gene counting matrix
#' @param k the number of top PCs in PCA to keep for following analysis
#' @description TODO
#' 
#' @return an SCData object contains \code{reducedExprDimension}
#' 
#' @import destiny
#' @export
#' @examples TODO
#' 
sc_cluster = function(scd, k){
  exprs_mat <- switch(scd@useForExprs,
                      exprs = exprs(scd),
                      tpm = tpm(scd),
                      cpm = cpm(scd),
                      fpkm = fpkm(scd),
                      counts = counts(scd))
  pca_out = prcomp(t(exprs_mat))
}