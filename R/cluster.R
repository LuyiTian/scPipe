
#' dimensionality reduction based on PCA, diffusion map and tsne
#'
#' @param scd an SCData object which contains gene counting matrix or
#' a gene counting matrix
#' @param k the number of top PCs in PCA to keep for following analysis
#' @param n the number of dimensions after
#' @description TODO
#'
#' @return an SCData object contains \code{reducedExprDimension} or
#' a matrix contains t-SNE output
#'
#' @import destiny Rtsne
#' @export
#' @examples TODO
#'
sc_dim_reduction = function(scd,
                            k=20,
                            n=2){
  # check format:
  if (is(scd, "SCData")){
    exprs_mat <- switch(scd@useForExprs,
                        exprs = exprs(scd),
                        tpm = tpm(scd),
                        cpm = cpm(scd),
                        fpkm = fpkm(scd),
                        counts = counts(scd))
  }
  else if (is.matrix(scd)){
    exprs_mat = scd
  }
  else{
    stop("scd must be an SCESet object or a matrix.")
  }

  pca_out = prcomp(t(exprs_mat))
  # k = pca_permutation(pca_out,n=10) TODO
  dif = DiffusionMap(pca_out$x[,1:k],n.eigs = k)
  tsne_out = Rtsne(dif@eigenvectors,dim=n)
  reduced_dim = tsne_out$Y
  rownames(reduced_dim) = colnames(exprs_mat)
  colnames(reduced_dim) = paste("Dim",1:n,sep="")


}
