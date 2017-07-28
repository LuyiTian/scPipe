
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
#' @importFrom stats prcomp
#' @importFrom destiny DiffusionMap
#' @importFrom Rtsne Rtsne
#' 
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
#' scd = calQCMetrics(scd)
#' scd = detect_outlier(scd)
#' scd_afterqc = remove_outliers(scd)
#' scd_afterqc = sc_dim_reduction(scd_afterqc,n=2)
#'
sc_dim_reduction = function(scd,
                            k=20,
                            n=2) {
  # check format:
  if (is(scd, "SCData")) {
    exprs_mat <- switch(scd@useForExprs,
                        exprs = exprs(scd),
                        tpm = tpm(scd),
                        cpm = cpm(scd),
                        fpkm = fpkm(scd),
                        counts = counts(scd))
  }
  else if (is.matrix(scd)) {
    exprs_mat = scd
  }
  else {
    stop("scd must be an SCData object or a matrix.")
  }

  pca_out = prcomp(t(exprs_mat))
  # k = pca_permutation(pca_out, n=10) TODO
  dif = DiffusionMap(pca_out$x[, 1:k], n.eigs = k)
  tsne_out = Rtsne(dif@eigenvectors, dim=n)
  reduced_dim = tsne_out$Y
  rownames(reduced_dim) = colnames(exprs_mat)
  colnames(reduced_dim) = paste("Dim", 1:n, sep="")
  DimReducedExpr(scd) = reduced_dim
  return(scd)
}
