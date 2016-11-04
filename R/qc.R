
#' Detect outliers based on robust linear regression of QQ plot
#'
#' @param x a vector of mahalanobis distance
#' @param df degree of freedom for chi-square distribution
#' @param conf confidence for linear regression
#'
#' @import MASS
#'
#' @return cell names of outliers
#' or both (`both`)
.qq_outliers_robust = function(x, df, conf){
  n <- length(x)
  P = ppoints(n)
  z = qchisq(P, df = df)
  ord.x <- x[order(x)]
  coef <- coef(rlm(ord.x ~ z,maxit = 200))
  a <- coef[1]
  b <- coef[2]
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dchisq(z, df))*sqrt(P*(1 - P)/n)
  fit.value <- a + b*z
  upper <- fit.value + zz*SE
  thr = max(ord.x[ord.x < upper])
  return(names(x[x>thr]))
}


#' Detect outliers based on QC metrics
#'
#' @param scd an SCData object containing expression values and
#' QC metrics.
#' @param comp the number of component used in GMM. can be one or two.
#' @param plot logical, whether to generate a pairwise plot after
#' @param type only looking at low quality cells (`low`) or possible doublets (`high`)
#' or both (`both`)
#' @param conf confidence interval for linear regression at lower and upper tails.
#' Usually this is smaller for lower tail because we hope to pick out more
#' low quality cells than doublets.
#' @details detect outlier using Mahalanobis distances TODO
#'
#' @return an updated SCData object with an outlier column in \code{QualityControlInfo}
#'
#' @import mclust robustbase
#'
#' @export
#' @examples
#' TODO
#'
detect_outlier = function(scd,
                          comp=1,
                          plot=TRUE,
                          type = c("both","low","high"),
                          conf = c(0.9,0.99)){
  # check format:
  if (is(scd, "SCData")){
    x = pData(QC_metrics(scd))
    x = x[,!("outliers" == colnames(x))]
  }
  else{
    stop("scd must be an SCESet object.")
  }

  dist = mahalanobis(x, center=colMeans(x), cov=cov(x))
  keep = !(dist>qchisq(0.99, ncol(x)))
  mod = Mclust(x[keep,],
               G=comp,
               modelNames="EEE")
  #print(plot(mod,what="classification"))
  if(comp == 1){
    covr <- covMcd(x, alpha=0.7)
    dist = mahalanobis(x,
                       center=covr$center,
                       cov=covr$cov)
    mean_diff = sign(t(x)-covr$center)
    QC_sign = c(-1,1)[as.factor(apply(mean_diff,2,function(t){sum(t)>0}))]
    neg_dist = dist[QC_sign == -1]
    pos_dist = dist[QC_sign == 1]
    if(type == "both"){
      outlier_cells = .qq_outliers_robust(neg_dist, ncol(x), conf[1])
      outlier_cells = c(outlier_cells,
                        .qq_outliers_robust(pos_dist, ncol(x), conf[2]))
    }
    else if(type == "low"){
      outlier_cells = .qq_outliers_robust(neg_dist, ncol(x), conf[1])
    }
    else if(type == "high"){
      outlier_cells = .qq_outliers_robust(pos_dist, ncol(x), conf[2])
    }
  }
  else{
    ord_fst = c(1:comp)[order(mod$parameters$mean[1,],decreasing = TRUE)]
    poor_comp = ord_fst[2:comp]
    good_comp = ord_fst[1]
    keep1 = rep(TRUE,nrow(x))
    keep1[keep][mod$classification %in% poor_comp] = FALSE
    keep1[!keep] = FALSE
    sub_x = x[keep1,]
    covr <- covMcd(sub_x, alpha=0.7)
    sub_dist = mahalanobis(sub_x,
                           center=covr$center,
                           cov=covr$cov)

    mean_diff = sign(t(sub_x)-covr$center)
    QC_sign = c(-1,1)[as.factor(apply(mean_diff,2,function(t){sum(t)>0}))]
    neg_dist = sub_dist[QC_sign == -1]
    pos_dist = sub_dist[QC_sign == 1]
    outlier_cells = .qq_outliers_robust(neg_dist, ncol(sub_x), conf[1])
    outlier_cells = c(outlier_cells,
                      .qq_outliers_robust(pos_dist, ncol(sub_x), conf[2]))
    outlier_cells = c(outlier_cells,rownames(x[!keep1,]))
    if(!(type == "both")){
      mean_diff = sign(t(x)-mod$parameters$mean[,good_comp])
      QC_sign = c(-1,1)[as.factor(apply(mean_diff,2,function(t){sum(t)>0}))]
      if (type == "low"){
        outlier_cells = rownames(x)[(rownames(x) %in% outlier_cells) & (QC_sign == -1)]
      }
      else if (type == "high"){
        outlier_cells = rownames(x)[(rownames(x) %in% outlier_cells) & (QC_sign == 1)]
      }
    }
  }
  outliers = as.factor(rownames(x)  %in% outlier_cells)
  x$outliers = outliers
  QC_metrics(scd) = x
  return(scd)
}



#' get QC metrics using gene counting matrix
#'
#' @param scd an SCData object containing count
#' @details get QC metrics using gene counting matrix
#' the QC statistics added are
#' `total_count_per_cell`,
#' `non_mt_percent`: 1- percent of mitochondrial gene counts
#' mitochondrial genes are retrived by GO term GO:0005739
#' `exon_to_ERCC_ratio`: ratio of exon counts to ERCC counts
#' `non_ribo_percent`: 1- percent of ribosomal gene counts
#' ribosomal genes are retrived by GO term GO:0005840
#' @return no return
#'
#'
#' @export
#' @examples
#' TODO
#'
calculate_QC_metrics = function(scd){
  if (is(scd, "SCData")){
    exprs_mat <- switch(scd@useForExprs,
                        exprs = exprs(scd),
                        tpm = tpm(scd),
                        cpm = cpm(scd),
                        fpkm = fpkm(scd),
                        counts = counts(scd))
  }
  else{
    stop("require a SCData object.")
  }
  QC_met = pData(QC_metrics(scd))
  # get ERCC ratio
  spikein = fData(scd)$isSpike
  exon_count = colSums(exprs_mat[!spikein,])
  if(any(spikein)){
    ERCC_count = colSums(exprs_mat[spikein,])
    QC_met$exon_to_ERCC_ratio = exon_count/ERCC_count
  }
  else{
    print("cannot detect ERCC Spikeins from data. skip `exon_to_ERCC_ratio`.")
  }
  # get mt percentage
  mt_genes = get_genes_by_GO(returns=gene_id_type(scd),
                             dataset=organism(scd),
                             go=c("GO:0005739"))
  mt_count = colSums(exprs_mat[rownames(exprs_mat) %in% mt_genes,])
  QC_met$non_mt_percent = (exon_count-mt_count)/exon_count

  # get ribosomal percentage
  ribo_genes = get_genes_by_GO(returns=gene_id_type(scd),
                             dataset=organism(scd),
                             go=c("GO:0005840"))
  ribo_count = colSums(exprs_mat[rownames(exprs_mat) %in% ribo_genes,])
  QC_met$non_ribo_percent = (exon_count-ribo_count)/exon_count
  QC_metrics(scd) = QC_met
  return(scd)
}


#' plot QC statistics for SCData object
#' @param scd an SCData object
#' @import GGally ggplot2
#' @export
#'
plotQC = function(scd){
  if (is(scd, "SCData")){
    x = pData(QC_metrics(scd))
  }
  else{
    stop("scd must be an SCESet object.")
  }

  if ("outliers" %in% colnames(x)){
    return(ggpairs(x,
            mapping=ggplot2::aes(colour = outliers)))
  }
  else{
    return(ggpairs(x))
  }
}

#' remove outliers for SCData
#' @param scd an SCData object
#' @export
#'
remove_outliers = function(scd){
  if (!is(scd, "SCData")){
    stop("scd must be an SCESet object.")
  }
  if(!("outliers" %in% colnames(pData(QC_metrics(scd))))){
    stop("no outlier information. please run `detect_outlier()` first.")
  }
  out_cell = pData(QC_metrics(scd))$outliers == FALSE
  return(scd[,out_cell])
}



