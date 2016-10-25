
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
#' QC metrics. or just a matrix contains QC metrics.
#' @param type only looking at low quality cells (`low`) or possible doublets (`high`)
#' or both (`both`)
#' @param conf confidence interval for linear regression at lower and upper tails.
#' Usually this is smaller for lower tail because we hope to pick out more
#' low quality cells than doublets.
#' @details detect outlier using Mahalanobis distances TODO
#'
#' @return boolean vector of outliers
#'
#' @import mclust
#'
#' @export
#' @examples
#' TODO
#'
detect_outlier = function(scd,
                          type = c("both","low","high"),
                          conf = c(0.9,0.99)){
  # check format:
  if (is(scd, "SCData")){
    x = pData(QC_metrics(scd))
  }
  else if (is.matrix(scd)){
    x = scd
  }
  else{
    stop("scd must be an SCESet object or a matrix.")
  }

  dist = mahalanobis(x, center=colMeans(x), cov=cov(x))
  keep = !(dist>qchisq(0.99, ncol(x)))
  mod = Mclust(x[keep,],
               G=1:2,
               modelNames="EEE")
  print(summary(mod))
  print(plot(mod,what="classification"))
  if(mod$G == 1){
    mean_diff = sign(t(x)-colMeans(x))
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
    dif = sum(sign(mod$parameters$mean[,1] - mod$parameters$mean[,2]))
    if (dif < 0){
      poor_comp = 1
      good_comp = 2
    }
    else if(dif > 0){
      poor_comp = 2
      good_comp = 1
    }
    else{
      print(summary(mod))
      stop("we cannot decide which component in mixture model corresponding to the good quality cells.")
    }
    keep1 = rep(TRUE,nrow(x))
    keep1[keep][mod$classification == poor_comp] = FALSE
    sub_x = x[keep1,]
    sub_dist = mahalanobis(sub_x,
                           center=mod$parameters$mean[,good_comp],
                           cov=mod$parameters$variance$sigma[,,good_comp])

    mean_diff = sign(t(sub_x)-mod$parameters$mean[,good_comp])
    QC_sign = c(-1,1)[as.factor(apply(mean_diff,2,function(t){sum(t)>0}))]
    neg_dist = sub_dist[QC_sign == -1]
    pos_dist = sub_dist[QC_sign == 1]
    outlier_cells = .qq_outliers_robust(neg_dist, ncol(sub_x), conf[1])
    outlier_cells = c(outlier_cells,
                      .qq_outliers_robust(pos_dist, ncol(sub_x), conf[2]))
    outlier_cells = c(outlier_cells,rownames(x[!keep1,]))
    if(type == "both"){
      return(outlier_cells)
    }
    else{
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
  return(outlier_cells)
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

  # get ERCC ratio
  spikein = fData(scd)$isSpike
  exon_count = colSums(exprs_mat[!spikein,])
  if(any(spikein)){
    ERCC_count = colSums(exprs_mat[spikein,])
    pData(QC_metrics(scd))$exon_to_ERCC_ratio = exon_count/ERCC_count
  }
  else{
    print("cannot detect ERCC Spikeins from data. skip `exon_to_ERCC_ratio`.")
  }

  # get mt percentage
  mt_genes = get_genes_by_GO(returns=gene_id_type(scd),
                             dataset=organism(scd),
                             go=c("GO:0005739"))
  mt_count = colSums(exprs_mat[rownames(exprs_mat) %in% mt_genes,])
  pData(QC_metrics(scd))$non_mt_percent = (exon_count-mt_count)/exon_count

  # get ribosomal percentage
  ribo_genes = get_genes_by_GO(returns=gene_id_type(scd),
                             dataset=organism(scd),
                             go=c("GO:0005840"))
  ribo_count = colSums(exprs_mat[rownames(exprs_mat) %in% ribo_genes,])
  pData(QC_metrics(scd))$non_ribo_percent = (exon_count-ribo_count)/exon_count
}
