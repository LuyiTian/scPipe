
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
#' @param x an SCESet object containing expression values and
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
#' @references P. Filzmoser, R.G. Garrett, and C. Reimann.
#' Multivariate outlier detection in exploration geochemistry.
#' \emph{Computers & Geosciences}, 31:579-587, 2005.
#'
detect_outlier = function(x, 
                          type = c("both","low","high"),
                          conf = c(0.9,0.99)){
  dist = mahalanobis(x, center=colMeans(x), cov=cov(x))
  keep = !(dist>qchisq(0.99, ncol(x)))
  mod = Mclust(x[keep,],
               G=1:2,
               modelNames="EEE")
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