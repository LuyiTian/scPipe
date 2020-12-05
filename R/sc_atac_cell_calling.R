#' sc_atac_bam_tagging()
#'
#' @return 
#'
#' @examples
#' \dontrun{
#' 
#' 
#' }
#'
#' @export
#'
sc_atac_cell_calling = function (mat, cell_calling, output_folder, genome_size = NULL, qc_per_bc_file = NULL){
  
  if(cell_calling == 'emptydrops'){
    
    library(DropletUtils)
    library(data.table)
    library(Matrix)
    
    set.seed(2019)
    cell.out <- emptyDrops2(mat)
    
    filter.out <- cell.out[S4Vectors::complete.cases(cell.out), ]
    
    saveRDS(filter.out, file = paste0(output_folder, '/EmptyDrop_obj.rds'))
    
    if(length(filter.out$FDR) > 0){
      cat("Empty filter.out\n")
      filter.out = filter.out[filter.out$FDR <= fdr, ]
    }
    
    select.cells = rownames(filter.out)
    
    out_mat = mat[, colnames(mat) %in% select.cells]
    barcodes     <- colnames(out_mat)
    features     <- rownames(out_mat)
    
    if(length(filter.out$FDR) > 0){
      cat("Empty filter.out\n")
      writeMM(out_mat, file = paste0(output_folder, '/matrix.mtx'))
      write.table(barcodes, file = paste0(output_folder, '/non_empty_barcodes.txt'), sep = '\t',
                  row.names = FALSE, quote = FALSE, col.names = FALSE)
      write.table(features, file = paste0(output_folder, '/non_empty_features.txt'), sep = '\t',
                  row.names = FALSE, quote = FALSE, col.names = FALSE)
    }
    
    
  } # end emptydrops
  
  
  if(cell_calling == 'cellranger'){
    cellranger_cell_caller(mat, output_folder, genome_size, qc_per_bc_file)
  } # end cellranger
  
  
  
  if(cell_calling == 'filter'){
    filter_barcodes(mtx = mat, output_folder = output_folder)
  } # end filter
  
  
}








testEmptyDrops2 = function (m, lower = 100, niters = 10000, test.ambient = FALSE, 
                            ignore = NULL, alpha = NULL, BPPARAM = SerialParam()) 
{
  discard <- rowSums(m) == 0
  m <- m[!discard, , drop = FALSE]
  ncells <- ncol(m)
  umi.sum <- as.integer(round(colSums(m)))
  ambient <- umi.sum <= lower
  ambient.cells <- m[, ambient, drop = F]
  ambient.prof <- rowSums(ambient.cells)
  if (sum(ambient.prof) == 0) {
    stop("no counts available to estimate the ambient profile")
  }
  ambient.prop <- edgeR::goodTuringProportions(ambient.prof)
  if (!test.ambient) {
    keep <- !ambient
  }
  else {
    keep <- umi.sum > 0L
  }
  if (!is.null(ignore)) {
    keep <- keep & umi.sum > ignore
  }
  obs <- m[, keep, drop = FALSE]
  obs.totals <- umi.sum[keep]
  if (is.null(alpha)) {
    alpha <- DropletUtils:::.estimate_alpha(m[, ambient, drop = FALSE], 
                                            ambient.prop, umi.sum[ambient])
  }
  obs.P <-DropletUtils:::.compute_multinom_prob_data(obs, ambient.prop, alpha = alpha)
  rest.P <-DropletUtils:::.compute_multinom_prob_rest(obs.totals, alpha = alpha)
  n.above <-DropletUtils:::.permute_counter(totals = obs.totals, probs = obs.P, 
                                            ambient = ambient.prop, iter = niters, BPPARAM = BPPARAM, 
                                            alpha = alpha)
  limited <- n.above == 0L
  pval <- (n.above + 1)/(niters + 1)
  all.p <- all.lr <- all.exp <- rep(NA_real_, ncells)
  all.lim <- rep(NA, ncells)
  all.p[keep] <- pval
  all.lr[keep] <- obs.P + rest.P
  all.lim[keep] <- limited
  output <- DataFrame(Total = umi.sum, LogProb = all.lr, PValue = all.p, 
                      Limited = all.lim, row.names = colnames(m))
  metadata(output) <- list(lower = lower, niters = niters, 
                           ambient = ambient.prop, alpha = alpha)
  output
}





emptyDrops2 = function (m, lower = 100, retain = NULL, barcode.args = list(), 
                        ...) 
{
  m <- DropletUtils:::.rounded_to_integer(m)
  stats <- testEmptyDrops2(m, lower = lower, ...)
  tmp <- stats$PValue
  if (is.null(retain)) {
    br.out <- do.call(barcodeRanks, c(list(m, lower = lower), 
                                      barcode.args))
    retain <- metadata(br.out)$knee
  }
  always <- stats$Total >= retain
  tmp[always] <- 0
  metadata(stats)$retain <- retain
  stats$FDR <- p.adjust(tmp, method = "BH")
  return(stats)
}
















cellranger_cell_caller = function(mat, output_folder, genome_size, qc_per_bc_file){
  # https://github.com/wbaopaul/scATAC-pro/blob/master/scripts/src/cellranger_cell_caller.R
  
  library(data.table)
  library(Matrix)
  library(flexmix)
  library(countreg)  ##install.packages("countreg", repos="http://R-Forge.R-project.org")
  
  # args = commandArgs(T)
  # 
  # input_mtx_file = args[1]
  # output_folder = args[2]
  # genome_size = as.numeric(args[3])
  # qc_per_bc_file = args[4]
  # 
  # 
  # ## read matrix data
  # input_mtx_dir = dirname(input_mtx_file)
  # mat = readMM(input_mtx_file)
  # 
  # barcodes = fread(paste0(input_mtx_dir, '/barcodes.txt'), header = F)
  # 
  # colnames(mat) = barcodes$V1
  
  
  ## filter barcodes by frac_in_peak
  qc_per_bc = fread(qc_per_bc_file)
  peak_cov_frac = min(0.05, nrow(mat) * 1000/genome_size)
  qc_sele_bc = qc_per_bc[frac_peak >= peak_cov_frac]
  
  ## subtract counts due to contamination (rate 0.02)
  CN = max(1, round(median(qc_sele_bc$total_frags)* 0.02))
  qc_sele_bc[, 'total_frags' := total_frags -CN]
  qc_sele_bc = qc_sele_bc[total_frags >= 0]
  
  ## fit two NB mixture model & using signal to noisy ratio to select cells
  n_in_peak = Matrix::colSums(mat)
  n_in_peak = n_in_peak[names(n_in_peak) %in% qc_sele_bc$bc]
  flexmix(n_in_peak ~ 1, k = 2, model = FLXMRnegbin(theta = 1))
  fm0 <- flexmix(n_in_peak ~ 1, k = 2, model = FLXMRnegbin())
  prob1 = posterior(fm0)[, 1]
  prob2 = posterior(fm0)[, 2]
  mus = parameters(fm0)[1, ]
  
  if(mus[1] > mus[2]){
    #odd = prob1/prob2
    odd = prob1
  }else{
    #odd = prob2/prob1
    odd = prob2
  }
  aa = which(odd == 1)
  select.cells = names(n_in_peak)[aa]
  length(select.cells)
  
  out_mat = mat[, colnames(mat) %in% select.cells]
  barcodes = colnames(out_mat)
  # dim(out_mat)
  
  
  # system(paste('mkdir -p', output_folder))
  writeMM(out_mat, file = paste0(output_folder, '/matrix.mtx'))  
  write.table(barcodes, file = paste0(output_folder, '/barcodes.txt'), 
              sep = '\t', row.names = F, quote = F, col.names = F)
  # system(paste0('cp ', input_mtx_dir, '/features.txt ',  output_folder, '/'))
}











filter_barcodes = function(
  mtx, 
  output_folder,
  bc_stat_file = NULL,
  min_uniq_frags = 3000,
  max_uniq_frags = 50000,
  min_frac_peak = 0.05,
  min_frac_tss = 0,
  min_frac_enhancer = 0,
  min_frac_promoter = 0,
  max_frac_mito = 0.2){
  
  # https://github.com/wbaopaul/scATAC-pro/blob/master/scripts/src/filter_barcodes.R
  
  ## call cell by filtering barcodes given some qc stats cutoffs
  library(data.table)
  library(Matrix)
  
  qc_bc_stat = fread(bc_stat_file)
  
  # Change names of variables from optparse package
  cut.min.frag = min_uniq_frags
  cut.max.frag = max_uniq_frags
  cut.mito = max_frac_mito
  cut.peak = min_frac_peak
  cut.tss = min_frac_tss
  cut.promoter = min_frac_promoter
  cut.enh = min_frac_enhancer
  
  qc_sele = qc_bc_stat[total_frags >= cut.min.frag & total_frags <= cut.max.frag &
                         frac_mito <= cut.mito &
                         frac_peak >= cut.peak &
                         frac_tss >= cut.tss &
                         frac_promoter >= cut.promoter &
                         frac_enhancer >= cut.enh]
  
  # mtx = readMM(mtx_file)
  # input_mtx_dir = dirname(mtx_file)
  # colnames(mtx) = fread(paste0(input_mtx_dir, '/barcodes.txt'), header = F)$V1
  mtx = mtx[, colnames(mtx) %in% qc_sele$bc]
  saveRDS(mtx, file = paste0(output_folder, '/matrix.rds'))
  
  # mtx.dir = dirname(mtx_file)
  # system(paste('mkdir -p', output_folder))
  writeMM(mtx, file = paste0(output_folder, '/matrix.mtx'))
  write.table(colnames(mtx), file = paste0(output_folder, '/barcodes.txt'), 
              sep = '\t', row.names = F, quote = F, col.names = F)
  # system(paste0('cp ', dirname(mtx_file), '/features.txt ', output_folder, '/features.txt'))
}




