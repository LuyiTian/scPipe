#######################################
# Cell Calling on a matrix (scATAC-Seq)
#######################################

#' @name sc_atac_cell_calling
#' @title identifying true vs empty cells
#' @description the methods to call true cells are of various ways.
#' implement (i.e. \code{filtering} from \code{scATAC-Pro} and \code{cellranger approach}).
#' @param mat the feature by cell matrix. 
#' @param cell_calling the cell calling approach, possible options are "emptydrops" , "cellranger" and "filter".
#' @param output_folder output directory for the cell called matrix.
#' @param genome_size genome size for the data in feature by cell matrix.
#' @param qc_per_bc_file quality per barcode file for the barcodes in the matrix if using the \code{cellranger} or \code{filter} options.
#'
#' @param lower the lower threshold for the data if using the \code{emptydrops} function for cell calling.
#' 
#' @param min_uniq_frags The minimum number of required unique fragments required for a cell (used for \code{filter} cell calling)
#' @param max_uniq_frags The maximum number of required unique fragments required for a cell (used for \code{filter} cell calling)
#' @param min_frac_peak The minimum proportion of fragments in a cell to overlap with a peak (used for \code{filter} cell calling)
#' @param min_frac_tss The minimum proportion of fragments in a cell to overlap with a tss (used for \code{filter} cell calling)
#' @param min_frac_enhancer The minimum proportion of fragments in a cell to overlap with a enhancer sequence (used for \code{filter} cell calling)
#' @param min_frac_promoter The minimum proportion of fragments in a cell to overlap with a promoter sequence (used for \code{filter} cell calling)
#' @param max_frac_mito The maximum proportion of fragments in a cell that are mitochondrial (used for \code{filter} cell calling)
#'
#' @examples
#' \dontrun{
#' sc_atac_cell_calling <- function(mat, 
#'  cell_calling, 
#'  output_folder, 
#'  genome_size    = NULL, 
#'  qc_per_bc_file = NULL, 
#'  lower          = NULL)
#' }
#'
#'
#' @export
#'
sc_atac_cell_calling <- function(mat, 
                                 cell_calling, 
                                 output_folder, 
                                 genome_size    = NULL, 
                                 qc_per_bc_file = NULL, 
                                 lower          = NULL,
                                 min_uniq_frags = 3000,
                                 max_uniq_frags = 50000,
                                 min_frac_peak = 0.3,
                                 min_frac_tss = 0,
                                 min_frac_enhancer = 0,
                                 min_frac_promoter = 0.1,
                                 max_frac_mito = 0.15){
  
  cat("calling `", cell_calling, "` function for cell calling ... \n")

  selected_cells <- tryCatch(
    {
      if(cell_calling == "emptydrops") {
        sc_atac_emptydrops_cell_calling(mat = mat, output_folder = output_folder, lower = lower)
      }
      else if(cell_calling == 'cellranger') {
        sc_atac_cellranger_cell_calling(mat = mat, genome_size = genome_size, qc_per_bc_file = qc_per_bc_file)
      } 
    },
    error = function(e) {
      message(e)
      message("\nWill now default to the filter method with the default qc cutoffs.")
      message("min_uniq_frags = ", min_uniq_frags)
      message("max_uniq_frags = ", max_uniq_frags)
      message("min_frac_peak = ", min_frac_peak)
      message("min_frac_tss = ", min_frac_tss)
      message("min_frac_enhancer = ", min_frac_enhancer)
      message("min_frac_promoter = ", min_frac_promoter)
      message("max_frac_mito = ", max_frac_mito)
      message("Running the filter method...")
      sc_atac_filter_cell_calling(mtx = mat, 
                                  qc_per_bc_file = qc_per_bc_file,
                                  min_uniq_frags = min_uniq_frags,
                                  max_uniq_frags = max_uniq_frags,
                                  min_frac_peak = min_frac_peak,
                                  min_frac_tss = min_frac_tss,
                                  min_frac_enhancer = min_frac_enhancer,
                                  min_frac_promoter = min_frac_promoter,
                                  max_frac_mito = max_frac_mito)
    }
  )
  
  if(cell_calling == 'filter') {
    selected_cells <- sc_atac_filter_cell_calling(mtx = mat, 
                                                  qc_per_bc_file = qc_per_bc_file,
                                                  min_uniq_frags = min_uniq_frags,
                                                  max_uniq_frags = max_uniq_frags,
                                                  min_frac_peak = min_frac_peak,
                                                  min_frac_tss = min_frac_tss,
                                                  min_frac_enhancer = min_frac_enhancer,
                                                  min_frac_promoter = min_frac_promoter,
                                                  max_frac_mito = max_frac_mito)
  } 
  else if(isFALSE(cell_calling)) {
    cat("No cell calling method was selected.\n")
  }
  else if (!cell_calling %in% c("emptydrops", "cellranger", "filter")) { # no legitimate cell calling method chosen, so just return the original matrix
    cat(cell_calling, "was not an implemented cell calling method\n")
  }
  out_mat <- mat
  barcodes     <- colnames(out_mat)
  features     <- rownames(out_mat)
  if (cell_calling %in% c("emptydrops", "cellranger", "filter")) {
    # Only keep selected cells in matrix
    out_mat      <- mat[, colnames(mat) %in% selected_cells]
    barcodes     <- colnames(out_mat)
    features     <- rownames(out_mat)
    
    # cat("Number of called barcodes: ")
    # cat(length(selected_cells))
    # 
    # cat("\nNumber of columns: ")
    # cat(length(barcodes))
    # cat("\n")
    if (length(barcodes) == 0) {
      stop("No cells were called...")
    }
  }
  
  # Store output matrix
  # Matrix::writeMM(Matrix::Matrix(out_mat), file =file.path(output_folder, 'cell_called_matrix.mtx'))
  cat("cell called and stored in ", output_folder, "\n")
  utils::write.table(barcodes, file = file.path(output_folder, 'non_empty_barcodes.txt'), sep = '\t',
              row.names = FALSE, quote = FALSE, col.names = FALSE)
  utils::write.table(features, file = file.path(output_folder, 'non_empty_features.txt'), sep = '\t',
              row.names = FALSE, quote = FALSE, col.names = FALSE)

  
  return(out_mat)
}


# ############## Modified testEmptyDrops function ######################
# testEmptyDrops2 <- function (m, lower = 100, niters = 10000, test.ambient = FALSE, 
#                             ignore = NULL, alpha = NULL, BPPARAM = SerialParam()) 
# {
#   discard      <- rowSums(m) == 0
#   m            <- m[!discard, , drop = FALSE]
#   ncells       <- ncol(m)
#   umi.sum      <- as.integer(round(colSums(m)))
#   ambient      <- umi.sum <= lower
#   ambient.cells<- m[, ambient, drop = FALSE]
#   ambient.prof <- rowSums(ambient.cells)
#   
#   if (sum(ambient.prof) == 0) {
#     stop("no counts available to estimate the ambient profile")
#   }
#   ambient.prop <- edgeR::goodTuringProportions(ambient.prof)
#   if (!test.ambient) {
#     keep <- !ambient
#   } else {
#     keep <- umi.sum > 0L
#   }
#   if (!is.null(ignore)) {
#     keep <- keep & umi.sum > ignore
#   }
#   obs <- m[, keep, drop = FALSE]
#   obs.totals <- umi.sum[keep]
#   if (is.null(alpha)) {
#     alpha <- DropletUtils:::.estimate_alpha(m[, ambient, drop = FALSE], 
#                                             ambient.prop, umi.sum[ambient])
#   }
#   obs.P <-DropletUtils:::.compute_multinom_prob_data(obs, ambient.prop, alpha = alpha)
#   rest.P <-DropletUtils:::.compute_multinom_prob_rest(obs.totals, alpha = alpha)
#   n.above <-DropletUtils:::.permute_counter(totals = obs.totals, probs = obs.P, 
#                                             ambient = ambient.prop, iter = niters, BPPARAM = BPPARAM, 
#                                             alpha = alpha)
#   limited <- n.above == 0L
#   pval <- (n.above + 1)/(niters + 1)
#   all.p <- all.lr <- all.exp <- rep(NA_real_, ncells)
#   all.lim <- rep(NA, ncells)
#   all.p[keep] <- pval
#   all.lr[keep] <- obs.P + rest.P
#   all.lim[keep] <- limited
#   output <- DataFrame(Total = umi.sum, LogProb = all.lr, PValue = all.p, 
#                       Limited = all.lim, row.names = colnames(m))
#   metadata(output) <- list(lower = lower, niters = niters, 
#                            ambient = ambient.prop, alpha = alpha)
#   output
# }
# 
# ############## Rounded to integer function #####################
# .rounded_to_integer <- function(m, round=TRUE) {
#   if (round) {
#     m <- round(m)
#   }
#   m
# }
# 
# ############## Modified emptyDrops function ######################
# emptyDrops2 <- function (m, lower = 100, retain = -1, barcode.args = list(), 
#                         ...) 
# {
#   m <- .rounded_to_integer(m)
#   stats <- testEmptyDrops2(m, lower = lower, ...)
#   tmp <- stats$PValue
#   if (is.null(retain)) {
#     br.out <- do.call(barcodeRanks, c(list(m, lower = lower), 
#                                       barcode.args))
#     retain <- metadata(br.out)$knee
#   }
#   always <- stats$Total >= retain
#   tmp[always] <- 0
#   metadata(stats)$retain <- retain
#   stats$FDR <- p.adjust(tmp, method = "BH")
#   return(stats)
# }

#' @name sc_atac_cellranger_cell_calling
#' @title cellranger cell calling
#' @description use the cellranger cell calling algorithm
#' 
#' @param mat The input matrix
#' @param qc_per_bc_file A file containing qc statistics for each cell
#' @param genome_size The size of the genome
#' 
#' @importFrom data.table :=
#' 
#' @export
#' 
sc_atac_cellranger_cell_calling <- function(mat, qc_per_bc_file, genome_size){
  # https://github.com/wbaopaul/scATAC-pro/blob/master/scripts/src/cellranger_cell_caller.R
  # https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/algorithms/overview
  
  frac_peak <- total_frags <- NULL
  
  # first filter barcodes by frac_in_peak
  qc_per_bc <- fread(qc_per_bc_file)
  peak_cov_frac <- min(0.05, nrow(mat) * 1000/genome_size)
  qc_sele_bc <- qc_per_bc[frac_peak >= peak_cov_frac]
  
  # subtract counts due to contamination (rate 0.02)
  CN <- max(1, round(stats::median(qc_sele_bc$total_frags)* 0.02))
  qc_sele_bc[, 'total_frags' := total_frags -CN]
  qc_sele_bc <- qc_sele_bc[total_frags >= 0]
  
  # fit two NB mixture model & using signal to noisy ratio to select cells
  n_in_peak <- Matrix::colSums(mat)
  n_in_peak <- n_in_peak[names(n_in_peak) %in% qc_sele_bc$bc]
  
  # Need countreg for FLXMRnegbin
  if (!requireNamespace("countreg", quietly = TRUE)) 
      utils::install.packages("countreg", repos="http://R-Forge.R-project.org")

  fm0 <- flexmix::flexmix(n_in_peak ~ 1, k = 2, model = countreg::FLXMRnegbin())
  prob1 <- flexmix::posterior(fm0)[, 1]
  prob2 <- flexmix::posterior(fm0)[, 2]
  mus <- flexmix::parameters(fm0)[1, ]
  
  if(mus[1] > mus[2]){
    odd = prob1
  } else {
    odd = prob2
  }

  aa = which(odd == 1)

  select.cells = names(n_in_peak)[aa]

  return(select.cells)
}


#' @name sc_atac_filter_cell_calling
#' @title filter cell calling
#' @description specify various qc cutoffs to select the desired cells
#' 
#' @param mtx The input matrix
#' @param qc_per_bc_file A file containing qc statistics for each cell
#' @param min_uniq_frags The minimum number of required unique fragments required for a cell
#' @param max_uniq_frags The maximum number of required unique fragments required for a cell
#' @param min_frac_peak The minimum proportion of fragments in a cell to overlap with a peak
#' @param min_frac_tss The minimum proportion of fragments in a cell to overlap with a tss
#' @param min_frac_enhancer The minimum proportion of fragments in a cell to overlap with a enhancer sequence
#' @param min_frac_promoter The minimum proportion of fragments in a cell to overlap with a promoter sequence
#' @param max_frac_mito The maximum proportion of fragments in a cell that are mitochondrial
#' 
#' 
#' @export
#' 
sc_atac_filter_cell_calling <- function(
  mtx, 
  qc_per_bc_file,
  min_uniq_frags = 0,
  max_uniq_frags = 50000,
  min_frac_peak = 0.05,
  min_frac_tss = 0,
  min_frac_enhancer = 0,
  min_frac_promoter = 0,
  max_frac_mito = 0.2) {
  
  # https://github.com/wbaopaul/scATAC-pro/blob/master/scripts/src/filter_barcodes.R

  total_frags <- frac_mito <- frac_peak <- frac_tss <- frac_promoter <- frac_enhancer <- NULL
    
  qc_bc_stat <- data.table::fread(qc_per_bc_file)

  qc_sele <- qc_bc_stat[total_frags >= min_uniq_frags & total_frags <= max_uniq_frags &
                        frac_mito <= max_frac_mito &
                        frac_peak >= min_frac_peak &
                        frac_tss >= min_frac_tss &
                        frac_promoter >= min_frac_promoter &
                        frac_enhancer >= min_frac_enhancer]

  return(qc_sele$bc)
}


#' @name sc_atac_emptydrops_cell_calling
#' @title empty drops cell calling
#' @description The empty drops cell calling method
#' 
#' @param mat The input matrix
#' @param output_folder The path of the output folder
#' @param lower The lower threshold for the data if using the \code{emptydrops} function for cell calling.
#' 
#' @export
#' 
sc_atac_emptydrops_cell_calling <- function(
  mat, 
  output_folder,
  lower = NULL) {
  
  set.seed(2019)
  
  # generating the knee plot
  my.counts <- Matrix::Matrix(mat)
  br.out    <- DropletUtils::barcodeRanks(my.counts)
  
  # Making a plot
  while (!is.null(grDevices::dev.list()))  grDevices::dev.off()
  grDevices::png(file=paste0(output_folder, "/scPipe_atac_stats/knee_plot.png"))
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  graphics::lines(br.out$rank[o], br.out$fitted[o], col="red")
  
  graphics::abline(h=S4Vectors::metadata(br.out)$knee, col="dodgerblue", lty=2)
  graphics::abline(h=S4Vectors::metadata(br.out)$inflection, col="forestgreen", lty=2)
  graphics::legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  grDevices::dev.off()
  
  if(is.null(lower)){
    #lower <- floor(0.1*ncol(mat))
    lower <- 1
  }
  cell.out <- DropletUtils::emptyDrops(mat, lower = lower)
  # cell.out <- emptyDrops2(mat, lower = lower)
  
  filter.out <- cell.out[S4Vectors::complete.cases(cell.out), ]
  
  #saveRDS(filter.out, file = paste0(output_folder, '/EmptyDrop_obj.rds'))
  #cat("Empty cases are removed and saved in ", output_folder, "\n")
  
  if(length(filter.out$FDR) > 0){
    fdr <- 0.01
    cat("FDR of 0.01 is assigned... \n")
    filter.out <- filter.out[filter.out$FDR <= fdr, ]
  } else{
    message("insufficient unique points for computing knee/inflection points ... Use the output matrices with caution! \n")
  }
  
  select.cells <- rownames(filter.out)
  
  return(select.cells)
}







