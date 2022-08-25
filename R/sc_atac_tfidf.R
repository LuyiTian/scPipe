
###########################################################
# TF-IDF normalisation and UMAP for scATAC-Seq data
###########################################################

#' @name sc_atac_tfidf
#'
#' @title generating the UMAPs for sc-ATAC-Seq preprocessed data
#' @description Takes the binary matrix and generate a TF-IDF so the clutering can take place on the reduced dimentions. 
#' 
#' @param binary.mat The final, filtered feature matrix in binary format 
#' @param output_folder The path of the output folder 
#'
#' @examples
#' \dontrun{
#' sc_atac_tfidf(binary.mat = final_binary_matrix) 
#' }
#' 
#' @export
#' 

sc_atac_tfidf <- function(binary.mat, output_folder = NULL) {

  # Check if output directory exists
  if (is.null(output_folder)) {
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory does not exist. Created directory: ", output_folder, "\n")
  }
  
  #Create log folder/file
  log_and_stats_folder <- file.path(output_folder, "scPipe_atac_stats")
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file            <- file.path(log_and_stats_folder, "log_file.txt")
  if(!file.exists(log_file)) file.create(log_file)
  cat(
    paste0(
      "sc_atac_tfidf starts at ",
      as.character(Sys.time()),
      "\n"
    ),
    file = log_file, append = TRUE)
  
    
  TF.IDF.custom <- function(binary.mat, verbose = TRUE) {
    if (class(x = binary.mat) == "data.frame") {
      binary.mat <- as.matrix(x = binary.mat)
    }
    if (class(x = binary.mat) != "dgCMatrix") {
      binary.mat <- as(object = binary.mat, Class = "dgCMatrix")
    }
    if (verbose) {
      message("Performing TF-IDF normalization")
    }
    
    npeaks       <- Matrix::colSums(x = object)
    tf           <- Matrix::tcrossprod(x = as.matrix(object), y = Matrix::Diagonal(x = 1 / npeaks))
    rsums        <- Matrix::rowSums(x = object)
    idf          <- ncol(x = object) / rsums
    norm.data    <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
    scale.factor <- 1e4
    slot(object = norm.data, name = "x") <- log1p(x = slot(object = norm.data, name = "x") * scale.factor)
    norm.data[which(x = is.na(x = norm.data))] <- 0
    return(norm.data)
  }
  
  
  binary.mat <- TF.IDF.custom(binary.mat)
  
  library(irlba)
  set.seed(123)
  mat.lsi          <- irlba(binary.mat, 50)
  d_diagtsne       <- matrix(0, 50, 50)
  diag(d_diagtsne) <- mat.lsi$d
  mat_pcs          <- t(d_diagtsne %*% t(mat.lsi$v))
  rownames(mat_pcs)<- colnames(binary.mat)
  
  # clustering in the PCA space using KNN --------------
  
  library(RANN)
  knn.info<- RANN::nn2(mat_pcs, k = 30)
  
  ## convert to adjacency matrix
  knn           <- knn.info$nn.idx
  adj           <- matrix(0, nrow(mat_pcs), nrow(mat_pcs))
  rownames(adj) <- colnames(adj) <- rownames(mat_pcs)
  for(i in seq_len(nrow(mat_pcs))) {
    adj[i,rownames(mat_pcs)[knn[i,]]] <- 1
  }
  
  ## convert to graph
  library(igraph)
  g <- igraph::graph.adjacency(adj, mode="undirected")
  g <- simplify(g) ## remove self loops
  
  # identify communities, many algorithums. Use the Louvain clustering ------------
  km         <- igraph::cluster_louvain(g)
  com        <- km$membership
  names(com) <- km$names
  
  # running UMAP ------------------------------
  
  library(umap)
  library(ggplot2)
  library(tibble)
  set.seed(345)
  
  norm.data.umap    <- umap::umap(mat_pcs)
  
  df_umap           <- as.data.frame(norm.data.umap$layout)
  colnames(df_umap) <- c("UMAP1", "UMAP2")
  df_umap$barcode   <- rownames(mat_pcs)
  
  df_umap           <- dplyr::left_join(df_umap, enframe(com), by = c("barcode" = "name")) %>%
    dplyr::rename(cluster = value) %>%
    dplyr::mutate(cluster = as.factor(cluster))
  
  
  ggplot(df_umap, aes(x = UMAP1, y = UMAP2)) + 
    geom_point(aes(col = cluster), size = 0.5) +
    theme_bw(base_size = 14)
   
}
