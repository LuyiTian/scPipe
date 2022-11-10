
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
  
  saveRDS(binary.mat, file = file.path(output_folder, "tf-idf_mat.rds"))
}

#' @name sc_interactive_umap_plot
#' @title Produces an interactive UMAP plot via Shiny
#' @description Can colour the UMAP by any of the colData columns in the SCE object
#' @param sce The SingleCellExperiment object
#' @export
sc_interactive_umap_plot <- function(sce) {
  umap_data <- sc_get_umap_data(sce)
  
  library(shiny)
  
  ui <- fluidPage(
    titlePanel("Interactive UMAP plot"),
    
    sidebarLayout(
      
      sidebarPanel( # INPUT
        
        selectInput("by", "Colour by",
                    colnames(umap_data)[!colnames(umap_data) %in% c("barcode", "UMAP1", "UMAP2")]),
      ),
      mainPanel( # OUTPUT      
        plotlyOutput(outputId = "umapPlot")
      )
    )
  )
  
  # Define server logic required to draw a histogram ----
  server <- function(input, output) {
    output$umapPlot <- renderPlotly({
      g <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, text = paste("barcode: ", barcode))) +
        geom_point(aes(col = .data[[input$by]]), size = 2, alpha = 0.5) +
        theme_bw(base_size = 14)
      if (is.numeric(umap_data[[input$by]])) {
        g <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, text = paste("barcode: ", barcode))) +
          geom_point(aes(col = .data[[input$by]]), size = 2, alpha = 0.5) +
          scale_colour_gradientn(colours = c("green", "red")) +
          theme_bw(base_size = 14)
      }
      plotly::ggplotly(g)
    })
  }
  
  shinyApp(ui = ui, server = server)
}


#' @name sc_get_umap_data
#' @title Generates UMAP data from sce object
#' @description Produces a DataFrame containing the UMAP dimensions, as well as all the colData of the sce object for each cell
#' @param mae The SingleCellExperiment object
#' @param n_neighbours No. of neighbours for KNN 
#' @export
sc_get_umap_data <- function(sce,
                             n_neighbours = 30) {
  set.seed(123)
  
  counts <- assay(sce)
  bin_mat <- as.matrix((counts>0)+0)
  
  TF.IDF.custom <- function(binary.mat, verbose = TRUE) {
    object <- binary.mat
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
  
  binary.mat <- TF.IDF.custom(bin_mat)
  n_bcs <- max(min(50, ncol(binary.mat), nrow(binary.mat))-1,0)
  mat.lsi          <- irlba(binary.mat, n_bcs)
  d_diagtsne       <- matrix(0, n_bcs, n_bcs)
  diag(d_diagtsne) <- mat.lsi$d
  mat_pcs          <- t(d_diagtsne %*% t(mat.lsi$v))
  rownames(mat_pcs)<- colnames(binary.mat)
  
  # clustering in the PCA space using KNN --------------
  knn.info<- RANN::nn2(mat_pcs, k = n_neighbours)
  
  ## convert to adjacency matrix
  knn           <- knn.info$nn.idx
  adj           <- matrix(0, nrow(mat_pcs), nrow(mat_pcs))
  rownames(adj) <- colnames(adj) <- rownames(mat_pcs)
  for(i in seq_len(nrow(mat_pcs))) {
    adj[i,rownames(mat_pcs)[knn[i,]]] <- 1
  }
  
  ## convert to graph
  g <- igraph::graph.adjacency(adj, mode="undirected")
  g <- simplify(g) ## remove self loops
  
  # identify communities, many algorithms. Use the Louvain clustering ------------
  km         <- igraph::cluster_louvain(g)
  com        <- km$membership
  names(com) <- km$names
  
  # running UMAP ------------------------------
  norm.data.umap    <- umap::umap(mat_pcs)
  
  df_umap           <- as.data.frame(norm.data.umap$layout)
  colnames(df_umap) <- c("UMAP1", "UMAP2")
  df_umap$barcode   <- rownames(mat_pcs)
  
  df_umap           <- dplyr::left_join(df_umap, enframe(com), by = c("barcode" = "name")) %>%
    dplyr::rename(cluster = value) %>%
    dplyr::mutate(cluster = as.factor(cluster))
  
  
  sce_coldata <- colData(sce)
  df_umap <- base::merge(df_umap, sce_coldata, by.x = "barcode", by.y = "row.names", all.x = TRUE) 
  
  
  return(df_umap)
}
