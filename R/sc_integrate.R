#####################################################
# Integrate multi-omic scRNA-Seq and scATAC-Seq data
#####################################################

#' @name sc_integrate
#' @title Integrate multi-omic scRNA-Seq and scATAC-Seq data into a MultiAssayExperiment
#' @description Generates an integrated SCE object with scRNA-Seq and scATAC-Seq data produced by the scPipe pipelines
#' @param sce_list A list of SCE objects, named with the corresponding technologies
#' @param sce_column_to_barcode_files A list of files containing the barcodes for each tech (if not needed then give a `NULL` entry)
#' @param cell_line_info A list of files, each of which contains 2 columns corresponding to the barcode and cell line for each tech (if not needed then provide a `NULL` entry)
#' @param barcode_match_file A .csv file with columns corresponding to the barcodes for each tech
#' @param output_folder The path to the output folder
#' 
#' @examples
#' \dontrun{
#' sc_integrate(
#'    sce_list = list("RNA" = sce.rna, "ATAC" = sce.atac),
#'    barcode_match_file = bc_match_file,
#'    sce_column_to_barcode_files = list("RNA" = rna_bc_anno, "ATAC" = NULL),
#'    rev_comp = list("RNA" = FALSE, "ATAC" = TRUE),
#'    cell_line_info = list("RNA" = rna_cell_line_info, "ATAC" = atac_cell_line_info,)
#'    output_folder = output_folder
#'    )
#' }  
#' 
#' @export
#'
sc_integrate <- function(sce_list,
                         barcode_match_file = NULL,
                         sce_column_to_barcode_files = NULL,
                         rev_comp = NULL,
                         cell_line_info = NULL,
                         output_folder = NULL) {
  
  if(is.null(output_folder)) {
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory is not provided. Created directory: ", output_folder, "\n")
  }
  techs <- names(sce_list)
  
  # Apply reverse complement if specified
  if (!is.null(rev_comp)) {
    cat("Applying reverse complement to columns of SCE objects where required.\n")
    sce_list <- lapply(seq_along(sce_list), function(i) {
      sce <- sce_list[[i]]
      if (isTRUE(rev_comp[[i]]))
        colnames(sce) <- Biostrings::reverseComplement(DNAStringSet(colnames(sce))) %>% as.character()
      sce
    })
  }
  names(sce_list) <- techs # add names back
  
  if (!is.null(sce_column_to_barcode_files)) {
    cat("Updating columns of SCE objects to barcodes where required.\n")
    sce_list <- lapply(seq_along(sce_list), function(i) {
      tech <- names(sce_list)[[i]]
      sce <- sce_list[[i]]
      column_to_bc_file <- sce_column_to_barcode_files[tech][[1]]
      if(!is.null(column_to_bc_file)) {
        column_to_bc <- read.csv(column_to_bc_file, header = TRUE)
        colnames(column_to_bc)[1:2] <- c("cell_name", "barcode_sequence")
        if (!all(colnames(sce) %in% column_to_bc$cell_name)) 
          stop("Columns of SCE object not present in annotation file!")
        colnames(sce) <- column_to_bc$barcode_sequence[match(colnames(sce), column_to_bc$cell_name)]
      }
      sce
    })
    names(sce_list) <- techs # add names back
  }
  
  if (!is.null(barcode_match_file)) {
    barcode_match_df <- read.csv(barcode_match_file, header = TRUE)
    # Check if the barcodes match up with those in the match file, and if not, try taking the reverse complement to see if it helps
    sce_list <- lapply(seq_along(sce_list), function(i) {
      sce <- sce_list[[i]]
      tech <- names(sce_list)[[i]]
      matched <- na.omit(match(colnames(sce), barcode_match_df[[tech]]))
      if (length(matched)/length(colnames(sce)) < 0.2) {
        cat("Less than 20% of", tech, "barcodes match with the barcodes in the barcode match file. Trying the reverse complement.\n")
        rev_comp <- Biostrings::reverseComplement(DNAStringSet(colnames(sce))) %>% as.character()
        new_matched <- na.omit(match(rev_comp, barcode_match_df[[tech]]))
        if (length(new_matched)/length(colnames(sce)) <= length(matched)/length(colnames(sce))) {
          cat("Reverse complement didn't yield higher proportion of matched barcodes so using existing barcodes.\n")
        } else {
          cat("Reverse complement yielded higher proportion of matched barcodes so using reverse complement\n")
          colnames(sce) <- rev_comp
        }
      }
      sce
    })
    names(sce_list) <- techs # add names back
  }
  
  
  # Left join the cell line info
  if (!is.null(cell_line_info)) {
    cat("Attaching cell line information provided.\n")
    sce_list <- lapply(seq_along(sce_list), function(i) {
      sce <- sce_list[[i]]
      tech <- names(sce_list)[[i]]
      cell_line_file <- cell_line_info[[i]]
      
      if (!is.null(cell_line_file)) {
        cli <- read.csv(cell_line_file, header=FALSE)
        colnames(cli) <- c("barcode", "cell_line")
        merged_cell_line_info <- base::merge(data.frame(barcode = colnames(sce)), cli, all.x = TRUE) %>% column_to_rownames("barcode")
        reordered_cli <- merged_cell_line_info[colnames(sce), , drop=FALSE]
        colData(sce) <- cbind(colData(sce), cell_line = reordered_cli$cell_line)
      } 
      sce
    })
    names(sce_list) <- techs # add names back
  }
  
  # Outer join the column data
  if (!is.null(barcode_match_file)) { 
    cat("Merging qc metrics\n")
    qc_dfs <- lapply(seq_along(sce_list), function(i) {
      sce <- sce_list[[i]]
      tech <- techs[[i]]
      qc <- data.frame(colData(sce)) %>% rename_with( ~ paste0(tech, "_", .x))
      base::merge(qc, barcode_match_df, by.x = "row.names", by.y = techs[[i]], all.x = TRUE) %>% rename(!!tech := "Row.names")
    })
    merged_qc <- Reduce(function(x, y) base::merge(x, y, all = TRUE), qc_dfs)
  }
  
  # Mark in coldata of each experiment if the barcode is shared
  shared_bc <- na.omit(merged_qc)
  for (i in 1:length(sce_list)) {
    tech <- techs[[i]]
    colData(sce_list[[tech]])[, "shared"] <- rownames(colData(sce_list[[tech]])) %in% shared_bc[[tech]]
  }
  
  # Create MultiAssayExperiment
  cat("Creating MultiAssayExperiment\n")
  mae <- MultiAssayExperiment(experiments = sce_list)
  mae@metadata$scPipe$version <- packageVersion("scPipe") 
  if (!is.null(barcode_match_file)) mae@metadata$scPipe$integrated_qc <- merged_qc
  
  # Save the MAE object
  # saveRDS(mae, file.path(output_folder, "scPipe_MAE_object.rds"))
  
  cat("sc_integrate_complete.\n")
  return(mae)
}

#' @name sc_mae_plot_umap
#' @title Generates UMAP of multiomic data
#' @description Uses feature count data from multiple experiment objects to produce a UMAP 
#' @param mae The MultiAssayExperiment object
#' @param by What to colour the points by. Needs to be in colData of all experiments.
#' @param output_file The path of the output file
#' @export
sc_mae_plot_umap <- function(mae,
                             by = NULL,
                             output_file = NULL) {
  set.seed(123)
  
  umap_dfs <- lapply(seq_along(mae), function(i) {
    tech <- names(mae)[[i]]
    counts <- assays(mae)[[tech]]
    bin_mat <- as.matrix((counts>0)+0)
    binary.mat <- TF.IDF.custom(bin_mat)
    n_bcs <- max(min(50, ncol(binary.mat), nrow(binary.mat))-1,0)
    mat.lsi          <- irlba(binary.mat, n_bcs)
    d_diagtsne       <- matrix(0, n_bcs, n_bcs)
    diag(d_diagtsne) <- mat.lsi$d
    mat_pcs          <- t(d_diagtsne %*% t(mat.lsi$v))
    rownames(mat_pcs)<- colnames(binary.mat)
    
    # clustering in the PCA space using KNN --------------
    knn.info<- RANN::nn2(mat_pcs, k = 30)
    
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
    
    df_umap$source <- tech
    
    if (!is.null(by) && !is.null(colData(experiments(mae)[[tech]])[[by]])) {
      sce_coldata <- colData(experiments(mae)[[tech]])[, c(by), drop=FALSE]
      df_umap <- base::merge(df_umap, sce_coldata, by.x = "barcode", by.y = "row.names", all.x = TRUE) 
    } else if (!is.null(by)) {
      df_umap[[by]] <- NA
    }
    df_umap
  })
  names(umap_dfs) <- names(mae)
  
  umap_data <- do.call(rbind, umap_dfs)
  if (!is.null(by)) {
    g <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(col = .data[[by]]), size = 0.5) +
      theme_bw(base_size = 14)
    if (is.numeric(umap_data[[by]])) {
      g <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(col = .data[[by]]), size = 0.5) +
        scale_colour_gradientn(colours=c("green","black")) +
        theme_bw(base_size = 14)
    }
  } else {
    g <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(col = source), size = 0.5) +
      theme_bw(base_size = 14)
  }
  
  if (!is.null(output_file)) {
    if (file.exists(output_file)) {
      ggsave(output_file)
      cat("Saved plot to", output_file, "\n")
    } else {
      cat("The supplied output file path was invalid.\n")
    }
  }
  
  # plotly::ggplotly(g)
  
  return(g)
}