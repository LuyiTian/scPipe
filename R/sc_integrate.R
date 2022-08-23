#####################################################
# Integrate multi-omic scRNA-Seq and scATAC-Seq data
#####################################################

#' @name sc_integrate
#' @title Integrate multi-omic scRNA-Seq and scATAC-Seq data into a MultiAssayExperiment
#' @description Generates an integrated SCE object with scRNA-Seq and scATAC-Seq data produced by the scPipe pipelines
#' @param sce_list A list of SCE objects, named with the corresponding technologies
#' @param sce_column_to_barcode_files A list of files containing the barcodes for each tech (if not needed then give a `NULL` entry)
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
#'    output_folder = output_folder
#'    )
#' }  
#' 
#' @export
#'
sc_integrate <- function(sce_list,
                         barcode_match_file,
                         sce_column_to_barcode_files = NULL,
                         rev_comp = NULL,
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
  
  # Check if the barcodes match up with those in the match file, and if not, try taking the reverse complement to see if it helps
  sce_list <- lapply(seq_along(sce_list), function(i) {
    sce <- sce_list[[i]]
    tech <- names(sce_list)[[i]]
    matched <- na.omit(match(colnames(sce), barcode_match_df[[tech]]))
    if (length(matched)/length(colnames(sce)) < 0.2) {
      # Try taking reverse comp
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
  
  # Outer join the column data
  cat("Merging qc metrics\n")
  barcode_match_df <- read.csv(barcode_match_file, header = TRUE)
  first_tech <- names(sce_list)[[1]]
  
  qc_dfs <- lapply(seq_along(sce_list), function(i) {
    sce <- sce_list[[i]]
    tech <- techs[[i]]
    qc <- data.frame(colData(sce)) %>% rename_with( ~ paste0(tech, "_", .x))
    base::merge(qc, barcode_match_df, by.x = "row.names", by.y = techs[[i]], all.x = TRUE) %>% rename(!!tech := "Row.names")
  })
  merged_qc <- Reduce(function(x, y) base::merge(x, y, all = TRUE), qc_dfs)
  
  # Create MultiAssayExperiment
  cat("Creating MultiAssayExperiment\n")
  mae <- MultiAssayExperiment(experiments = sce_list)
  mae@metadata$scPipe$version <- packageVersion("scPipe") 
  mae@metadata$scPipe$integrated_qc <- merged_qc
  
  # Save the MAE object
  saveRDS(mae, file.path(output_folder, "scPipe_MAE_object.rds"))
  
  cat("sc_integrate_complete.\n")
  return(mae)
}
