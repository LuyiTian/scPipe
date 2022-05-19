
#########################################################################
# Generating the Feature Matrix from a Given BAM file and a Feature FIle
#########################################################################

#' @name sc_atac_feature_counting
#' @title generating the feature by cell matrix 
#' @description feature matrix is created using a given demultiplexed BAM file and 
#' a selected feature type
#' @param insortedbam The input bam
#' @param feature_input The feature input data e.g. the .narrowPeak file for peaks of a bed file format
#' @param bam_tags The BAM tags
#' @param feature_type The type of feature
#' @param organism The organism type (contains hg19, hg38, mm10)
#' @param cell_calling The desired cell calling method; either \code{cellranger}, \code{emptydrops} or  \code{filter}.
#' @param pheno_data The phenotypic data as a data frame
#' 
#' @param promoters_file The path of the promoter annotation file (if the specified organism isn't recognised).
#' @param tss_file The path of the tss annotation file (if the specified organism isn't recognised).
#' @param enhs_file The path of the enhs annotation file (if the specified organism isn't recognised).
#' @param gene_anno_file The path of the gene annotation file (gtf or gff3 format).
#' @param sample_name The sample name to identify which is the data is analysed for.
#' 
#' @param bin_size The size of the bins
#' @param yieldsize The yield size
#' @param mapq The minimum MAPQ score
#' @param n_filter_cell_counts An integer value to filter the feature matrix on the number of reads per cell (default = 200)
#' @param n_filter_feature_counts An integer value to filter the feature matrix on the number of reads per feature (default = 10).
#' @param exclude_regions Whether or not the regions (specified in the file) should be excluded
#' @param excluded_regions_filename The filename of the file containing the regions to be excluded
#' @param fix_chr Whether chr should be fixed or not
#' 
#' @param lower the lower threshold for the data if using the \code{emptydrops} function for cell calling
#' @param genome_size The size of the genome (used for the \code{cellranger} cell calling method)
#'
#' @param min_uniq_frags The minimum number of required unique fragments required for a cell (used for \code{filter} cell calling)
#' @param max_uniq_frags The maximum number of required unique fragments required for a cell (used for \code{filter} cell calling)
#' @param min_frac_peak The minimum proportion of fragments in a cell to overlap with a peak (used for \code{filter} cell calling)
#' @param min_frac_tss The minimum proportion of fragments in a cell to overlap with a tss (used for \code{filter} cell calling)
#' @param min_frac_enhancer The minimum proportion of fragments in a cell to overlap with a enhancer sequence (used for \code{filter} cell calling)
#' @param min_frac_promoter The minimum proportion of fragments in a cell to overlap with a promoter sequence (used for \code{filter} cell calling)
#' @param max_frac_mito The maximum proportion of fragments in a cell that are mitochondrial (used for \code{filter} cell calling)
#' @param output_folder The output folder
#' @param create_report Logical value to say whether to create the report or not (default = TRUE).

#' @importFrom BiocGenerics start end which strand start<- end<- as.data.frame
#' @importFrom rlang .data
#' @import dplyr
#' @import tidyr
#' 
#' @examples
#' \dontrun{
#' sc_atac_feature_counting(
#'    insortedbam = tiny_tagged_sorted_bam,
#'    cell_calling = "filter",
#'    exclude_regions = TRUE,
#'    feature_input = feature_file)
#' }    
#' @export
#' 

sc_atac_feature_counting <- function(
  insortedbam, 
  feature_input             = NULL, 
  bam_tags                  = list(bc="CB", mb="OX"), 
  feature_type              = "peak", 
  organism                  = "hg38", 
  cell_calling              = "filter", # either c("cellranger", "emptydrops", "filter")
  sample_name               = "",
  genome_size               = NULL, # this is optional but needed if the cell_calling option is cellranger AND organism in NULL
  promoters_file            = NULL,
  tss_file                  = NULL,
  enhs_file                 = NULL,
  gene_anno_file            = NULL,
  pheno_data                = NULL, 
  bin_size                  = NULL, 
  yieldsize                 = 1000000,
  mapq                      = 30,
  n_filter_cell_counts      = 200,
  n_filter_feature_counts   = 10,
  exclude_regions           = FALSE, 
  excluded_regions_filename = NULL,
  output_folder             = NULL,
  fix_chr                   = "none", # should be either one of these: c("none", "excluded_regions", "feature", "both")
  lower                     = NULL,
  min_uniq_frags            = 3000,
  max_uniq_frags            = 50000,
  min_frac_peak             = 0.3,
  min_frac_tss              = 0,
  min_frac_enhancer         = 0,
  min_frac_promoter         = 0.1,
  max_frac_mito             = 0.15,
  create_report             = TRUE
) {
  
  . <- V1 <- V2 <- V3 <- init <- peakStart <- seqnames <- peakEnd <- CB <- chromosome <- chrS <- feature <- barcode <- NULL
  
  init_time = Sys.time()
  
  available_organisms = c("hg19",
                          "hg38",
                          "mm10")
  
  if(!is.null(organism)) stopifnot(organism %in% available_organisms)
  
  stopifnot(fix_chr %in% c("none", "excluded_regions", "feature", "both"))
  
  if(is.null(output_folder)) {
    output_folder <- file.path(getwd(), "scPipe-atac-output")
    cat("Output directory has not been provided. Saving output in\n", output_folder, "\n")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory does not exist. Creating...", output_folder, "\n")
  }
  
  
  # initiate log file and location for stats in scPipe_atac_stats
  log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file             <- paste0(log_and_stats_folder, "log_file.txt")
  stats_file           <- paste0(log_and_stats_folder, "stats_file_align.txt")
  if(!file.exists(log_file)) file.create(log_file)
  # file.create(stats_file)
  
  # timer
  cat(
    paste0(
      "sc_atac_feature_counting starts at ",
      as.character(Sys.time()),
      "\n"
    ), 
    file = log_file, append = TRUE)

  
  
  ############# feature type is genome_bin ####################
  
  if(feature_type == 'genome_bin'){
    # TODO: test the format of the fasta file here and stop if not proper format
    cat("`genome bin` feature type is selected for feature input. reading the genome fasta file ...", "\n")
    
    # index for feature_input is created if not found in the same directory as the feature_input
    if(!file.exists(paste0(feature_input,".fai"))){
      Rsamtools::indexFa(feature_input)
      cat("Index for ", feature_input, " is being created... ", "\n")
    }
    
    if(is.null(bin_size)){
      bin_size <- 2000
      cat("Default bin size of 2000 is selected", "\n")
    }
    
    out_bed_filename <- paste0(output_folder, "/", sub('\\..[^\\.]*$', '', basename(feature_input)), ".bed") # remove extension and append output folder
    
    if(file.exists(out_bed_filename)) {
      message(out_bed_filename, " file already exists. Replacing!")
      file.remove(out_bed_filename)
    }
    
    if(!is.null(feature_input)){
      
      rcpp_fasta_bin_bed_file(feature_input, out_bed_filename, bin_size)
      
      # Check if file was created
      if(file.exists(out_bed_filename)) {
        cat("Generated the genome bin file:\n", out_bed_filename, "\n")
      }
      else {
        stop("File ", out_bed_filename, "file was not created.")
      }
      
      ## End if is.null
    } else{
      if(!is.null(organism)){
        
        # Use organism data sizes saved in repository
        organism_files     <- list.files(system.file("extdata/annotations/", package = "scPipe", mustWork = TRUE))
        sizes_filename_aux <- grep(pattern = organism, x = organism_files, value = TRUE) %>% 
          grep(pattern     <- "size", x = ., value = TRUE)
        
        sizes_filename <- system.file(paste0("extdata/annotations/", sizes_filename_aux), package = "scPipe", mustWork = TRUE)
        
        sizes_df <- utils::read.table(sizes_filename, header = FALSE, col.names = c("V1", "V2"))
        
        sizes_df_aux <- sizes_df %>% 
          dplyr::group_by(V1) %>% 
          dplyr::mutate(V3 = floor(V2/bin_size),
                        V4 = bin_size*V3) %>% 
          dplyr::ungroup()
        
        
        
        out_bed <- purrr::map_df(1:nrow(sizes_df_aux), function(i){
          
          aux_i        <- sizes_df_aux %>% dplyr::slice(i)
          seq_aux_end  <- seq(from = bin_size, to = aux_i$V4, by = bin_size)
          seq_aux_init <- seq_aux_end - bin_size + 1
          
          if(seq_aux_end[length(seq_aux_end)] == aux_i$V2){
            last_row_init <- c()
            last_row_end  <- c()
          } else{
            last_row_init <- seq_aux_end[length(seq_aux_end)] + 1
            last_row_end  <- aux_i$V2
          }
          
          options(scipen = 999) # To prevent scientific notation
          out_df <- dplyr::tibble(init = c(seq_aux_init, last_row_init), 
                                  end = c(seq_aux_end, last_row_end)
          ) %>% 
            dplyr::mutate(name = aux_i$V1, .before = init)
          
          return(out_df)
        })
        
        utils::write.table(out_bed, file = out_bed_filename, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
        
        # Check if file was created
        if(file.exists(out_bed_filename)) {
          cat("Generated the genome bin file:\n", out_bed_filename, "\n")
        }
        else {
          stop("File ", out_bed_filename, "file was not created.")
        }
      } else{
        stop("Either organism or feature_input should be provided.")
      }
    } # End else is.null(genome_size)
    
    # feature_input <- genome_bin
  }
  
  ############################## feature type is peak #####################
  if(feature_type == 'peak' || feature_type == 'tss' || feature_type == 'gene'){
    cat("`peak`, `tss` or `gene` feature_type is selected for feature input", "\n")
    
    ############################### fix_chr in feature file
    if(fix_chr %in% c("feature", "both")){
      
      out_bed_filename_feature <- paste0(output_folder, "/", 
                                         sub('\\..[^\\.]*$', 
                                             '', 
                                             basename(feature_input)
                                         ), "_fixedchr.bed"
      ) # remove extension and append output folder
      
      # Try to read first 5 rows of feature_input file to see if the format is correct
      feature_head <- utils::read.table(feature_input, nrows = 5)
      if(ncol(feature_head) < 3){
        warning("Feature file provided does not contain 3 columns. Cannot append chr")
        break;
      }
      
      
      rcpp_append_chr_to_bed_file(feature_input, out_bed_filename_feature)
      
      # Check if file was created
      if(file.exists(out_bed_filename_feature)) {
        cat("Appended 'chr' to feature file and output created in:", out_bed_filename_feature, "\n")
        feature_input <- out_bed_filename_feature
      }
      else {
        stop("File ", out_bed_filename_feature, "file was not created.")
      }
      
    }
    
    
    
    ############################### fix_chr in excluded regions
    if(fix_chr %in% c("excluded_regions", "both")){
      if(!is.null(excluded_regions_filename)){
        
        out_bed_filename_excluded_regions <- paste0(output_folder, "/", 
                                                    sub('\\..[^\\.]*$', 
                                                        '', 
                                                        basename(excluded_regions_filename)
                                                    ), "_fixedchr.bed"
        ) # remove extension and append output folder
        
        # Try to read first 5 rows of excluded_regions file to see if the format is correct
        excluded_regions_head <- utils::read.table(excluded_regions_filename, nrows = 5)
        if(ncol(excluded_regions_head) < 3){
          warning("excluded_regions file provided does not contain 3 columns. Cannot append chr")
          break
        }
        
        rcpp_append_chr_to_bed_file(excluded_regions_filename, out_bed_filename_excluded_regions)
        
        # Check if file was created
        if(file.exists(out_bed_filename_excluded_regions)) {
          cat("Appended 'chr' to excluded_regions file and output created in:", out_bed_filename_excluded_regions, "\n")
          excluded_regions_filename <- out_bed_filename_excluded_regions
        }
        else {
          stop("File ", out_bed_filename_excluded_regions, "file was not created.")
        }
      }
    
    } # end if excluded_regions
  }
  
  ############################### end fix_chr
  
  ###################### generate the GAlignment objects from BAM and features ######################
  

  ############## generate the GAalignment file from the feature_input file
  cat("Creating Galignment object for the feature input...\n")
  if(feature_type != 'genome_bin'){
    feature.gr <- rtracklayer::import(feature_input)
  } else {
    feature.gr <- rtracklayer::import(out_bed_filename)
  }
  
  ############################### exclude regions 
  
  number_of_lines_to_remove <- 0
  
  if(exclude_regions){
    
    if(is.null(excluded_regions_filename) & !is.null(organism)){
      ## If excluded_regions_filename is null but organism is not, then read the file from system
      
      organism_files                <- list.files(system.file("extdata/annotations/", package = "scPipe", mustWork = TRUE))
      excluded_regions_filename_aux <- grep(pattern = organism, x = organism_files, value = TRUE) %>% 
        grep(pattern = "blacklist", x = ., value = TRUE)
      excluded_regions_filename     <- system.file(paste0("extdata/annotations/", excluded_regions_filename_aux), package = "scPipe", mustWork = TRUE)
    } 
    
    if(!is.null(excluded_regions_filename)){
      excluded_regions.gr               <- rtracklayer::import(excluded_regions_filename)
      overlaps_excluded_regions_feature <- IRanges::findOverlaps(excluded_regions.gr, feature.gr, maxgap = -1L, minoverlap = 0L) # Find overlaps
      lines_to_remove                   <- as.data.frame(overlaps_excluded_regions_feature)$subjectHits # Lines to remove in feature file
      number_of_lines_to_remove         <- length(lines_to_remove)
      
      if(number_of_lines_to_remove > 0){ # If there are lines to remove
        feature.gr.df            <- as.data.frame(feature.gr)
        lines_to_keep            <- setdiff(1:nrow(feature.gr.df), lines_to_remove)
        feature.gr               <- feature.gr[lines_to_keep, ]
      } # End if(number_of_lines_to_remove > 0)
      
    } else{
      warning("Parameter exclude_regions was TRUE but no known organism or excluded_regions_filename provided. Proceding without excluding regions.")
    }
    
    
  } # End if(exclude_regions)
  
  # default values of findoverlaps (type= "any") works and we don't need to compute below paramters
  #median_feature_overlap <- stats::median(GenomicAlignments::ranges(feature.gr)@width)
  #minoverlap             <- 0.1*median_feature_overlap
  #maxgap                 <- 0.1*median_feature_overlap
  

  
  ############## read in the aligned and demultiplexed BAM file
  
  
  add_matrices <- function(...) {
    a <- list(...)
    cols <- sort(unique(unlist(lapply(a, colnames))))
    rows <- sort(unique(unlist(lapply(a, rownames))))

    nrows <- length(rows)
    ncols <- length(cols)
    newms <- lapply(a, function(m) {
      b <- as(as(m, "dgCMatrix"), "dgTMatrix")
      s <- cbind.data.frame(i = b@i + 1, j = b@j + 1, x = b@x)
     

      i <- match(rownames(m), rows)[s$i]
      j <- match(colnames(m), cols)[s$j]
     
      Matrix::sparseMatrix(i=i,
                   j=j,
                   x=s$x,
                   dims=c(nrows, ncols),
                   dimnames=list(rows, cols))
    })
    Reduce(`+`, newms)
  }
  
  
  param <- Rsamtools::ScanBamParam(tag = as.character(bam_tags),  mapqFilter=mapq)
  bamfl <- open(Rsamtools::BamFile(insortedbam, yieldSize=yieldsize))
  
  cat("Iterating over chunks of the BAM file and gnerating the feature by barcode count matrix\n")
  
  feature_matrix <- NULL
  iter <- 1
  total_reads <- 0 # total reads in BAM file in that satisfy the specified parameter
  all_bcs <- c() # all the unique barcodes in the BAM file
  last_time <- NULL

  while(length(yld <- GenomicAlignments::readGAlignments(bamfl, use.names = TRUE, param = param))) {
    yld.gr <- GenomicRanges::makeGRangesFromDataFrame(yld, keep.extra.columns=TRUE) 
    cat("Chunk", iter " completed in ")
    if (iter > 1) {
      cat(Sys.time() - last_time, "seconds\n")
    } else {
      cat("\n")
    }
    last_time <- Sys.time()
    iter <- iter+1
    
    bcs <- yld.gr$CB
    all_bcs <- unique(append(all_bcs, bcs)) 
    total_reads <- total_reads + length(yld.gr)
    
    ############# Adjusting for the Tn5 cut site
    isMinus <- which(strand(yld.gr) == "-")
    isOther <- which(strand(yld.gr) != "-")
    #Forward
    start(yld.gr)[isOther] <- start(yld.gr)[isOther] - 5
    end(yld.gr)[isOther] <- end(yld.gr)[isOther] + 4
    #Reverse
    end(yld.gr)[isMinus] <- end(yld.gr)[isMinus] + 5
    start(yld.gr)[isMinus] <- start(yld.gr)[isMinus] - 4
    
    ############### Overlaps
    #cat ("Finding overlaps between alignments and features\n")
    overlaps               <- GenomicAlignments::findOverlaps(query         = feature.gr, # feature.gr,
                                                              subject       = yld.gr, # yld.gr, #
                                                              type          = "any", 
                                                              ignore.strand = TRUE)
    # generate the matrix using this overlap results above.
    
    GenomicRanges::mcols(yld.gr)[S4Vectors::subjectHits(overlaps), "peakStart"] <- start(GenomicAlignments::ranges(feature.gr)[S4Vectors::queryHits(overlaps)])
    GenomicRanges::mcols(yld.gr)[S4Vectors::subjectHits(overlaps), "peakEnd"]   <- end(GenomicAlignments::ranges(feature.gr)[S4Vectors::queryHits(overlaps)])
    
    #is removing NAs here the right thing to do?
    overlap.df <- data.frame(yld.gr) %>% filter(!is.na(peakStart)) %>% select(seqnames, peakStart, peakEnd, CB)
    overlap.df <- overlap.df %>% 
      group_by(seqnames, peakStart, peakEnd, CB) %>% 
      summarise(count = n()) %>% 
      purrr::set_names(c("chromosome","start","end","barcode","count")) %>% 
      unite("chrS", chromosome:start, sep=":") %>%
      unite("feature", chrS:end, sep="-")
    
    matrixData <- overlap.df %>%
      group_by(feature, barcode) %>% 
      spread(barcode, count)
  
    matrixData           <- as.data.frame(matrixData)
    
    rownames(matrixData) <- matrixData$feature
  
    matrixData.old       <- matrixData
    matrixData           <- matrixData %>%
      dplyr::select(-1) %>%
      data.table::as.data.table() %>%
      Biostrings::as.matrix() %>%
      replace(is.na(.), 0)

    # add dimensions of the matrix
    dimnames(matrixData)  <-  list(matrixData.old[1] %>% rownames(), matrixData.old %>% dplyr::select(-1) %>% colnames())
    if (is.null(feature_matrix)) {
      feature_matrix <- Matrix::Matrix(matrixData)
    } else {
      feature_matrix <- add_matrices(feature_matrix, Matrix::Matrix(matrixData))
    }
  }

  # Calculate average no. of reads per cellular barcode  
  average_number_of_lines_per_CB <- total_reads/length(all_bcs)
  cat("Average no. of reads per barcode: ", average_number_of_lines_per_CB, "\n")

  
  saveRDS(feature_matrix, file = file.path(output_folder, "unfiltered_feature_matrix.rds"))
  cat("Raw feature matrix generated: ", file.path(output_folder, "unfiltered_feature_matrix.rds") , "\n")
  
  
  ################ Initiate log file
  log_and_stats_folder       <- paste0(output_folder, "/scPipe_atac_stats/")
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file                   <- paste0(log_and_stats_folder, "log_file.txt")
  if(!file.exists(log_file)) file.create(log_file)

  cat("Average number of reads per cell barcode:", average_number_of_lines_per_CB, "\n",
      file = log_file, append = TRUE)

  cat("Number of regions removed from feature input due to being invalid:", number_of_lines_to_remove, "\n",
      file = log_file, append = TRUE)

  ########## Check if organism is pre-recognized and if so then use the package's annotation files
  if (organism %in% c("hg19", "hg38", "mm10")) {
    cat(organism, "is a recognized organism. Using annotation files in repository.\n")
    anno_paths <- system.file("extdata/annotations/", package = "scPipe", mustWork = TRUE)
    
    promoters_file <- file.path(anno_paths, paste0(organism, "_promoter.bed.gz"))
    tss_file <- file.path(anno_paths, paste0(organism, "_tss.bed.gz"))
    enhs_file <- file.path(anno_paths, paste0(organism, "_enhancer.bed.gz"))
  }
  else if (!all(file.exists(c(promoters_file, tss_file, enhs_file)))) {
    stop("One of the annotation files could not be located. Please make sure their paths are valid.")
  }

  cat(
    paste0(
      "Raw matrix generated at ",
      as.character(Sys.time()),
      "\n"
    ),
    file = log_file, append = TRUE)

  ########### generate quality control metrics for cells
  sc_atac_create_qc_per_bc_file(inbam      = insortedbam,
                                frags_file = file.path(output_folder, "fragments.bed"),
                                peaks_file = feature_input,
                                promoters_file = promoters_file,
                                tss_file = tss_file,
                                enhs_file = enhs_file,
                                output_folder = output_folder)

  qc_per_bc_file <- file.path(output_folder, "qc_per_bc_file.txt")

  cat(
    paste0(
      "Cell QC metrics generated at ",
      as.character(Sys.time()),
      "\n"
    ),
    file = log_file, append = TRUE)

  # Cell calling
  matrixData <- sc_atac_cell_calling(mat = feature_matrix,
                                     cell_calling = cell_calling,
                                     output_folder = output_folder,
                                     genome_size = genome_size,
                                     qc_per_bc_file = qc_per_bc_file,
                                     lower = lower,
                                     min_uniq_frags = min_uniq_frags,
                                     max_uniq_frags = max_uniq_frags,
                                     min_frac_peak = min_frac_peak,
                                     min_frac_tss = min_frac_tss,
                                     min_frac_enhancer = min_frac_enhancer,
                                     min_frac_promoter = min_frac_promoter,
                                     max_frac_mito = max_frac_mito)

  
  cat(
    paste0(
      "Cell calling completed at ",
      as.character(Sys.time()),
      "\n"
    ),
    file = log_file, append = TRUE)
  
  # filter the matrix based on counts per cell
  if (n_filter_cell_counts > 0) {
    message(paste0("cells with less than ", n_filter_cell_counts, " counts are filtered out."))
    cat(
      paste0(
        "cells with less than ",
        n_filter_cell_counts," counts are filtered out.",
        "\n"
      ),
      file = log_file, append = TRUE)
    
    filtered_indices <- (colSums(matrixData, na.rm=TRUE) > n_filter_cell_counts)
    matrixData        <- matrixData[, filtered_indices] # all the remaining columns
  } else {
    message("no cells were filtered out based on counts.")
  }
  
  
  
  # filter the matrix based on counts per feature
  if (n_filter_feature_counts > 0) {
    message(paste0("features with less than ", n_filter_feature_counts, " counts are filtered out."))
    cat(
      paste0(
        "features with less than ",
        n_filter_feature_counts," counts are filtered out.",
        "\n"
      ),
      file = log_file, append = TRUE)
    
    filtered_indices <- (rowSums(matrixData, na.rm=TRUE) > n_filter_feature_counts)
    matrixData        <- matrixData[filtered_indices,] # all the remaining rows
  } else {
    message("no features were filtered out based on counts.")
  }
  

  # converting the NAs to 0s if the sparse option to create the sparse Matrix properly
  sparseM <- Matrix::Matrix(matrixData, sparse=TRUE)
  # # add dimensions of the sparse matrix if available
  # if(cell_calling != FALSE){
  #   barcodes <- utils::read.table(paste0(output_folder, '/non_empty_barcodes.txt'))$V1
  #   features <- utils::read.table(paste0(output_folder, '/non_empty_features.txt'))$V1
  #   dimnames(sparseM) <- list(features, barcodes)
  # }

  cat("Sparse matrix generated", "\n")
  saveRDS(sparseM, file = paste(output_folder, "/sparse_matrix.rds", sep = ""))
  # Matrix::writeMM(obj = sparseM, file=paste(output_folder, "/sparse_matrix.mtx", sep =""))
  cat("Sparse count matrix is saved in\n", paste(output_folder,"/sparse_matrix.mtx",sep = "") , "\n")

  # generate and save jaccard matrix
  # jaccardM <- locStra::jaccardMatrix(sparseM)
  # cat("Jaccard matrix generated", "\n")
  # saveRDS(jaccardM, file = paste(output_folder,"/jaccard_matrix.rds",sep = ""))
  # cat("Jaccard matrix is saved in\n", paste(output_folder,"/jaccard_matrix.rds",sep = "") , "\n")

  # generate and save the binary matrix
  matrixData[matrixData>0] <- 1
  saveRDS(matrixData, file = paste(output_folder,"/binary_matrix.rds",sep = ""))
  cat("Binary matrix is saved in:\n", paste(output_folder,"/binary_matrix.rds",sep = "") , "\n")
  
  ###### calculate teh TF-IDF matrices here : TO-DO

  # following can be used to plot the stats and load it into sce object ########
  # (from https://broadinstitute.github.io/2020_scWorkshop/data-wrangling-scrnaseq.html)
  counts_per_cell    <- Matrix::colSums(sparseM)
  counts_per_feature <- Matrix::rowSums(sparseM)
  features_per_cell  <- Matrix::colSums(sparseM>0)
  cells_per_feature  <- Matrix::rowSums(sparseM>0)

  # generating the data frames for downstream use
  info_per_cell <- data.frame(counts_per_cell = counts_per_cell) %>%
    tibble::rownames_to_column(var = "cell") %>%
    full_join(
      data.frame(features_per_cell = features_per_cell) %>%
        tibble::rownames_to_column(var = "cell"),
      by = "cell"
    )

  info_per_feature <- data.frame(counts_per_feature = counts_per_feature) %>%
    tibble::rownames_to_column(var = "feature") %>%
    full_join(
      data.frame(cells_per_feature = cells_per_feature) %>%
        tibble::rownames_to_column(var = "feature"),
      by = "feature"
    )

  # Add annotation overlap information to the feature information data frame
  features_in_matrix <- unique(feature.gr)[paste(seqnames(unique(feature.gr)), GenomicAlignments::ranges(unique(feature.gr)), sep=":") %in% info_per_feature$feature]

  pro.gr <- rtracklayer::import(promoters_file)
  enhs.gr <- rtracklayer::import(enhs_file)
  tss_df <- data.table::fread(tss_file, select=c(1:3), header = F, col.names = c("chr", "start", "end"))
  tss.gr <- GenomicRanges::makeGRangesFromDataFrame(tss_df)

  pro.overlaps <- GenomicRanges::findOverlaps(query = features_in_matrix,
                                                   subject = pro.gr,
                                                   type = "any",
                                                   ignore.strand = TRUE)
  enhs.overlaps <- GenomicRanges::findOverlaps(query = features_in_matrix,
                                              subject = enhs.gr,
                                              type = "any",
                                              ignore.strand = TRUE)
  tss.overlaps <- GenomicRanges::findOverlaps(query = features_in_matrix,
                                               subject = tss.gr,
                                               type = "any",
                                               ignore.strand = TRUE)


  pro.hits <- seq(length(features_in_matrix)) %in% S4Vectors::queryHits(pro.overlaps)
  enhs.hits <- seq(length(features_in_matrix)) %in% S4Vectors::queryHits(enhs.overlaps)
  tss.hits <- seq(length(features_in_matrix)) %in% S4Vectors::queryHits(tss.overlaps)

  info_per_feature <- cbind(info_per_feature,
                            promoter_overlaps = pro.hits,
                            enhancer_overlaps = enhs.hits,
                            tss_overlaps = tss.hits)

  if (!is.null(gene_anno_file) && file.exists(gene_anno_file)) {
    gene_anno.gr <- rtracklayer::import(gene_anno_file)
    gene.overlaps <- GenomicRanges::findOverlaps(query = features_in_matrix,
                                                 subject = gene_anno.gr,
                                                 type = "any",
                                                 ignore.strand = TRUE)
    gene.hits <- seq(length(features_in_matrix)) %in% S4Vectors::queryHits(gene.overlaps)
    intergenic <- !(pro.hits | enhs.hits | tss.hits | gene.hits)
    info_per_feature <- cbind(info_per_feature, gene_overlaps = gene.hits, intergenic = intergenic)
  }


  cat(
    paste0(
      "Feature quality control metrics produced at ",
      as.character(Sys.time()),
      "\n"
    ),
    file = log_file, append = TRUE)

  utils::write.csv(info_per_cell, paste0(log_and_stats_folder, "filtered_stats_per_cell.csv"), row.names = FALSE)
  utils::write.csv(info_per_feature, paste0(log_and_stats_folder, "filtered_stats_per_feature.csv"), row.names = FALSE)
  
  end_time = Sys.time()
  
  print(end_time - init_time)
  
  cat(
    paste0(
      "sc_atac_counting finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = TRUE)
  
  # create the sce object
  sc_atac_create_sce(input_folder = output_folder,
                            organism     = organism,
                            feature_type = feature_type,
                            pheno_data   = pheno_data,
                            report       = FALSE)
  
  # create the report if the option is given
  if(create_report){
    inputfolder <- dirname(insortedbam)
    sc_atac_create_report(input_folder = inputfolder,
                          output_folder= file.path(output_folder, "scPipe_atac_stats"),
                          sample_name  = sample_name,
                          organism     = organism,
                          feature_type = feature_type)
  }
  
}


