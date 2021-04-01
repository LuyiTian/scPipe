
#########################################################################
# Generating the Feature Matrix from a Given BAM file and a Feature FIle
#########################################################################

#' @name sc_atac_feature_counting
#' @title generating the feature by cell matrix 
#' @description feature matrix is created using a given demultiplexed BAM file and 
#' a selected feature type
#' @param 
#'
#' @export
#' 

sc_atac_feature_counting <- function(
  insortedbam, 
  feature_input  = NULL, 
  bam_tags       = list(bc="CB", mb="OX"), 
  feature_type   = "peak", 
  organism       = NULL,
  cell_calling   = FALSE, # either c("cellranger", "emptydrops", "filter")
  genome_size    = NULL, # this is optional but needed if the cell_calling option is cellranger AND organism in NULL
  qc_per_bc_file = NULL, # this is optional but needed if the cell_calling option is cellranger
  bin_size       = NULL, 
  yieldsize      = 1000000,
  mapq           = 0,
  exclude_regions= FALSE, 
  excluded_regions_filename = NULL,
  output_folder  = "",
  fix_chr        = "none" # should be either one of these: c("none", "excluded_regions", "feature", "both")
){
  
  init_time = Sys.time()
  
  available_organisms = c("hg19",
                          "hg38",
                          "mm10")
  
  if(!is.null(organism)) stopifnot(organism %in% available_organisms)
  
  stopifnot(fix_chr %in% c("none", "excluded_regions", "feature", "both"))
  
  if(output_folder == ''){
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
      "sc_atac_counting starts at ",
      as.character(Sys.time()),
      "\n"
    ), 
    file = log_file, append = TRUE)
  
  
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
        # Do something
        # Use organism data sizes saved in repository
        
        organism_files     <- list.files(system.file("extdata/organism_data/", package = "scPipe", mustWork = TRUE))
        sizes_filename_aux <- grep(pattern = organism, x = organism_files, value = TRUE) %>% 
          grep(pattern     <- "size", x = ., value = TRUE)
        
        sizes_filename <- system.file(paste0("extdata/organism_data/", sizes_filename_aux), package = "scPipe", mustWork = TRUE)
        
        sizes_df <- read.table(sizes_filename, header = FALSE, col.names = c("V1", "V2"))
        
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
        
        write.table(out_bed, file = out_bed_filename, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
        
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
  
  # read in the aligned and demltiplexed BAM file
  
  cat("Creating GAlignment object for the sorted BAM file...\n")
  
  param <- Rsamtools::ScanBamParam(tag = as.character(bam_tags),  mapqFilter=mapq)
  bamfl <- Rsamtools::BamFile(insortedbam, yieldSize = yieldsize)
  open(bamfl)
  
  yld                            <- GenomicAlignments::readGAlignments(bamfl,use.names = TRUE, param = param)
  yld.gr                         <- makeGRangesFromDataFrame(yld,keep.extra.columns=TRUE) 
  average_number_of_lines_per_CB <- length(yld.gr$CB)/length(unique(yld.gr$CB))
  
  cat("Adjusting for the Tn5 cut site...\n")
  # to do
  
  saveRDS(yld.gr, file = paste(output_folder,"/BAM_GAlignmentsObject.rds",sep = ""))
  cat("GAlignment object is created and saved in \n", paste(output_folder,"/BAM_GAlignmentsObject.rds",sep = "") , "\n")
  
  # generate the GAalignment file from the feature_input file
  cat("Creating Galignment object for the feature input...\n")
  if(feature_type != 'genome_bin'){
  feature.gr <- rtracklayer::import(feature_input)
  } else {
    feature.gr <- rtracklayer::import(out_bed_filename)
  }
  
  ############################### fix_chr
  
  if(feature_type == 'peak'){
    cat("`peak` feature_type is selected for feature input", "\n")
    
    if(fix_chr %in% c("feature", "both")){
      
      out_bed_filename_feature <- paste0(output_folder, "/", 
                                         sub('\\..[^\\.]*$', 
                                             '', 
                                             basename(feature_input)
                                         ), "_fixedchr.bed"
      ) # remove extension and append output folder
      
      # Try to read first 5 rows of excluded_regions file to see if the format is correct
      feature_head <- read.table(feature_input, nrows = 5)
      if(ncol(feature_head) != 3){
        warning("Feature file provided does not contain 3 columns. Cannot append chr")
        break;
      }
      
      
      rcpp_append_chr_to_bed_file(feature_input, out_bed_filename_feature)
      
      # Check if file was created
      if(file.exists(out_bed_filename_feature)) {
        cat("Appended 'chr' to feature file and output created in:", out_bed_filename_feature, "\n")
      }
      else {
        stop("File ", out_bed_filename_feature, "file was not created.")
      }
      
    }
    
    
    if(fix_chr %in% c("excluded_regions", "both")){
      if(!is.null(excluded_regions_filename)){
        
        out_bed_filename_excluded_regions <- paste0(output_folder, "/", 
                                                    sub('\\..[^\\.]*$', 
                                                        '', 
                                                        basename(excluded_regions_filename)
                                                    ), "_fixedchr.bed"
        ) # remove extension and append output folder
        
        # Try to read first 5 rows of excluded_regions file to see if the format is correct
        excluded_regions_head <- read.table(excluded_regions_filename, nrows = 5)
        if(ncol(excluded_regions_head) != 3){
          warning("excluded_regions file provided does not contain 3 columns. Cannot append chr")
          break
        }
        
        rcpp_append_chr_to_bed_file(excluded_regions_filename, out_bed_filename_excluded_regions)
        
        # Check if file was created
        if(file.exists(out_bed_filename_excluded_regions)) {
          cat("Appended 'chr' to excluded_regions file and output created in:", out_bed_filename_excluded_regions, "\n")
        }
        else {
          stop("File ", out_bed_filename_excluded_regions, "file was not created.")
        }
      }
      
    } # end if excluded_regions
  }
  
  ############################### end fix_chr
  
  number_of_lines_to_remove <- 0
  
  if(exclude_regions){
    
    if(is.null(excluded_regions_filename) & !is.null(organism)){
      ## If excluded_regions_filename is null but organism is not, then read the file from system
      
      organism_files                <- list.files(system.file("data/extdata/annotations/", package = "scPipe", mustWork = TRUE))
      excluded_regions_filename_aux <- grep(pattern = organism, x = organism_files, value = TRUE) %>% 
        grep(pattern = "blacklist", x = ., value = TRUE)
      excluded_regions_filename     <- system.file(paste0("data/extdata/annotations/", excluded_regions_filename_aux), package = "scPipe", mustWork = TRUE)
    } 
    
    if(!is.null(excluded_regions_filename)){
      excluded_regions.gr               <- rtracklayer::import(excluded_regions_filename)
      overlaps_excluded_regions_feature <- findOverlaps(excluded_regions.gr, feature.gr, maxgap = -1L, minoverlap = 0L) # Find overlaps
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
  
  
  # Log file
  log_and_stats_folder       <- paste0(output_folder, "/scPipe_atac_stats/")
  dir.create(log_and_stats_folder, showWarnings = FALSE)
  log_file                   <- paste0(log_and_stats_folder, "log_file.txt")
  if(!file.exists(log_file)) file.create(log_file)
  
  cat("Average number of reads per CB:", average_number_of_lines_per_CB, "\n", 
      file = log_file, append = TRUE)
  
  cat("Number of regions removed from feature_input:", number_of_lines_to_remove, "\n", 
      file = log_file, append = TRUE)
  
  # Overlaps
  median_feature_overlap <- median(ranges(feature.gr)@width)
  minoverlap             <- 0.51*median_feature_overlap
  maxgap                 <- 0.51*median_feature_overlap
  
  overlaps               <- GenomicAlignments::findOverlaps(query         = feature.gr, # feature.gr,
                                                            subject       = yld.gr, # yld.gr, #
                                                            type          = "equal", 
                                                            maxgap        = maxgap, 
                                                            #minoverlap   = minoverlap, 
                                                            ignore.strand = TRUE)
  
  # generate the matrix using this overlap results above.
  
  mcols(yld.gr)[subjectHits(overlaps), "peakStart"] <- start(ranges(feature.gr)[queryHits(overlaps)])
  mcols(yld.gr)[subjectHits(overlaps), "peakEnd"]   <- end(ranges(feature.gr)[queryHits(overlaps)])
  
  #is removing NAs here the right thing to do?
  overlap.df <- data.frame(yld.gr) %>% filter(!is.na(peakStart)) %>% dplyr::select(seqnames, peakStart, peakEnd, CB)
  
  overlap.df <- overlap.df %>% 
    dplyr::group_by(seqnames, peakStart, peakEnd, CB) %>% 
    dplyr::summarise(count = n()) %>% 
    purrr::set_names(c("chromosome","start","end","barcode","count")) %>% 
    unite("chrS", chromosome:start, sep=":") %>%
    unite("feature", chrS:end, sep="-")
  
  matrixData <- overlap.df %>%
    dplyr::group_by(feature,barcode) %>% 
    tidyr::spread(barcode, count)
  
  matrixData           <- as.data.frame(matrixData)
  rownames(matrixData) <- matrixData$feature
  
  matrixData.old       <- matrixData
  matrixData           <- matrixData %>%
    dplyr::select(-1) %>%
    data.table::as.data.table() %>%
    as.matrix() %>%
    replace(is.na(.), 0)
  
  # add dimensions of the matrix
  dimnames(matrixData)  <-  list(matrixData.old[1] %>% rownames(), matrixData.old %>% dplyr::select(-1) %>% colnames())
  
  # call sc_atac_cell_callling.R here : emptyDrops() function currently implemented
  if(cell_calling =="emptydrops" || cell_calling =="cellranger" || cell_calling =="filter"){
    cat("calling `EmptyDrops` function for cell calling ... \n")
    sc_atac_cell_calling(mat = matrixData, cell_calling = 'emptydrops', output_folder = output_folder)
  }
  
  saveRDS(matrixData, file = paste(output_folder,"/feature_matrix.rds",sep = ""))
  cat("Feature matrix generated: ", paste(output_folder,"/feature_matrix.rds",sep = "") , "\n")
  
  # converting the NAs to 0s if the sparse option to create the sparse Matrix properly
  sparseM <- Matrix(matrixData, sparse=TRUE)
  cat("Sparse matrix generated", "\n")
  writeMM(obj = sparseM, file=paste(output_folder,"/sparse_matrix.mtx", sep =""))
  cat("Sparse count matrix is saved in\n", paste(output_folder,"/sparse_matrix.mtx",sep = "") , "\n")

  # generate and save jaccard matrix
  jaccardM <- jaccardMatrix(sparseM)
  cat("Jaccard matrix generated", "\n")
  saveRDS(jaccardM, file = paste(output_folder,"/jaccard_matrix.rds",sep = ""))
  cat("Jaccard matrix is saved in\n", paste(output_folder,"/jaccard_matrix.rds",sep = "") , "\n")
  
  # generate and save the binary matrix
  matrixData[matrixData>0] <- 1
  saveRDS(matrixData, file = paste(output_folder,"/binary_matrix.rds",sep = ""))
  cat("Binary matrix is saved in:\n", paste(output_folder,"/binary_matrix.rds",sep = "") , "\n")
  
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
  
  info_per_feature = data.frame(counts_per_feature = counts_per_feature) %>% 
    tibble::rownames_to_column(var = "feature") %>% 
    full_join(
      data.frame(cells_per_feature = cells_per_feature) %>% 
        tibble::rownames_to_column(var = "feature"),
      by = "feature"
    )
  

  write.csv(info_per_cell, paste0(log_and_stats_folder, "filtered_stats_per_cell.csv"), row.names = FALSE)
  write.csv(info_per_feature, paste0(log_and_stats_folder, "filtered_stats_per_feature.csv"), row.names = FALSE)
  
  end_time = Sys.time()
  
  print(end_time - init_time)
  
  cat(
    paste0(
      "sc_atac_counting finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = TRUE)
  
}

