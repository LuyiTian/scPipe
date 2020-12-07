library(rtracklayer)
library(GenomicAlignments)
library(Rsamtools)
library(SummarizedExperiment) 
library(DelayedArray)
library(BiocParallel)
library(matrixStats)
library(Biobase)
library(GenomicRanges)
library(GenomeInfoDb)
library(IRanges)            
library(S4Vectors)
library(BiocGenerics)  









#' TODO: documentation
#'
#' @export
#' 
sc_atac_feature_counting <- function(
  insortedbam, 
  feature_input, 
  bam_tags       = list(bc="CB", mb="OX"), 
  feature_type   = "peak", 
  organism       = NULL,
  cell_calling   = "cellranger", # either c("cellranger", "emptydrops", "filter")
  genome_size    = NULL, # this is optional but needed if the cell_calling option is cellranger AND organism in NULL
  qc_per_bc_file = NULL, # this is optional but needed if the cell_calling option is cellranger
  bin_size       = NULL, 
  yieldsize      = 1000000,
  mapq           = 0,
  blacklist      = NULL, # use a better term, e.g. excluded_regions
  output_folder  = "",
  fix_chr        = "none" # should be either one of these: c("none", "blacklist", "feature", "both")
){
  
  init_time = Sys.time()
  
  stopifnot(fix_chr %in% c("none", "blacklist", "feature", "both"))
  
  if(output_folder == ''){
    output_folder <- file.path(getwd(), "scPipe-atac-output")
  }
  
  if (!dir.exists(output_folder)){
    dir.create(output_folder,recursive=TRUE)
    cat("Output directory does not exist. Saving output in\n", output_folder, "\n")
  }
  
  
  # initate log file and location for stats in scPipe_atac_stats
  log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
  dir.create(log_and_stats_folder, showWarnings = F)
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
    file = log_file, append = T)
  

  if(feature_type == 'genome_bin'){
    # TODO: test the format of the fasta file here and stop if not proper format
    cat("`genome bin` feature type is selected for feature input. reading the genome fasta file ...", "\n")
    
    # TODO: execute following if the index for fasta is not found...
    Rsamtools::indexFa(feature_input)
    cat("Index for ", feature_input, " not found. Creating one now... ", "\n")
    
    if(is.null(bin_size)){
      bin_size <- 2000
      cat("Default bin size of 2000 is selected", "\n")
    }
    
    out_bed_filename <- paste0(output_folder, "/", sub('\\..[^\\.]*$', '', basename(feature_input)), ".bed") # remove extension and append output folder
    
    if(file.exists(out_bed_filename)) {
      warning(out_bed_filename, " file already exists. Replacing!")
      file.remove(out_bed_filename)
    }
    
    rcpp_fasta_bin_bed_file(feature_input, out_bed_filename, bin_size)
    
    # Check if file was created
    if(file.exists(out_bed_filename)) {
      cat("Generated the genome bin file:\n", out_bed_filename, "\n")
    }
    else {
      stop("File ", out_bed_filename, "file was not created.")
      }
    # feature_input <- genome_bin
    }
  
  cat("Creating GAlignment object for the sorted BAM file...\n")

  param <- Rsamtools::ScanBamParam(tag = as.character(bam_tags),  mapqFilter=mapq)
  bamfl <- Rsamtools::BamFile(insortedbam, yieldSize = yieldsize)
  open(bamfl)
  
  yld                            <- GenomicAlignments::readGAlignments(bamfl,use.names = TRUE, param = param)
  yld.gr                         <- makeGRangesFromDataFrame(yld,keep.extra.columns=TRUE) 
  #yld.gr                        <- as(GRanges(yld), "GAlignments")
  average_number_of_lines_per_CB <- length(yld.gr$CB)/length(unique(yld.gr$CB))
  
  saveRDS(yld.gr, file = paste(output_folder,"/BAM_GAlignmentsObject.rds",sep = ""))
  cat("GAlignment object is created and saved in \n", paste(output_folder,"/BAM_GAlignmentsObject.rds",sep = "") , "\n")
  
  # generate the GAalignment file from the feature_input file
  cat("Creating Galignment object for the feature input...\n")
  feature.gr <- rtracklayer::import(feature_input)
  
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
      
      # Try to read first 5 rows of blacklist file to see if the format is correct
      feature_head <- read.table(feature_input, nrows = 5)
      if(ncol(feature_head) != 3){
        warning("Feature file provided does not contain 3 columns. Cannot append chr")
        break
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
  
  
  if(fix_chr %in% c("blacklist", "both")){
    if(!is.null(blacklist)){
      
      out_bed_filename_blacklist <- paste0(output_folder, "/", 
                                          sub('\\..[^\\.]*$', 
                                              '', 
                                              basename(blacklist)
                                              ), "_fixedchr.bed"
                                          ) # remove extension and append output folder
     
       # Try to read first 5 rows of blacklist file to see if the format is correct
      blacklist_head <- read.table(blacklist, nrows = 5)
      if(ncol(blacklist_head) != 3){
        warning("Blacklist file provided does not contain 3 columns. Cannot append chr")
        break
      }
      
      rcpp_append_chr_to_bed_file(blacklist, out_bed_filename_blacklist)
      
      # Check if file was created
      if(file.exists(out_bed_filename_blacklist)) {
        cat("Appended 'chr' to blacklist file and output created in:", out_bed_filename_blacklist, "\n")
      }
      else {
        stop("File ", out_bed_filename_blacklist, "file was not created.")
        }
    }
    
  } # end if blacklist
  }
  
  ############################### end fix_chr
  
  # ______________ need to add the blacklist filtering here , relevant param has been added to the function i.e. blacklist.
  # this param should take either a tab delimited file in the format of feature_input or keywords "hs", "mm" which would load a stored GAlignment file related to them
  
  number_of_lines_to_remove    <- 0
  if(!is.null(blacklist)){
    blacklist.gr               <- rtracklayer::import(blacklist)
    
    overlaps_blacklist_feature <- findOverlaps(blacklist.gr, feature.gr, maxgap = -1L, minoverlap = 0L) # Find overlaps
    lines_to_remove            <- as.data.frame(overlaps_blacklist_feature)$subjectHits # Lines to remove in feature file
    number_of_lines_to_remove  <- length(lines_to_remove)
    
    if(number_of_lines_to_remove > 0){ # If there are lines to remove
      feature.gr.df            <- as.data.frame(feature.gr)
      lines_to_keep            <- setdiff(1:nrow(feature.gr.df), lines_to_remove)
      feature.gr               <- feature.gr[lines_to_keep, ]
    }
  }
  
  # Log file
  log_file_name <- paste0(output_folder, "/log_file.txt")
  file.create(log_file_name)
  
  cat("Average number of reads per CB:", average_number_of_lines_per_CB, "\n", 
      file = log_file_name, append = T)
  
  cat("Number of regions removed from feature_input:", number_of_lines_to_remove, "\n", 
      file = log_file_name, append = T)
  
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
  
  #yld.gr[-queryHits(findOverlaps(yld.gr, feature.gr, type="any", ignore.strand = TRUE)),] 
  #feature.gr[-queryHits(findOverlaps(feature.gr, yld.gr, type="any", ignore.strand = TRUE)),] 
  # mergeByOverlaps(feature.gr, yld.gr)
  
  #mcols(yld.gr)[queryHits(overlaps), "peakStart"] <- start(ranges(feature.gr)[subjectHits(overlaps)])
  #mcols(yld.gr)[queryHits(overlaps), "peakEnd"]   <- end(ranges(feature.gr)[subjectHits(overlaps)])
  
  mcols(yld.gr)[subjectHits(overlaps), "peakStart"] <- start(ranges(feature.gr)[queryHits(overlaps)])
  mcols(yld.gr)[subjectHits(overlaps), "peakEnd"]   <- end(ranges(feature.gr)[queryHits(overlaps)])
  
  #is removing NAs here the right thing to do?
  #overlap.df <- data.frame(yld.gr) %>% dplyr::select(seqnames, peakStart, peakEnd, CB)
  overlap.df <- data.frame(yld.gr) %>% filter(!is.na(peakStart)) %>% dplyr::select(seqnames, peakStart, peakEnd, CB)
  
  # # this method still does nto assign the dims to the matrix properly, hence comented out.
  # matrix_rownames = overlap.df %>% 
  #   dplyr::select(1:3) %>% 
  #   distinct() %>% 
  #   mutate(name = paste0(seqnames, ":", peakStart, "-", peakEnd)) %>% 
  #   pull(name)
  # 
  # matrixData <- overlap.df %>% 
  #   dplyr::group_by(seqnames, peakStart, peakEnd, CB) %>% 
  #   dplyr::summarise(count = n()) %>% 
  #   purrr::set_names(c("chromosome","start","end","barcode","count")) %>% 
  #   unite("chrS", chromosome:start, sep=":") %>%
  #   unite("feature", chrS:end, sep="-") %>% 
  #   dplyr::group_by(feature,barcode) %>% 
  #   dplyr::mutate(grouped_id = row_number()) %>% 
  #   tidyr::spread(barcode, count) %>% 
  #   dplyr::select(-grouped_id) %>% 
  #   # TODO: sanity check to see if all features are unique here
  #   #if(!unique(matrixData$feature)) {stop("There are duplicate values in the feature input. Please check and rerun this step again")}
  #   as_tibble(rownames = "feature") %>% 
  #   dplyr::select(-1) %>% 
  #   as.matrix(.)
  # 
  # row.names(matrixData) = matrix_rownames
  
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
    as.matrix()
  
  # it still bugs me how there are no dimnames in the Matrix. Trying to add them manually here
  dimnames(matrixData)  <-  list(matrixData.old[1] %>% rownames(), matrixData.old %>% dplyr::select(-1) %>% colnames())

  # call sc_atac_cell_callling.R here ... still ongoing
  sc_atac_cell_calling(mat = matrixData, cell_calling = cell_calling, output_folder = output_folder)
  
  saveRDS(matrixData, file = paste(output_folder,"/feature_matrix.rds",sep = ""))
  cat("Feature matrix generated: ", paste(output_folder,"/feature_matrix.rds",sep = "") , "\n")
  
  # here you need to convert the NAs to 0s if the sparse option to create the sparse Matrix properly
  sparseM <- Matrix(matrixData%>%replace(is.na(.), 0), sparse=TRUE)
  cat("Sparse matrix generated", "\n")
  writeMM(obj = sparseM, file=paste(output_folder,"/sparse_matrix.mtx", sep =""))
  cat("Sparse count matrix is saved in\n", paste(output_folder,"/sparse_matrix.mtx",sep = "") , "\n")
  #saveRDS(, file = paste(output_folder,"/sparse_matrix.rds",sep = ""))
  
  jaccardM <- jaccardMatrix(sparseM)
  cat("Jaccard matrix generated", "\n")
  saveRDS(jaccardM, file = paste(output_folder,"/jaccard_matrix.rds",sep = ""))
  cat("Jaccard matrix is saved in\n", paste(output_folder,"/jaccard_matrix.rds",sep = "") , "\n")
  
  matrixData[matrixData>0] <- 1
  saveRDS(matrixData, file = paste(output_folder,"/binary_matrix.rds",sep = ""))
  cat("Binary matrix is saved in:\n", paste(output_folder,"/binary_matrix.rds",sep = "") , "\n")
  
  # following can be used to plot the stas and load it into sce object ########
  # (from https://broadinstitute.github.io/2020_scWorkshop/data-wrangling-scrnaseq.html)
  counts_per_cell    <- Matrix::colSums(sparseM)
  counts_per_feature <- Matrix::rowSums(sparseM)
  features_per_cell  <- Matrix::colSums(sparseM>0)
  cells_per_feature  <- Matrix::rowSums(sparseM>0)
  
  # plots
  #hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
  #hist(log10(features_per_cell+1), main='features per cell', col='wheat')
  #hist(log10(counts_per_feature+1), main='counts per feature', col='wheat')
  
  #plot(counts_per_cell, features_per_cell, log='xy', col='wheat')
  #title('counts vs features per cell')
  
  # plot(sort(features_per_cell), xlab='cell', log='y', main='features per cell (ordered)')
  
  # plotting examples end here #############
  
  end_time = Sys.time()

  print(end_time - init_time)
  
  cat(
    paste0(
      "sc_atac_counting finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = T)
  
  }

