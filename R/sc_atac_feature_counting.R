
#########################################################################
# Generating the Feature Matrix from a Given BAM file and a Feature FIle
#########################################################################

#' @name sc_atac_feature_counting
#' @title generating the feature by cell matrix 
#' @description feature matrix is created using a given demultiplexed BAM file and 
#' a selected feature type
#' @param fragment_file The fragment file
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
#' 
#' @returns None (invisible `NULL`)
#' 
#'
#' @examples
#' \dontrun{
#' sc_atac_feature_counting(
#'    fragment_file = fragment_file,
#'    cell_calling = "filter",
#'    exclude_regions = TRUE,
#'    feature_input = feature_file)
#' }    
#' @export
#' 
sc_atac_feature_counting <- function(
    fragment_file, 
    feature_input             = NULL, 
    bam_tags                  = list(bc="CB", mb="OX"), 
    feature_type              = "peak", 
    organism                  = "hg38", 
    cell_calling              = "filter",
    sample_name               = "",
    genome_size               = NULL, 
    promoters_file            = NULL,
    tss_file                  = NULL,
    enhs_file                 = NULL,
    gene_anno_file            = NULL,
    pheno_data                = NULL, 
    bin_size                  = NULL, 
    yieldsize                 = 1000000,
    n_filter_cell_counts      = 200,
    n_filter_feature_counts   = 10,
    exclude_regions           = FALSE, 
    excluded_regions_filename = NULL,
    output_folder             = NULL,
    fix_chr                   = "none", 
    lower                     = NULL,
    min_uniq_frags            = 3000,
    max_uniq_frags            = 50000,
    min_frac_peak             = 0.3,
    min_frac_tss              = 0,
    min_frac_enhancer         = 0,
    min_frac_promoter         = 0.1,
    max_frac_mito             = 0.15,
    create_report             = FALSE
    ) {
  
    . <- V1 <- V2 <- V3 <- init <- peakStart <- seqnames <- peakEnd <- CB <- chromosome <- chrS <- feature <- barcode <- NULL
    
    init_time <- Sys.time()
    
    available_organisms <- c("hg19",
                            "hg38",
                            "mm10")
    
    if(!is.null(organism)) stopifnot(organism %in% available_organisms)
    
    stopifnot(fix_chr %in% c("none", "excluded_regions", "feature", "both"))
    
    if(is.null(output_folder)) {
        output_folder <- file.path(getwd(), "scPipe-atac-output")
        message("Output directory has not been provided. Saving output in\n", output_folder)
    }
    
    if (!dir.exists(output_folder)){
        dir.create(output_folder,recursive=TRUE)
        message("Output directory does not exist. Creating...", output_folder)
    }
    
    
    # initiate log file and location for stats in scPipe_atac_stats
    log_and_stats_folder <- paste0(output_folder, "/scPipe_atac_stats/")
    dir.create(log_and_stats_folder, showWarnings = FALSE)
    log_file             <- paste0(log_and_stats_folder, "log_file.txt")
    stats_file           <- paste0(log_and_stats_folder, "stats_file_align.txt")
    if(!file.exists(log_file)) file.create(log_file)
    # file.create(stats_file)
    
    # timer
    write(
        c(
        "sc_atac_feature_counting starts at ",
        as.character(Sys.time()),
        "\n"
        ), 
        file = log_file, append = TRUE)

    
    
    ############# feature type is genome_bin ####################

    if(feature_type == 'genome_bin'){
        # TODO: test the format of the fasta file here and stop if not proper format
        message("`genome bin` feature type is selected for feature input. reading the genome fasta file ...")
        
        # index for feature_input is created if not found in the same directory as the feature_input
        if(!file.exists(paste0(feature_input,".fai"))){
            Rsamtools::indexFa(feature_input)
            message("Index for ", feature_input, " is being created... ")
        }
        
        if(is.null(bin_size)){
            bin_size <- 2000
            message("Default bin size of 2000 is selected")
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
                message("Generated the genome bin file:\n", out_bed_filename)
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
                
                
                
                out_bed <- purrr::map_df(seq_len(nrow(sizes_df_aux)), function(i){
                    
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
                message("Generated the genome bin file:\n", out_bed_filename)
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
        message("`peak`, `tss` or `gene` feature_type is selected for feature input")
        
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
                stop("Feature file provided does not contain 3 columns. Cannot append chr")
                # warning("Feature file provided does not contain 3 columns. Cannot append chr")
                # break;
            }
            
            
            rcpp_append_chr_to_bed_file(feature_input, out_bed_filename_feature)
            
            # Check if file was created
            if(file.exists(out_bed_filename_feature)) {
                message("Appended 'chr' to feature file and output created in:", out_bed_filename_feature)
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
                    stop("excluded_regions file provided does not contain 3 columns. Cannot append chr")
                #   warning("excluded_regions file provided does not contain 3 columns. Cannot append chr")
                #   break
                }
                
                rcpp_append_chr_to_bed_file(excluded_regions_filename, out_bed_filename_excluded_regions)
                
                # Check if file was created
                if(file.exists(out_bed_filename_excluded_regions)) {
                    message("Appended 'chr' to excluded_regions file and output created in:", out_bed_filename_excluded_regions)
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
    message("Creating Galignment object for the feature input ...")
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
                lines_to_keep            <- setdiff(seq_len(nrow(feature.gr.df)), lines_to_remove)
                feature.gr               <- feature.gr[lines_to_keep, ]
            } # End if(number_of_lines_to_remove > 0)
        
        } else{
            warning("Parameter exclude_regions was TRUE but no known organism or excluded_regions_filename provided. Proceding without excluding regions.")
        }
        
        
    } # End if(exclude_regions)
  

    # Check if organism is pre-recognized and if so then use the package's annotation files
    if (organism %in% c("hg19", "hg38", "mm10")) {
        message(organism, "is a recognized organism. Using annotation files in repository.")
        anno_paths <- system.file("extdata/annotations/", package = "scPipe", mustWork = TRUE)
        
        promoters_file <- file.path(anno_paths, paste0(organism, "_promoter.bed.gz"))
        tss_file <- file.path(anno_paths, paste0(organism, "_tss.bed.gz"))
        enhs_file <- file.path(anno_paths, paste0(organism, "_enhancer.bed.gz"))
    } else if (!all(file.exists(c(promoters_file, tss_file, enhs_file)))) {
        stop("One of the annotation files could not be located. Please make sure their paths are valid.")
    }

    # Create bins used for TSS enrichment plot
    tss_df <- data.table::fread(tss_file, select=c(seq_len(3)), header = FALSE, col.names = c("chr", "start", "end"))
    range <- 4000
    bin_size <- 100
    n_bins <- range/bin_size-1
    bins_df_all <- get_all_TSS_bins(tss_df, range, bin_size)
    bins_df_all.gr <- GenomicRanges::makeGRangesFromDataFrame(bins_df_all, keep.extra.columns=TRUE) 
    bin_hits <- c() # used to store all overlaps with bins
    
    ############## read in the aligned and demultiplexed BAM file
    
    # Helper utility function for adding matrices
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
    
    unique_feature.gr <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(dplyr::distinct(data.frame(feature.gr), seqnames, start, end)))
    #   unique_feature.gr <- data.frame(feature.gr) %>% dplyr::distinct(seqnames, start, end) %>% as.data.frame() %>% GenomicRanges::makeGRangesFromDataFrame()
    min_feature_width <- min(GenomicAlignments::ranges(feature.gr)@width)
    
    fragments <- data.table::fread(fragment_file, select=seq_len(5), header = FALSE, col.names = c("seqnames", "start", "end", "barcode", "count")) 
    fragments.gr <- fragments %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # Compute overlaps with feature
    peak.overlaps <- GenomicRanges::findOverlaps(query = fragments.gr,
                                                subject = unique_feature.gr,
                                                type = "any",
                                                ignore.strand = TRUE)
    
    # Compute overlaps with bins for TSS enrichment plot
    message("Generating TSS plot data")
    median_feature_overlaptss <- stats::median(GenomicAlignments::ranges(bins_df_all.gr)@width)
    maxgaptss                 <- 0.51*median_feature_overlaptss
    tss_bin_overlaps <- GenomicAlignments::findOverlaps(query = bins_df_all.gr,
                                                        subject = fragments.gr,
                                                        type = "any",
                                                        ignore.strand = TRUE)
    bin_hits <- S4Vectors::queryHits(tss_bin_overlaps) 
    
    num_overlaps <- length(peak.overlaps)
    
    barcodes <- fragments[S4Vectors::queryHits(peak.overlaps), ]$barcode
    features <- data.frame(unique_feature.gr[S4Vectors::subjectHits(peak.overlaps), ]) %>% 
        dplyr::select(seqnames, start, end) %>% 
        tidyr::unite("range", start:end, sep="-") %>% 
        tidyr::unite("feature", seqnames:range, sep=":") 
    
    # Chunk size
    chunk_size <- yieldsize
    i <- 1
    
    message("Generating feature-barcode fragment count matrix")
    # Iterate over chunks
    feature_matrix <- NULL
    while (i <= num_overlaps) {
        temp_mat <- as.matrix(table(feature = features$feature[i:(i + chunk_size-1)], 
                                    barcode = barcodes[i:(i + chunk_size-1)]))
        if (is.null(feature_matrix)) {
            feature_matrix <- Matrix::Matrix(temp_mat)
        } else {
            feature_matrix <- add_matrices(feature_matrix, Matrix::Matrix(temp_mat))
        }
        i <- i + chunk_size
    }
    
    average_number_of_lines_per_CB <- sum(feature_matrix)/ncol(feature_matrix)

    message("Average no. of fragments per barcode: ", average_number_of_lines_per_CB, "\n")

    
    saveRDS(feature_matrix, file = file.path(output_folder, "unfiltered_feature_matrix.rds"))
    message("Raw feature matrix generated: ", file.path(output_folder, "unfiltered_feature_matrix.rds") )
    
    
    ################ Initiate log file
    log_and_stats_folder       <- paste0(output_folder, "/scPipe_atac_stats/")
    dir.create(log_and_stats_folder, showWarnings = FALSE)
    log_file                   <- paste0(log_and_stats_folder, "log_file.txt")
    if(!file.exists(log_file)) file.create(log_file)

    write(c("Average number of reads per cell barcode:", average_number_of_lines_per_CB, "\n"),
        file = log_file, append = TRUE)

    write(c("Number of regions removed from feature input due to being invalid:", number_of_lines_to_remove, "\n"),
        file = log_file, append = TRUE)

    write(
        c(
        "Raw matrix generated at ",
        as.character(Sys.time()),
        "\n"
        ),
        file = log_file, append = TRUE)
    
    message("Calculating TSS enrichment scores")
    # Create a matrix where the rows are TSSs and the columns are the bins
    mat <- matrix(0, nrow(tss_df), n_bins)
    for (i in seq_len(length(bin_hits))) { # populate matrix with hits
        tss_index <- (bin_hits[i]-1) %/% n_bins
        bin <- (bin_hits[i]-1) %% n_bins
        mat[tss_index+1, bin+1] <- mat[tss_index+1, bin+1]+1
    }
    
    mat_non_zero <- mat[rowSums(mat) != 0,]
    
    # Calculate TSS enrichment scores
    tss_dists <- seq(-range/2+bin_size, range/2-bin_size, bin_size)
    flank_read_depth <- (mat_non_zero[,1]+mat_non_zero[,ncol(mat_non_zero)])/2 # Used to normalise
    norm_read_depths <- stats::na.omit(mat_non_zero/flank_read_depth) # also ignore rows where flanks have no overlaps
    
    aggregate_tss_scores <- colMeans(norm_read_depths)
    tsse <- max(aggregate_tss_scores)
    tss_plot_data <- data.frame(dists = tss_dists, 
                                agg_tss_scores = aggregate_tss_scores)
    
    utils::write.csv(tss_plot_data, file = file.path(log_and_stats_folder, "tss_plot_data.csv"))
    
    # generate quality control metrics for cells
    message("Generating QC metrics for cells\n")
    sc_atac_create_cell_qc_metrics(frags_file = fragment_file,
                                    peaks_file = feature_input,
                                    promoters_file = promoters_file,
                                    tss_file = tss_file,
                                    enhs_file = enhs_file,
                                    output_folder = output_folder)


    cell_qc_metrics_file <- file.path(output_folder, "cell_qc_metrics.csv")

    write(
        c(
        "Cell QC metrics generated at ",
        as.character(Sys.time()),
        "\n"
        ),
        file = log_file, append = TRUE)

    # Cell calling
    matrixData <- sc_atac_cell_calling(mat                  = feature_matrix,
                                        cell_calling         = cell_calling,
                                        output_folder        = output_folder,
                                        genome_size          = genome_size,
                                        cell_qc_metrics_file = cell_qc_metrics_file,
                                        lower                = lower,
                                        min_uniq_frags       = min_uniq_frags,
                                        max_uniq_frags       = max_uniq_frags,
                                        min_frac_peak        = min_frac_peak,
                                        min_frac_tss         = min_frac_tss,
                                        min_frac_enhancer    = min_frac_enhancer,
                                        min_frac_promoter    = min_frac_promoter,
                                        max_frac_mito        = max_frac_mito)

    
    write(
        c(
        "Cell calling completed at ",
        as.character(Sys.time()),
        "\n"
        ),
        file = log_file, append = TRUE)
    
    # filter the matrix based on counts per cell
    if (n_filter_cell_counts > 0) {
        message("Cells with less than ", n_filter_cell_counts, " counts are filtered out.")
        write(
            c(
                "cells with less than ",
                n_filter_cell_counts," counts are filtered out.",
                "\n"
            ),
            file = log_file, append = TRUE)
        filtered_indices  <- base::colSums(Matrix::as.matrix(matrixData), na.rm=TRUE) > n_filter_cell_counts
        message("Number of cells to remove:", sum(!filtered_indices))
        if (length(filtered_indices[filtered_indices == TRUE]) >= 10) { # only use if resulting matrix isn't too small
            matrixData <- as.matrix(matrixData[, filtered_indices]) # all the remaining columns
        } else {
            message("No cells were filtered out since otherwise there would be too few left.")
        }
    } else {
        message("No cells were filtered out based on counts.")
    }
    
    
    # filter the matrix based on counts per feature
    if (n_filter_feature_counts > 0) {
        message("features with less than ", n_filter_feature_counts, " counts are filtered out.")
        write(
            c(
                "features with less than ",
                n_filter_feature_counts," counts are filtered out.",
                "\n"
            ),
            file = log_file, append = TRUE)
        filtered_indices  <- base::rowSums(Matrix::as.matrix(matrixData), na.rm=TRUE) > n_filter_feature_counts
        message("Number of features to remove:", sum(!filtered_indices))
        if (length(filtered_indices[filtered_indices == TRUE]) >= 10) {
            matrixData <- matrixData[filtered_indices,] # all the remaining rows
        } else {
            message("No features were filtered out since otherwise there would be too few left.")
        }

    } else {
        message("No features were filtered out based on counts.")
    }

    # Update cell QC metrics to include whether the cell was kept or not
    cqc <- read.csv(file.path(output_folder, "cell_qc_metrics.csv"))
    cqc$cell_called <- cqc$bc %in% colnames(matrixData)
    write.csv(cqc, file.path(output_folder, "cell_qc_metrics.csv"), row.names = FALSE)
    
    
    # converting the NAs to 0s if the sparse option to create the sparse Matrix properly
    message("making sparse")
    sparseM <- Matrix::Matrix(matrixData, sparse=TRUE)
    # # add dimensions of the sparse matrix if available
    # if(cell_calling != FALSE){
    #   barcodes <- utils::read.table(paste0(output_folder, '/non_empty_barcodes.txt'))$V1
    #   features <- utils::read.table(paste0(output_folder, '/non_empty_features.txt'))$V1
    #   dimnames(sparseM) <- list(features, barcodes)
    # }

    message("Sparse matrix generated")
    saveRDS(sparseM, file = paste(output_folder, "/sparse_matrix.rds", sep = ""))
    # Matrix::writeMM(obj = sparseM, file=paste(output_folder, "/sparse_matrix.mtx", sep =""))
    message("Sparse count matrix is saved in\n", output_folder,"/sparse_matrix.mtx")

    # generate and save jaccard matrix
    # jaccardM <- locStra::jaccardMatrix(sparseM)
    # message("Jaccard matrix generated")
    # saveRDS(jaccardM, file = paste(output_folder,"/jaccard_matrix.rds",sep = ""))
    # message("Jaccard matrix is saved in\n", paste(output_folder,"/jaccard_matrix.rds",sep = "") )

    # generate and save the binary matrix
    matrixData[matrixData>0] <- 1
    saveRDS(matrixData, file = paste(output_folder,"/binary_matrix.rds",sep = ""))
    message("Binary matrix is saved in:\n", output_folder,"/binary_matrix.rds")
    
    ###### calculate the TF-IDF matrices here : TO-DO

    # following can be used to plot the stats and load it into sce object ########
    # (from https://broadinstitute.github.io/2020_scWorkshop/data-wrangling-scrnaseq.html)
    counts_per_cell    <- Matrix::colSums(sparseM)
    counts_per_feature <- Matrix::rowSums(sparseM)
    features_per_cell  <- Matrix::colSums(sparseM>0)
    cells_per_feature  <- Matrix::rowSums(sparseM>0)

    # generating the data frames for downstream use
    info_per_cell <- data.frame(counts_per_cell = counts_per_cell) %>%
        tibble::rownames_to_column(var = "cell") %>%
        dplyr::full_join(
        data.frame(features_per_cell = features_per_cell) %>%
            tibble::rownames_to_column(var = "cell"),
        by = "cell"
        )

    info_per_feature <- data.frame(counts_per_feature = counts_per_feature) %>%
        tibble::rownames_to_column(var = "feature") %>%
        dplyr::full_join(
        data.frame(cells_per_feature = cells_per_feature) %>%
            tibble::rownames_to_column(var = "feature"),
        by = "feature"
        )


    # Add annotation overlap information to the feature information data frame
    message("Computing feature QC metrics")
    features_in_matrix <- unique_feature.gr[paste(GenomicRanges::seqnames(unique_feature.gr), GenomicRanges::ranges(unique_feature.gr), sep=":") %in% info_per_feature$feature]

    pro.gr <- rtracklayer::import(promoters_file)
    enhs.gr <- rtracklayer::import(enhs_file)
    tss_df <- data.table::fread(tss_file, select=c(seq_len(3)), header = FALSE, col.names = c("chr", "start", "end"))
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
    

    write(
        c(
        "Feature quality control metrics produced at ",
        as.character(Sys.time()),
        "\n"
        ),
        file = log_file, append = TRUE)

    utils::write.csv(info_per_cell, paste0(log_and_stats_folder, "filtered_stats_per_cell.csv"), row.names = FALSE)
    utils::write.csv(info_per_feature, paste0(log_and_stats_folder, "filtered_stats_per_feature.csv"), row.names = FALSE)
    
    message("writing to csv")
    
    write(
        c(
        "sc_atac_feature_counting completed at ",
        as.character(Sys.time()),
        "\n\n"
        ), 
        file = log_file, append = TRUE)
    
    # create the sce object
    sc_atac_create_sce(input_folder = output_folder,
                                organism     = organism,
                                sample       = sample_name,
                                feature_type = feature_type,
                                pheno_data   = pheno_data,
                                report       = create_report)

    
    end_time <- Sys.time()
    message("sc_atac_feature_counting completed in ", difftime(end_time, init_time, units = "secs")[[1]], " seconds")
}
