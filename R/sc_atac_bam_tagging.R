#############################
# Demultiplxing FASTQ Reads
#############################

#' sc_atac_bam_tagging()
#' @name sc_atac_bam_tagging
#' @title BAM tagging
#' @description Demultiplexes the reads
#'
#' @param inbam The input BAM file
#' @param output_folder The path of the output folder
#' @param bc_length The length of the cellular barcodes
#' @param bam_tags The BAM tags
#' @param nthreads The number of threads
#'
#' @returns file path of the resultant demmultiplexed BAM file. 
#' @examples
#' r1 <- system.file("extdata", "small_chr21_R1.fastq.gz", package="scPipe")
#' r2 <- system.file("extdata", "small_chr21_R3.fastq.gz", package="scPipe")
#' barcode_fastq <- system.file("extdata", "small_chr21_R2.fastq.gz", package="scPipe")
#' out <- tempdir()
#'
#' sc_atac_trim_barcode(r1=r1, r2=r2, bc_file=barcode_fastq, output_folder=out)
#'
#' demux_r1 <- file.path(out, "demux_completematch_small_chr21_R1.fastq.gz")
#' demux_r2 <- file.path(out, "demux_completematch_small_chr21_R3.fastq.gz")
#' reference <- system.file("extdata", "small_chr21.fa", package="scPipe")
#'
#' aligned_bam <- sc_aligning(ref=reference, R1=demux_r1, R2=demux_r2, nthreads=6, output_folder=out)
#'
#' out_bam <- sc_atac_bam_tagging(
#'    inbam = aligned_bam,
#'    output_folder = out,
#'    nthreads = 6)
#'
#' @export
sc_atac_bam_tagging <- function(inbam,
                                output_folder = NULL,
                                bc_length = NULL,
                                bam_tags = list(bc="CB", mb="OX"),
                                nthreads = 1) {

    if (any(!file.exists(inbam))) {
        stop("At least one input bam file should be present")
    } else {
        inbam <- path.expand(inbam)
    }

    if(is.null(output_folder)){
        output_folder <- file.path(getwd(), "scPipe-atac-output")
        # output_folder <- paste0(sub('[/][^/]+$', '', inbam)) # same as output_folder <- basename(inbam)
    }

    if (!dir.exists(output_folder)){
        dir.create(output_folder,recursive=TRUE)
        message("Output directory does not exist. Creating ", output_folder)
    }

    get_filename_without_extension <- function(path) {
        vec <- strsplit(basename(path), "\\.")[[1]]
        return(paste(vec[seq_len(length(vec)-1)], collapse = "."))
    }
    fileNameWithoutExtension <- get_filename_without_extension(inbam)

    outbam                   <- file.path(output_folder, paste0(fileNameWithoutExtension, "_tagged_sorted.bam"))
    outsortedbam             <- file.path(output_folder, paste0(fileNameWithoutExtension, "_tagged_sorted")) # don't want extension

    log_and_stats_folder       <- file.path(output_folder, "scPipe_atac_stats")
    dir.create(log_and_stats_folder, showWarnings = FALSE)
    log_file                   <- file.path(log_and_stats_folder, "log_file.txt")
    stats_file                 <- file.path(log_and_stats_folder, "stats_file_bam_tagging.txt")
    if(!file.exists(log_file)) file.create(log_file)

    write(
        c(
            "sc_atac_tagging began at ",
            as.character(Sys.time()),
            "\n"
        ),
        file = log_file, append = TRUE)

    outbam <- path.expand(outbam)
    if(!file.exists(outbam)){
        file.create(outbam)

        rcpp_sc_atac_bam_tagging(inbam, outbam, bam_tags$bc, bam_tags$mb, nthreads)

        # Rsamtools::sortBam(outbam, outsortedbam, indexDestination = TRUE, maxMemory = 1024/32*max_memory)
        Rsamtools::indexBam(paste0(outsortedbam, ".bam"))

        message("Demultiplexed BAM file sorted and indexed ...")
    }
    #cat(paste0(outsortedbam, ".bam"))
    #cat("\n")


    if(is.null(bc_length)){
        message("Using default value for barcode length (bc_length = 16) ")
        bc_length <- 16
    }


    flag_defs <- tibble::tibble(
        type =
            c("one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "one_read_unmapped", "both_reads_unmapped", "both_reads_unmapped", "mapped", "mapped", "mapped", "mapped", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_wrong_orientation", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously", "mapped_ambigously")
        ,
        flag =
            c(73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137, 77, 141, 99, 147, 83, 163, 67, 131, 115, 179, 81, 161, 97, 145, 65, 129, 113, 177))

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

    message("Generating mapping statistics per barcode")
    message("Iterating through 5,000,000 reads at a time")
    bamfl <- open(Rsamtools::BamFile(outbam, yieldSize = 5000000))
    params <- Rsamtools::ScanBamParam(what=c("flag"), tag=c("CB"))
    iter <- 1
    full_matrix <- NULL

    # Iterate over the BAM file in chunks of 5 million reads
    while(length((bam0 <- Rsamtools::scanBam(bamfl, param = params)[[1]])$flag)) {
        message("chunk", iter)

        # Create a data.table object with each row representing a read, and the columns as the barcode and flag
        df <- data.table::setDT(data.frame(bam0$tag$CB, bam0$flag) %>% data.table::setnames(c("barcode", "flag")) %>% dplyr::left_join(flag_defs, by = "flag"))[, !"flag"]

        # Count the reads per barcode
        x <- df[, list(count=.N), names(df)]
        x <- data.table::dcast(x, barcode ~ type, value.var = "count")
        for (i in seq_along(x)) data.table::set(x, i=which(is.na(x[[i]])), j=i, value=0)

        if (is.null(full_matrix)) {
            full_matrix <- as.matrix(x, rownames=TRUE)
        } else {
            full_matrix <- add_matrices(full_matrix, as.matrix(x, rownames=TRUE))
        }
        iter <- iter+1
    }
    df <- data.table::setnames(data.table::setDT(as.data.frame(as.matrix(full_matrix)), keep.rownames = TRUE), "rn", "barcode")
    df[ ,count := rowSums(.SD), .SDcols = names(df[,!"barcode"])]
    mapping_stats_path <- file.path(log_and_stats_folder, "mapping_stats_per_barcode.csv")
    data.table::fwrite(df, mapping_stats_path)
    message("Completed generating mapping statistics per barcode, saved to ", file.path(log_and_stats_folder, "mapping_stats_per_barcode.csv"), "\n")
    write(
        c(
            "sc_atac_tagging() completed at ",
            as.character(Sys.time()),
            "\n\n"
        ),
        file = log_file, append = TRUE)

  
    message("The output tagged and sorted BAM file was sent to ", output_folder)
  
    return(outbam)
}
