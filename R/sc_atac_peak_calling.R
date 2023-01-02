###############################
# Peak calling Using MACS3
###############################

#' sc_atac_peak_calling()
#'
#' @param inbam The input tagged, sorted, duplicate-free input BAM file
#' @param ref The reference genome file
#' @param genome_size The size of the genome
#' @param output_folder The path of the output folder
#'
#' @returns None (invisible `NULL`)
#' @examples
#' \dontrun{
#' sc_atac_peak_calling(
#'     inbam,
#'     reference) 
#' }
#' @export
#' 
#' @importFrom basilisk basiliskRun basiliskStart basiliskStop
sc_atac_peak_calling <- function(inbam, 
                                ref = NULL, 
                                genome_size = NULL,
                                output_folder = NULL) {
  
    if (is.null(ref) && is.null(genome_size)) {
        stop("No genome or genome size was specified. Must specify at least one!")
    } else if (is.null(genome_size)) {
        if (is.null(ref) || !file.exists(ref)) {
            stop("Specified reference file doesn't exist!")
        }
        genome_size <- signif(file.info(ref)$size, 2)
        message("Using genome size of", genome_size)
    }
    
    
    if(is.null(output_folder)) {
        output_folder <- file.path(getwd(), "scPipe-atac-output")
    }
    
    if (!dir.exists(output_folder)){
        dir.create(output_folder,recursive=TRUE)
        message("Output directory does not exist. Created: ", output_folder)
    }
    
    # Create log folder/file
    log_and_stats_folder <- file.path(output_folder, "scPipe_atac_stats")
    dir.create(log_and_stats_folder, showWarnings = FALSE)
    log_file <- file.path(log_and_stats_folder, "log_file.txt")
    if(!file.exists(log_file)) file.create(log_file)
    write(
        c(
        "sc_atac_peak_calling starts at ",
        as.character(Sys.time()),
        "\n"
        ),
        file = log_file, append = TRUE)
    
    # MACSr::callpeak(inbam, 
    #                 nomodel         = TRUE, 
    #                 shift           = 100, 
    #                 extsize         = 200, 
    #                 gsize           = genome_size, 
    #                 call_summits    = TRUE, 
    #                 store_bdg       = TRUE,
    #                 do_SPMR         = TRUE,
    #                 cutoff_analysis = TRUE,
    #                 outdir          = output_folder)
    
    proc <- basilisk::basiliskStart(scPipe_env)
    on.exit(basilisk::basiliskStop(proc))

    basilisk::basiliskRun(proc, function(inbam, genome_size, output) {
        macs <- reticulate::import("MACS3")
        system2("macs3", c("callpeak", "-t", 
                        inbam, 
                        "--nomodel", 
                        "--shift", 100, 
                        "--extsize", 200, 
                        "--g", genome_size, 
                        "--outdir", output_folder,
                        "--bdg", 
                        "--SPMR", 
                        "--cutoff-analysis",
                        "--call-summits"))
    }, inbam=inbam, genome_size=genome_size, output=output_folder)
    
    write(
        c(
        "sc_atac_peak_calling finishes at ",
        as.character(Sys.time()),
        "\n"
        ), file = log_file, append = TRUE)
}