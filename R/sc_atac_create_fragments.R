#####################################################################
# Generating Fragments for Aligned and Demultiplexed scATAC-Seq Reads
#####################################################################

#' @name sc_atac_create_fragments
#' @title Generating the popular fragments for scATAC-Seq data using sinto
#' @description Takes in a tagged and sorted BAM file and outputs the associated fragments in a .bed file
#' 
#' @param inbam The tagged, sorted and duplicate-free input BAM file 
#' @param output_folder The path of the output folder 
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#'
#' @export
sc_atac_create_fragments <- function(inbam, output_folder = NULL) {
    # Check if output directory exists
    if(is.null(output_folder)) {
        output_folder <- file.path(getwd(), "scPipe-atac-output")
    }

    if (!dir.exists(output_folder)){
        dir.create(output_folder,recursive=TRUE)
        cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
    }

    # Create log folder/file
    log_and_stats_folder <- file.path(output_folder, "scPipe_atac_stats")
    dir.create(log_and_stats_folder, showWarnings = FALSE)
    log_file <- file.path(log_and_stats_folder, "log_file.txt")
    if(!file.exists(log_file)) file.create(log_file)
    cat(
      paste0(
        "sc_atac_create_fragments starts at ",
        as.character(Sys.time()),
        "\n"
      ),
      file = log_file, append = TRUE)
    
    output <- file.path(output_folder, "fragments.bed")
    
    # Need to use sinto to generate fragment file
    proc <- basiliskStart(scPipe_env)
    on.exit(basiliskStop(proc))

    basiliskRun(proc, function(inbam, out) {
        sin <- reticulate::import("sinto")
        #sin$fragments(inbam, out)
        system2("sinto", c("fragments", "-b", inbam, "-f", out))
    }, inbam=inbam, out=output) 
    
    cat(
      paste0(
        "sc_atac_create_fragments finishes at ",
        as.character(Sys.time()),
        "\n\n"
      ),
      file = log_file, append = TRUE)
  
}