
###############################################
# Create the final report for scATAC-Seq data
###############################################

#' @name sc_atac_create_report
#' @title HTML report generation
#' @description Generates a HTML report using the output folder produced by the pipeline
#'
#' @param input_folder The path of the folder produced by the pipeline
#' @param output_folder The path of the output folder to store the HTML report in
#' @param organism A string indicating the name of the organism being analysed
#' @param sample_name A string indicating the name of the sample
#' @param feature_type A string indicating the type of the feature (`genome_bin` or `peak`)
#' @param n_barcode_subset if you require only to visualise stats for a sample of barcodes to improve processing time (integer)
#'
#' @returns the path of the output file
#' @export
sc_atac_create_report <- function(input_folder, 
                                    output_folder    = NULL, 
                                    organism         = NULL,
                                    sample_name      = NULL,
                                    feature_type     = NULL,
                                    n_barcode_subset = 500){
    
    if (!dir.exists(input_folder)){
        stop("The input folder could not be found at ", input_folder);
    }

    if (!requireNamespace("rmarkdown", quietly=TRUE)) {
        stop("Install 'rmarkdown' to use this function")
    }
    
    if (!dir.exists(file.path(input_folder, "scPipe_atac_stats"))) {
        stop("The input stats folder could not be found at ", 
            file.path(input_folder, "scPipe_atac_stats"), 
            ". Either run the pipeline to generate it or create a `scPipe_atac_stats` folder that contains the required files.");
    }
    
    available_organisms <- c("hg19",
                            "hg38",
                            "mm10")


    if(is.null(output_folder)) {
        output_folder <- file.path(input_folder, "scPipe_atac_stats")
    }
    
    if (!dir.exists(output_folder)){
        message("Output folder could not be found at " , output_folder,  "Creating it now...");
        dir.create(output_folder, recursive=TRUE)
    }
    
    report.file <- system.file("extdata/rmd_report_skeleton.Rmd", package = "scPipe") 
    if (report.file == "") # Used if installed via the repo
        report.file <- system.file("inst/extdata/rmd_report_skeleton.Rmd", package = "scPipe", mustWork = TRUE)
    
    output_file <- file.path(output_folder, "scPipe_atac_report.html")
    rmarkdown::render(
        input  = report.file,
        params = list(
        input_folder     = input_folder,
        organism         = organism,
        sample           = sample_name,
        feature_type     = feature_type,
        n_barcode_subset = n_barcode_subset
        ),
        output_file = output_file 
    )

    return(output_file)
}
