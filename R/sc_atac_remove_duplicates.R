#' @name sc_atac_remove_duplicates
#'
#' @title Removing PCR duplicates using samtools
#' @description Takes in a BAM file and removes the PCR duplicates using the samtools markdup function. Requires samtools 1.10 or newer for statistics to be generated. 
#' 
#' @param inbam The tagged, sorted and duplicate-free input BAM file 
#' @param samtools_path The path of the samtools executable (if a custom installation is to be specified)
#' @param output_folder The path of the output folder 
#'
#' @returns file path to a bam file created from samtools markdup
#' @importFrom Rsamtools sortBam
#' @export
#' 
sc_atac_remove_duplicates <- function(inbam, 
                                        samtools_path = NULL,
                                        output_folder = NULL) {
    
    successfully_removed_duplicates <- FALSE
    # Check if samtools is installed
    if (!is.null(samtools_path)) {
        samtools <- samtools_path
    } else {
        samtools <- "samtools"
    }
    
    # Check if samtools is installed
    samtools.installed <- tryCatch(
        {
        system2(samtools, stdout = NULL, stderr = NULL)
        message("samtools was located!")
        TRUE
        },
        
        warning = function(w) {
            if (is.null(samtools_path)) {
                message("samtools was not located, so can't remove duplicates. Please make sure it is installed.")
            } else {
                message("samtools was not located via the specified path. Please make sure it is correct.")
            }
            return(FALSE)
        }
    )

    if (isTRUE(samtools.installed)) {
        return(tryCatch({
            # Check if file exists
            if (!file.exists(inbam)) {
                stop("Couldn't locate the supplied BAM file")
            }
            
            # Check for validity of file
            system2(samtools, c("quickcheck", inbam), stderr = NULL, stdout = NULL)
            
            
            # Check if output directory exists
            if (is.null(output_folder)) {
                output_folder <- file.path(getwd(), "scPipe-atac-output")
            }
            
            if (!dir.exists(output_folder)){
                dir.create(output_folder,recursive=TRUE)
                message("Output Directory Does Not Exist. Created Directory: ", output_folder)
            }
            
            # Create log folder/file
            log_and_stats_folder <- file.path(output_folder, "scPipe_atac_stats")
            dir.create(log_and_stats_folder, showWarnings = FALSE)
            log_file <- file.path(log_and_stats_folder, "log_file.txt")
            if(!file.exists(log_file)) file.create(log_file)
            write(c(
                "sc_atac_remove_duplicates starts at ",
                as.character(Sys.time()),
                "\n"
                ),
                file = log_file, append = TRUE)
            
            
            inbam.name <- substr(inbam, 0, nchar(inbam)-4)
            message("Sorting the BAM file by name")
            system2(samtools, c("sort", "-n", "-o", paste(inbam.name, "namesorted.bam", sep="_"), inbam))
            message("Running fixmate on the sorted BAM File")
            system2(samtools, c("fixmate", "-m", paste(inbam.name, "namesorted.bam", sep="_"), paste(inbam.name, "fixmate.bam", sep="_")))
            message("Sorting the BAM file by coordinates")
            # Remove existing copy if one exists
            output.bam <- paste0(output_folder, "/", basename(inbam.name), "_markdup.bam")
            if (file.exists(output.bam)) {
                system2("rm", output.bam)
            }
            
            Rsamtools::sortBam(paste(inbam.name, "fixmate.bam", sep="_"), paste(inbam.name, "positionsort", sep="_"), indexDestination = TRUE)

            # Note: the output bam file is originally created in the same directory as the input bam file
            
            message("Running samtools markdup now")
            # Check if the version of samtools is 1.10 or greater (to have the stats functionality)
            version.text <- strsplit(strsplit(system2(samtools, "--version", stdout=TRUE)[1], " ")[[1]][2], "\\.")[[1]]
            if (as.numeric(version.text[1]) < 1 || (as.numeric(version.text[1]) >= 1 && as.numeric(version.text[2]) < 10)) {
                message("Version of samtools isn't greater or equal to 1.10. Can't produce duplicate removal stats file.")
                system2(samtools, c("markdup", "-s", "-r", paste(inbam.name, "positionsort.bam", sep="_"), paste(inbam.name, "markdup.bam", sep="_")))
            } else {
                message("Detected samtools with version greater or equal to 1.10. Producing duplicate removal stats file.")
                system2(samtools, c("markdup", "-s", "-f", file.path(log_and_stats_folder, "duplicate_removal_stats.txt"), "-r", paste(inbam.name, "positionsort.bam", sep="_"), paste(inbam.name, "markdup.bam", sep="_")))
            }

            message("Indexing the BAM file")
            Rsamtools::indexBam(paste(inbam.name, "markdup.bam", sep="_"))
            message("Cleaning up intermediary files")
            system2("rm", paste(inbam.name, "namesorted.bam", sep="_"))
            system2("rm", paste(inbam.name, "positionsort.bam", sep="_"))
            system2("rm", paste(inbam.name, "fixmate.bam", sep="_"))
            
            if (file.exists(paste(inbam.name, "markdup.bam", sep="_"))) {
                message("The output BAM file was sent to", output_folder)
                
                # If the new file is already in the destination folder, don't need to move it!
                if (paste(inbam.name, "markdup.bam", sep="_") != output.bam)
                    system2("mv", c("--force", paste(inbam.name, "markdup.bam", sep="_"), output_folder))
                successfully_removed_duplicates <- output.bam
            } else {
                message("Couldn't remove duplicates from the input BAM file. Make sure it is a valid BAM file.")
            }
            write(
            c(
                "sc_atac_remove_duplicates finishes at ",
                as.character(Sys.time()),
                "\n\n"
            ),
            file = log_file, append = TRUE)
            successfully_removed_duplicates
        },
        warning = function(w) {w
            message(w)
            return(FALSE)
        },
        
        error = function(e) {
            message(e)
            return(FALSE)
        }))
    } else {
        return(FALSE)
    }
}
