#' @name sc_atac_remove_duplicates
#'
#' @title Removing PCR duplicates using samtools
#' @description Takes in a BAM file and removes the PCR duplicates using the samtools markdup function. Requires samtools 1.10 or newer. 
#' 
#' @param inbam The tagged, sorted and duplicate-free input BAM file 
#' @param samtools_path The path of the samtools executable (if a custom installation is to be specified)
#' @param output_folder The path of the output folder 
#'
#' @export
#' 

sc_atac_remove_duplicates <- function(inbam, 
                                      samtools_path = NULL,
                                      output_folder = NULL) {
  
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
      message("samtools was located")
      return(TRUE)
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
  
  if (samtools.installed)
    tryCatch(
      {
        # Check if file exists
        if (!file.exists(inbam)) {
          stop("Couldn't locate the supplied BAM file")
        }
        
        # Check for validity of file
        system2(samtools, c("quickcheck", inbam), stderr = NULL, stdout = NULL)
        
        
        # Check if output directory exists
        if (is.null(output_folder)) {
          output_folder = file.path(getwd(), "scPipe-atac-output")
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
            "sc_atac_remove_duplicates starts at ",
            as.character(Sys.time()),
            "\n"
          ),
          file = log_file, append = TRUE)
        
        
        inbam.name <- substr(inbam, 0, nchar(inbam)-4)
        cat("Sorting the BAM file by name\n")
        system2(samtools, c("sort", "-n", "-o", paste(inbam.name, "namesorted.bam", sep="_"), inbam))
        cat("Running fixmate on the sorted BAM File\n")
        system2(samtools, c("fixmate", "-m", paste(inbam.name, "namesorted.bam", sep="_"), paste(inbam.name, "fixmate.bam", sep="_")))
        cat("Sorting the BAM file by coordinates\n")
        # Remove existing copy if one exists
        output.bam <- paste0(output_folder, "/", basename(inbam.name), "_markdup.bam")
        if (file.exists(output.bam)) {
          system2("rm", output.bam)
        }
        
        Rsamtools::sortBam(paste(inbam.name, "fixmate.bam", sep="_"), paste(inbam.name, "positionsort", sep="_"), indexDestination = TRUE)

        # Note: the output bam file is originally created in the same directory as the input bam file
        
        cat("Running samtools markdup now")
        # Check if the version of samtools is 1.10 or greater (to have the stats functionality)
        version.text <- strsplit(strsplit(system2(samtools, "--version", stdout=TRUE)[1], " ")[[1]][2], "\\.")[[1]]
        if (as.numeric(version.text[1]) < 1 || (as.numeric(version.text[1]) >= 1 && as.numeric(version.text[2]) < 10)) {
          cat("Version of samtools isn't greater or equal to 1.10. Can't produce duplicate removal stats file.\n")
          system2(samtools, c("markdup", "-s", "-r", paste(inbam.name, "positionsort.bam", sep="_"), paste(inbam.name, "markdup.bam", sep="_")))
        } else {
          cat("Detected samtools with version greater or equal to 1.10. Producing duplicate removal stats file.\n")
          system2(samtools, c("markdup", "-s", "-f", file.path(log_and_stats_folder, "duplicate_removal_stats.txt"), "-r", paste(inbam.name, "positionsort.bam", sep="_"), paste(inbam.name, "markdup.bam", sep="_")))
        }

        cat("Indexing the BAM file\n")
        Rsamtools::indexBam(paste(inbam.name, "markdup.bam", sep="_"))
        cat("step 6\n")
        system2("rm", paste(inbam.name, "namesorted.bam", sep="_"))
        system2("rm", paste(inbam.name, "positionsort.bam", sep="_"))
        system2("rm", paste(inbam.name, "fixmate.bam", sep="_"))
        
        if (file.exists(paste(inbam.name, "markdup.bam", sep="_"))) {
          message(paste("The output BAM file was sent to", output_folder))
          
          # If the new file is already in the destination folder, don't need to move it!
          if (paste(inbam.name, "markdup.bam", sep="_") != output.bam)
            system2("mv", c("--force", paste(inbam.name, "markdup.bam", sep="_"), output_folder))
          
        } else {
          message("Couldn't remove duplicates from the input BAM file. Make sure it is a valid BAM file.")
        }
        cat(
          paste0(
            "sc_atac_remove_duplicates finishes at ",
            as.character(Sys.time()),
            "\n\n"
          ),
          file = log_file, append = TRUE)
      },
      warning = function(w) {w
        message(w)
      },
      
      error = function(e) {
        message(e)
      }
    )
}
