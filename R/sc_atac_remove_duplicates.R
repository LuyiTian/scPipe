#' sc_atac_remove_duplicates()
#'
#' @return 
#'
#' @examples
#' \dontrun{
#' 
#' 
#' }
#'
#' @export
#' 

sc_atac_remove_duplicates <- function(inbam, samtools_path = NULL, output_folder = ""){
  
  # Check if samtools is installed
  samtools.installed <<- TRUE
  
  
  if (!is.null(samtools_path)) {
    samtools <- samtools_path
  } else {
    samtools <- "samtools"
  }
  
  # Check if samtools is installed
  tryCatch(
    {
      system2(samtools, stdout = NULL, stderr = NULL)
      message("samtools was located")
    },
    
    warning = function(w) {
      samtools.installed <<- FALSE
      if (is.null(samtools_path)) {
        message("samtools was not located, so can't remove duplicates. Please make sure it is installed.")
      } else {
        message("samtools was not located via the specified path. Please make sure it is correct.")
      }
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
        if(output_folder == ''){
          output_folder = file.path(getwd(), "scPipe-atac-output")
        }
        
        if (!dir.exists(output_folder)){
          dir.create(output_folder,recursive=TRUE)
          cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
        }
        
        inbam.name <- substr(inbam, 0, nchar(inbam)-4)
        
        system2(samtools, c("collate", "-o", paste(inbam.name, "namecollate.bam", sep="_"), inbam))
        system2(samtools, c("fixmate", "-m", paste(inbam.name, "namecollate.bam", sep="_"), paste(inbam.name, "fixmate.bam", sep="_")))
        
        # Remove existing copy if one exists
        output.bam <- paste0(output_folder, "/", basename(inbam.name), "_markdup.bam")
        if (file.exists(output.bam)) {
          system2("rm", output.bam)
        }
        
        system2(samtools, c("sort", "-o", paste(inbam.name, "positionsort.bam", sep="_"), paste(inbam.name, "fixmate.bam", sep="_")))
        
        # Note: the output bam file is originally created in the same directory as the input bam file
        
        system2(samtools, c("markdup", "-s", "-r", paste(inbam.name, "positionsort.bam", sep="_"), paste(inbam.name, "markdup.bam", sep="_")))
        Rsamtools::indexBam(paste(inbam.name, "markdup.bam", sep="_"))
        
        system2("rm", paste(inbam.name, "namecollate.bam", sep="_"))
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
        
      },
      warning = function(w) {w
        message(w)
      },
      
      error = function(e) {
        message(e)
      }
    )
}
