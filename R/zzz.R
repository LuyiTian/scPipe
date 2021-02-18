.onAttach <- function(...) {
  # packageStartupMessage("Some startup message")
  
  
  # Check if required packages are installed
  list_of_packages <- c("tools", "tidyr", "locStra", "grid", "dplyr", "Matrix", "seqinr")
  bioconductor_packages <- c("BiocParallel", "Rsubread", "rtracklayer", "GenomicAlignments", "GenomicRanges", "GenomicFeatures", "Rsamtools", "SingleCellExperiment")
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")
  
  req_bioc_packages <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
  req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
  
  if(length(req_bioc_packages)) BiocManager::install(req_bioc_packages)
  if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")
  
  for (package in list_of_packages)
    library(package, character.only = TRUE)
  for (package in bioconductor_packages)
    library(package, character.only = TRUE)
  
  if (!"MACSr" %in% installed.packages()[,"Package"]) {
    if (paste0(R.version$major, ".", R.version$minor) > "4.1.0") {
      devtools::install_github("macs3-project/MACSr")
    } else {
      cat("R version 4.1.0 or greater is required for MACSr which is needed for sc_atac_call_peaks.\n")
    }
  }
  # Check for platform
  if(.Platform$OS.type != "unix") {
    packageStartupMessage("Windows platform detected, sc_atac_peak_calling() and sc_atac_remove_duplicates function may not operate. Please call peaks and remove duplicate reads outside the package.")  
  } else {
    packageStartupMessage("Linux/MacOS detected")
    
    # Check if samtools is installed
    message("Checking if samtools installed")
    tryCatch(
      {
        system2("samtools", stdout = NULL, stderr = NULL)
        message("samtools was located")
      },
      
      warning = function(w) {
        message("samtools was not located. Please make sure it is installed to run sc_atac_remove_duplicates")
      }
    )
  }
  
  
  
}