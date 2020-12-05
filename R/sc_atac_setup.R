#' Set up the background forscATAC-Seq pre-processing
#' sc_atac_setup
#' @return NULL
#'
#' @examples
#' \dontrun{
#' setpaths(macs2path,bedToolsPath)
#' 
#' }
#'
#' @export
#'


#********************** Packages (install if not found)

list_of_packages <- c("tools", "Rsubread", "tidyr", "locStra", "grid", "dplyr", "Matrix", "rtracklayer", "GenomicAlignments", "GenomicRanges", "GenomicFeatures", "Rsamtools", "SingleCellExperiment", "seqinr")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if("Rsubread" %in% req_packages) BiocManager::install("Rsubread")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

library("tools")
library("Rsubread")
library("dplyr")
library("tidyr")
library("locStra")
library("Matrix")
library("rtracklayer")
library("GenomicAlignments")
library("GenomicRanges")
library("GenomicFeatures")
library("Rsamtools")
library("SingleCellExperiment")
library("seqinr")

# Set the paths for subsequent usees ------
# need to change this to work only if the OS is not windows...
set_paths = function(macs2path,bedToolsPath){
  
  old_path <- Sys.getenv("PATH")
  new_path = paste(old_path, macs2Path, sep = ":")
  new_path = paste(old_path, bedToolsPath, sep = ":")
  
  Sys.setenv(PATH = new_path)
  
}