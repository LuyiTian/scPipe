#' sc_atac_bam_tagging()
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

sc_atac_bam_tagging = function(inbam, outfolder = "",
                               bam_tags = list(bc="CB", mb="OX"),
                               nthreads=1) {
  
  if (any(!file.exists(inbam))) {
    stop("At least one input bam file should be present")
  } else {
    cat("sc_atac_bam_tagging1\n")
    inbam = path.expand(inbam)
  }
  
  
  if(outfolder == ""){
    fileNameWithoutExtension <- strsplit(inbam, "\\.")[[1]][1]
    outbam                   <- paste(fileNameWithoutExtension, "_tagged.bam", sep = "")
    outsortedbam             <- paste(fileNameWithoutExtension, "_tagged_sorted", sep = "")
  } else{
    if(!dir.exists(outfolder)){
      cat(outfolder, "does not exist.\nCreating folder...")
      dir.create(outfolder)
      cat("Created.\n")
    }
    
    fileNameWithoutExtension <- basename(strsplit(inbam, "\\.")[[1]][1])
    outbam                   <- paste(outfolder, "/", fileNameWithoutExtension, "_tagged.bam", sep = "")
    outsortedbam             <- paste(outfolder, "/", fileNameWithoutExtension, "_tagged_sorted", sep = "")
  }
  
  outbam <- path.expand(outbam)
  if(!file.exists(outbam)){
    file.create(outbam)
  }
  
  # if(output_folder == ''){
  #   output_folder <- file.path(getwd(), "scPipe-atac-output")
  # }
  # 
  # if (!dir.exists(output_folder)){
  #   dir.create(output_folder,recursive=TRUE)
  #   cat("Output Directory Does Not Exist. Created Directory: ", output_folder, "\n")
  # }
  
  
  rcpp_sc_atac_bam_tagging(inbam, outbam, bam_tags$bc, bam_tags$mb,nthreads)
  
  cat("Tagged BAM file is located in: \n")
  cat(outbam)
  cat("\n")
  
  Rsamtools::sortBam(outbam, outsortedbam,indexDestination = TRUE)
  Rsamtools::indexBam(paste0(outsortedbam, ".bam"))
  
  cat("Sorted & indexed tagged BAM file is located in: \n")
  cat(paste0(outsortedbam, ".bam"))
  cat("\n")
  
  # generate the fragment file for the BAM file 
  # need bedtools v2.26.0 or later
  # system2("bedtools", c("bamToBed", "i", outsortedbam), "|", "awk", c(-F"#" '{print $1"\t"$2}'), stdout = paste(output_folder,"/fragments.bed",sep = ""))
}
