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
  
  
  log_and_stats_folder = paste0(output_folder, "/log_and_stats/")
  dir.create(log_and_stats_folder, showWarnings = F)
  log_file = paste0(log_and_stats_folder, "log_file.txt")
  stats_file = paste0(log_and_stats_folder, "stats_file_align.txt")
  if(!file.exists(log_file)) file.create(log_file)
  # file.create(stats_file)
  
  cat(
    paste0(
      "sc_atac_tagging starts at ",
      as.character(Sys.time()),
      "\n"
    ), 
    file = log_file, append = T)
  
  
  
  
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

  cat(
    paste0(
      "sc_atac_tagging finishes at ",
      as.character(Sys.time()),
      "\n\n"
    ), 
    file = log_file, append = T)
  
  }
