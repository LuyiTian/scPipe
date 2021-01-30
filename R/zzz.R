.onAttach <- function(...) {
  # packageStartupMessage("Some startup message")
  
  if(.Platform$OS.type != "unix") {
    packageStartupMessage("Windows platform detected, sc_atac_peak_calling() function may not operate. Please call peaks outside the package.")  
  } else{
    packageStartupMessage("Linux/MacOS detected, please make sure the latest MACS2 is in path if need to run sc_atac_peak_calling()")
  }
  
  
  
}