.onAttach <- function(...) {
    # packageStartupMessage("Some startup message")
    

    
    # Check for platform
    if(.Platform$OS.type != "unix") {
        packageStartupMessage("Windows platform detected, sc_atac_peak_calling() and sc_atac_remove_duplicates() function may not operate. Please call peaks and remove duplicate reads outside the package.")  
    } else {
        packageStartupMessage("Linux/MacOS detected! \n")
        
        # Check if samtools is installed
        packageStartupMessage("Checking if samtools is installed ... \n")
        tryCatch({
            system2("samtools", stdout = NULL, stderr = NULL)
            packageStartupMessage("samtools was located! \n")
        },
        
        warning = function(w) {
            packageStartupMessage("samtools was not located. Please specify the path of it when running sc_atac_remove_duplicates()")
        }
        )
    }
}