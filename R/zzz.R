.onAttach <- function(...) {
  packageStartupMessage("Some startup message")
  
  if(.Platform$OS.type != "unix") {
    packageStartupMessage("Some functions may not be working for Windows")  
  } 
  
  
  
}