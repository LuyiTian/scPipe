#' @importFrom basilisk BasiliskEnvironment
scPipe_env <- BasiliskEnvironment(envname="scPipe_env", 
    pkgname="scPipe",
    packages=c("python==3.9.0"),
    pip=c("sinto==0.7.2.2", "MACS3==3.0.0a6"))