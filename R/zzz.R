.onAttach <- function(lib, pkg){
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage("Package 'mclust' version ", version)
}
