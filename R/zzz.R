.onAttach <- function(lib, pkg)
{
  unlockBinding(".mclust", asNamespace("mclust")) 
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage("Package 'mclust' version ", version)
  # packageStartupMessage("----------------")
  # packageStartupMessage("By using mclust, invoked on its own or through another package,\nyou accept the license agreement in the mclust LICENSE file\nand at http://www.stat.washington.edu/mclust/license.txt\n")
  # packageStartupMessage("Type 'citation(\"mclust\")' for citing this R package in publications.")
  invisible()
}



  