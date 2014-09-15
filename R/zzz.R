.onAttach <- function(lib, pkg)
{
  unlockBinding(".mclust", asNamespace("mclust")) 
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  
  packageStartupMessage("Package 'mclust' version ", version)
  
  # figlet version obtained as 
  # > figlet -f slant mclust
  # packageStartupMessage("                   __           __
#   ____ ___  _____/ /_  _______/ /_
#  / __ `__ \\/ ___/ / / / / ___/ __/
# / / / / / / /__/ / /_/ (__  ) /_ 
#/_/ /_/ /_/\\___/_/\\__,_/____/\\__/

#Package version ", version)
  # > figlet -f slant MCLUST
#  packageStartupMessage("    __  ___________    __  _____________
#   /  |/  / ____/ /   / / / / ___/_  __/
#  / /|_/ / /   / /   / / / /\\__ \\ / /   
# / /  / / /___/ /___/ /_/ /___/ // /    
# /_/  /_/\\____/_____/\\____//____//_/    version ", version)

  # packageStartupMessage("----------------")
  # packageStartupMessage("By using mclust, invoked on its own or through another package,\nyou accept the license agreement in the mclust LICENSE file\nand at http://www.stat.washington.edu/mclust/license.txt\n")
  packageStartupMessage("Type 'citation(\"mclust\")' for citing this R package in publications.")
  invisible()
}



