# .onLoad <- function(libname, pkgname) 
# {
#   library.dynam("mclust", pkgname, libname)
# }

mclustStartupMessage <- function()
{
# Startup message obtained as 
# > figlet -f slant MCLUST
  msg <- c(paste0(
"    __  ___________    __  _____________
   /  |/  / ____/ /   / / / / ___/_  __/
  / /|_/ / /   / /   / / / /\\__ \\ / /   
 / /  / / /___/ /___/ /_/ /___/ // /    
/_/  /_/\\____/_____/\\____//____//_/    version ", 
packageVersion("mclust")),
"\nType 'citation(\"mclust\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # unlock .mclust variable allowing its modification
  unlockBinding(".mclust", asNamespace("mclust")) 
  # startup message
  msg <- mclustStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'mclust' version", packageVersion("mclust"))
  packageStartupMessage(msg)      
  invisible()
}
