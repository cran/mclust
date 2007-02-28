.First.lib <- function(lib, pkg) {
  cat("use of mclust requires a license agreement\n")
  cat("see http://www.stat.washington.edu/mclust/license.txt\n")
  library.dynam("mclust", pkg, lib)
}
