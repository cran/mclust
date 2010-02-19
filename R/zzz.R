.First.lib <- function(lib, pkg) {
cat("\nby using mclust, invoked on its own or through another package,\n")
cat("you accept the license agreement in the mclust LICENSE file\n")
cat("and at http://www.stat.washington.edu/mclust/license.txt\n\n")
  library.dynam("mclust", pkg, lib)
}
