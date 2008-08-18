.First.lib <- function(lib, pkg) {
cat("by using mclust, you accept the license agreement in the LICENSE file\n")
cat("and at http://www.stat.washington.edu/mclust/license.txt\n")
  library.dynam("mclust", pkg, lib)
}
