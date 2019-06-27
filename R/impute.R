imputeData <- function(data, categorical = NULL, 
                       seed = NULL, 
                       verbose = interactive()) 
{
  if(!requireNamespace("mix", quietly = TRUE)) 
     stop("imputeData function require 'mix' package to be installed!")
  
  fac <- apply(data, 2, is.factor)
  if(is.null(categorical)) 
    { categorical <- fac }
  else 
    { if(any(!categorical & fac)) 
        { stop("data has a factor that is not designated as categorical") }
      if(any(categorical | !fac)) 
        { warning("a categorical is not designated as a factor")
         for(i in which(categorical | !fac)) 
             data[[i]] <- as.factor(data[[i]])
      }
  }
  
  # remove categorical variables and add a dummy variable
  if(nocat <- !any(categorical)) 
    { data <- cbind(as.factor(1), data)
      categorical <- c(TRUE, categorical)
  }
  
  ord <- c(which(categorical), which(!categorical))
  
  # do the imputations
  s <- mix::prelim.mix(data[,ord], p = sum(categorical))
  if(is.null(seed)) 
    seed <- runif(1, min = .Machine$integer.max/1024,
                     max = .Machine$integer.max)
  # find ML estimate
  thetahat <- mix::em.mix(s, showits = verbose)
  # set random number generator seed
  mix::rngseed(seed) 
  # data augmentation from posterior
  newtheta <- mix::da.mix(s, thetahat, steps = 100, 
                          showits = verbose)
  # impute under newtheta
  dataImp <- mix::imp.mix(s, newtheta) 
  # there is a bug, so it needs to refix the seed and impute again
  mix::rngseed(seed) 
  dataImp <- mix::imp.mix(s, newtheta) 
  
  if(nocat) dataImp[,-1] else dataImp[,order(ord)]
}

imputePairs <- function(data, dataImp, 
                        symbols = c(1, 16), 
                        colors = c("black", "red"),
                        labels, panel = points, ...,  
                        lower.panel = panel, 
                        upper.panel = panel, 
                        diag.panel = NULL, 
                        text.panel = textPanel, 
                        label.pos = 0.5 + has.diag/3, 
                        cex.labels = NULL, font.labels = 1, 
                        row1attop = TRUE, gap = 0.2) 
{
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) 
    text(x, y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, oma, ...) 
  {
    if (side%%2 == 1) 
      Axis(x, side = side, xpd = NA, ...)
    else Axis(y, side = side, xpd = NA, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(data)) {
    data <- as.data.frame(data)
    for (i in seq_along(names(data))) {
      if (is.factor(data[[i]]) || is.logical(data[[i]])) 
        data[[i]] <- as.numeric(data[[i]])
      if (!is.numeric(unclass(data[[i]]))) 
        stop("non-numeric argument to 'pairs'")
    }
  }
  else if (!is.numeric(data)) 
    stop("non-numeric argument to 'pairs'")
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(data)
  if (nc < 2) 
    stop("only one column in the argument to 'pairs'")
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(data)
    if (is.null(labels)) 
      labels <- paste("var", 1:nc)
  }
  else if (is.null(labels)) 
    has.labs <- FALSE
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  else NULL
  main <- if ("main" %in% nmdots) 
    dots$main
  else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main)) 
      oma[3] <- 6
  }
  opar <- par(mfrow = c(nc, nc), mar = rep(gap/2, 4), oma = oma)
  on.exit(par(opar))
  for (i in if (row1attop) 1:nc else nc:1) 
    for (j in 1:nc) {
      localPlot(dataImp[, j], dataImp[, i], xlab = "", ylab = "", 
                axes = FALSE, type = "n", ...)
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        box()
        if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
          localAxis(1 + 2 * row1attop, dataImp[, j], dataImp[, i], 
                    ...)
        if (i == nc && (j%%2 || !has.upper || !has.lower)) 
          localAxis(3 - 2 * row1attop, dataImp[, j], dataImp[, i], 
                    ...)
        if (j == 1 && (!(i%%2) || !has.upper || !has.lower)) 
          localAxis(2, dataImp[, j], dataImp[, i], ...)
        if (j == nc && (i%%2 || !has.upper || !has.lower)) 
          localAxis(4, dataImp[, j], dataImp[, i], ...)
        mfg <- par("mfg")
        if (i == j) {
          if (has.diag) 
            localDiagPanel(as.vector(dataImp[, i]), ...)
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                       font = font.labels)
          }
        }
        else if (i < j) { 
          classification <- as.numeric(apply(data[,c(i,j)], 1, 
                                             function(x) any(is.na(x)))) + 1
          localLowerPanel(as.vector(dataImp[, j]), as.vector(dataImp[,i]), 
                          pch = symbols[classification], 
                          col = colors[classification], ...)
        }
        else {
          classification <- as.numeric(apply(data[,c(i,j)], 1, 
                                             function(x) any(is.na(x)))) + 1
          localUpperPanel(as.vector(dataImp[, j]), 
                          as.vector(dataImp[, i]), 
                          pch = symbols[classification], 
                          col = colors[classification], ...)
        }
        if (any(par("mfg") != mfg)) 
          stop("the 'panel' function made a new plot")
      }
      else par(new = FALSE)
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}

# LS: old to be removed
# matchCluster <- function(group, cluster)
# {
#   if(length(group) != length(cluster))
#     stop("arguments must be vector of the same length")
#   group <- as.factor(group)
#   cluster <- as.factor(cluster)  
#   tab <- table(group,cluster)
#   j <- apply(tab,2,which.max)
#   cluster <- factor(cluster, labels = levels(group)[j])
#   cluster <- as.character(cluster)
#   group <- as.character(group)
#   misclassified <- !(cluster == group)
#   out <- list(cluster = cluster, misclassified = misclassified, ord = j)
#   return(out)
# }

matchCluster <- function(group, cluster)
{
  if(length(group) != length(cluster))
    stop("arguments must be vector of the same length")
  group <- as.factor(group)
  cluster <- as.factor(cluster)  
  map <- mapClass(as.numeric(group), as.numeric(cluster))
  map1 <- unlist(map[[1]]); names(map1) <- NULL
  map2 <- unlist(map[[2]]); names(map2) <- NULL
  cl <- cluster
  levels(cl) <- map2
  cl <- as.character(levels(cl)[as.numeric(cl)])
  cl <- as.character(cl)
  group <- as.character(group)
  misclassified <- !(cluster == group)
  out <- list(cluster = cl, misclassified = misclassified, ord = map1)
  return(out)
}

majorityVote <- function(x)
{
  # local function to find the maximum position in a vector, 
  # breaking ties at random
  whichMax <- function (x) 
  {
    m <- seq_along(x)[x == max(x, na.rm = TRUE)]
    if(length(m) > 1) sample(m, size = 1) else m
  }
  x <- as.vector(x)
  tab <- table(x)
  m <- whichMax(tab)
  out <- list(table = tab, ind = m, majority = names(tab)[m])
  return(out)
}
