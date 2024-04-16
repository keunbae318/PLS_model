
Normalization100 <- function(X) {
  if (!any(class(X) %in% c("matrix", "data.frame"))) {
    stop("X must be a matrix or optionally a data.frame")
  }
  area <- apply(X, 1, sum)
  X <- sweep(X, 1, area, "/")
  X <- sweep(X, 1, 100, "*")
  as.matrix(X)
}

