rmOutlier <- function(x, s=3) {
  ans <- x[x < mean(x) + 3*sd(x) & ref > mean(ref) - 3*sd(ref)]
  ans
}