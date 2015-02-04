pamr.plotpredprob <- function(fit, testData, threshold) {
  
  x = testData$x
  y = testData$y
  samplelabel = testData$sample
  
  pp = pamr.predict(fit, x, threshold = threshold, type = "posterior")
  classPrediction = pamr.predict(fit, x, threshold = threshold, type = "class")
  o <- order(y)
  y <- factor(y[o])

  ss = samplelabel
  if(!is.null(ss)) {
    ss <- ss[o]
  }
  ppp <- pp[o,]
  n <- nrow(ppp)
  nc <- length(unique(classPrediction))

  plot(1:n, ppp[, 2], type = "n", xlab = "Sample", ylab = "Predicted Test Probabilities", ylim = c(0, 1.5), axes = FALSE)
  axis(1)
  axis(2, at=seq(0, 1.2, by=0.2), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0", ""))

  for(j in 1:nc) {
    points(1:n, ppp[,j], col = j +1)
    abline(v = cumsum(table(y))[j] + 0.5, lty = 2)
    text( c(0, cumsum(table(y)))[j] + 0.5 * table(y)[j], 1.1, label = levels(y)[j], col = j +1)
  }
  abline(h = 1)
  if(!is.null(ss)) {
    text(1:length(ss), 1.3, labels = ss, srt = 90, cex = 0.7)
  }

}
