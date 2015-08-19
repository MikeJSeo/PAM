pamr.plotTrainError = function(fit){
  ylength = length(fit$y)
  par(cex = 0.7)
  plot(fit$threshold, fit$error/ylength, ylim = c(0, 0.9), xlab = "Threshold", ylab = "Training Error", col = 4)
  axis(3, at = fit$threshold, labels = fit$nonzero)  
}

pamr.plotTestError = function(fit, testData){
  
  x = testData$x
  y = testData$y
  
  par(cex = 0.7)
  threshold = fit$threshold
  ylength = sum(!is.na(y))
  errorrate = vector(mode = "logical", length = length(threshold))
  
  for(i in 1: length(threshold)){
    
    predict = pamr.predict(fit, x, threshold = threshold[i], type = "class")
    errorrate[i] = sum(y != predict, na.rm = TRUE) / ylength
  }
  print(fit$threshold)
  print(errorrate)
  plot(fit$threshold, errorrate, ylim = c(0, 0.9), xlab = "Threshold", ylab = "Test Error", col = 4)
  axis(3, at = fit$threshold, labels = fit$nonzero)  
}
