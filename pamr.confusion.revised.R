pamr.confusion.revised = function(fit, threshold) {
  
  if(class(fit$yhat) == "data.frame"){
    predicted = fit$yhat[,1]
  }else{
    predicted = fit$yhat  
  }
  
  
  if(!is.null(fit$y)){    
    true <- fit$y[fit$sample.subset]
    tt <- table(true, predicted)
  }
  else{
    true <- fit$proby[fit$sample.subset,]
    ytemp<- apply(true,1,which.is.max)
    temp <- c(predicted,names(table(ytemp)))
    nams <- names(table(temp))
     
    Yhat <- model.matrix( ~ factor(temp) - 1, data = list(y = temp))
    Yhat <- Yhat[1:length(predicted),]
    tt <- matrix(NA,nrow=length(fit$prior),ncol=length(fit$prior))
  
    for(i in 1:length(fit$prior)){
      for(j in 1:length(fit$prior)){
        tt[i,j] <- sum(true[,i]*Yhat[,j])
      }
    }
    dimnames(tt) <- list(names(table(ytemp)),nams)
  }
  tt1 <- tt
  diag(tt1) <- 0
  tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
  dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
  
  ttlist = list(tt = tt, threshold = fit$threshold)
  
}
