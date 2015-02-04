options(shiny.maxRequestSize=10000*1024^2)

library(openxlsx)
library(pamr)
library(impute)
library(superpc)

source("pamr.confusion.revised.R")
source("pamr.plotcv.revised.R")
source("pamr.plotfdr.revised.R")
source("pamr.plotpredprob.R")
source("pamr.plotError.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  ############# Get data
  
  getData = reactive({
    
    objFile = input$iFile
    classLabel = input$classLabel
    sampleLabel = input$sampleLabel
    expressionStart = input$expressionStart
    batchLabel = input$batchLabel
    
    if(!is.null(objFile) && !is.na(classLabel) && !is.na(expressionStart)){
      
      dat = read.xlsx(objFile$datapath, 1, colNames = FALSE)   
      emptySpace = 1:(expressionStart-1)
      
      geneid = dat[-emptySpace,1]
      genenames = dat[-emptySpace,2]   
      
      if(!is.na(sampleLabel)){
        sample = as.matrix(dat[sampleLabel, c(-1,-2)])  
      }
      
      
      y = as.matrix(dat[classLabel,c(-1,-2)])
      y = factor(y)        

      x = as.matrix(dat[-emptySpace, c(-1,-2)])
      class(x) = "numeric"
      
      imputedX = NULL
      if(sum(is.na(x)) > 0){  
        imputedummy = impute.knn(x, k = input$numberOfNeighbors)
        x = imputedummy$data
        imputedX = cbind(geneid, genenames, x)
      }
      
      transformed = NULL
      if(!is.null(input$center) || !is.null(input$scale) || !is.null(input$cuberoot)){
        x = scale(x, center = input$center, scale = input$scale)  
        if(input$cuberoot){
          x = sign(x) * abs(x)^(1/3)
        }
        transformed = cbind(geneid, genenames, x)
      }
      
      if(!is.na(batchLabel)){
        batchLabels = as.matrix(dat[batchLabel, c(-1,-2)])
        mydata <- list(x=x,y=y,batchlabels=factor(batchLabels))
        mydata2 <- pamr.batchadjust(mydata)      
        x = mydata2$x  
      }
      
      data = list(x = x, y = y, genenames = genenames, geneid = geneid, dat = dat, sample = sample, transformed = transformed, imputedX = imputedX)
      data
    }
  })


  getTestData = reactive({
    
    objFile = input$testFile
    classLabel = input$classLabel
    sampleLabel= input$sampleLabel
    expressionStart = input$expressionStart
    
    if(!is.null(objFile) && !is.na(classLabel) && !is.na(expressionStart)){
      
      dat = read.xlsx(objFile$datapath, 1, colNames = FALSE)   
      emptySpace = 1:(expressionStart-1)
      
      y = as.matrix(dat[classLabel,c(-1,-2)])
      
      x = as.matrix(dat[-emptySpace, c(-1,-2)])
      class(x) = "numeric"
      
      if(!is.na(sampleLabel)){
        sample = as.matrix(dat[sampleLabel, c(-1,-2)])  
      }
      
      list(x = x, y = y, sample = sample)
    }
  })
  

  getPamrPredict = reactive({
    testData = getTestData()
    pamrTrain = getPamrTrain()
    if(!is.null(testData) && !is.na(input$threshold)){
      posterior = pamr.predict(pamrTrain, testData$x, threshold = input$threshold, type = "posterior")
      posterior = t(posterior)     
      Predicted = pamr.predict(pamrTrain, testData$x, threshold = input$threshold, type= "class")
      Predicted = as.vector(Predicted)
      
      newposterior = rbind(testData$y, Predicted, posterior)
      rownames(newposterior) = c("Class Labels", "Predicted Class Labels", rownames(posterior))
      
      if(!is.na(input$sampleLabel)){
        colnames(newposterior) = testData$sample
      }
      newposterior
    }
  })

  
  getPamrTrain = reactive({
    dat = getData()
    if(!is.null(dat)){
      data = list(x = dat$x, y = dat$y, genenames = dat$genenames, geneid = dat$geneid)    
      pamr.train(data, offset.percent = input$s0percentile, sign.contrast = input$sign.contrast, prior = getPrior())
    }
  })
  
  getPamrCv = reactive({
    dat = getData()
    if(!is.null(dat)){
      set.seed(input$randomSeed)
      fit = getPamrTrain()
      pamr.cv(fit, dat)
    }
  })

  
  getPamrConfusionTrain = reactive({
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain) && !is.na(input$threshold)){
      
      maxvalue = round(pamrTrain$threshold[length(pamrTrain$threshold)], 5)
      
      if(input$threshold <= maxvalue){
        pamrConfusion = pamr.confusion.revised(pamrTrain, threshold = input$threshold)
        pamrConfusion
      }
    }
  })
  
  getPamrConfusion = reactive({
    pamrCv = getPamrCv()
    if(!is.null(pamrCv) && !is.na(input$threshold)){
      
      pamrTrain = getPamrTrain()
      maxvalue = round(pamrTrain$threshold[length(pamrTrain$threshold)], 5)
            
      if(input$threshold <= maxvalue){
        pamrConfusion = pamr.confusion.revised(pamrCv, threshold = input$threshold)
        pamrConfusion
      }
    }
  })
  
  getPamrListGenes = reactive({
    
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain) && !is.na(input$threshold)){
      data = getData()
      pamr.listgenes(pamrTrain, data, threshold = input$threshold, genenames = TRUE)
    }
  })
  
  getFDR = reactive({
    data = getData()
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain)){
      fdr.obj = pamr.fdr(pamrTrain, data)
    }
  })
  
  getPrior = reactive({
    
    data = getData()
    if(!is.null(data)){
      y = data$y
      lengthy = length(unique(y)) 
      prior = NULL
      
      if(input$classPrior == "uniformprior"){
        prior = rep(1/lengthy, lengthy)
      }    
      if(input$classPrior == "customprior"){
        prior = as.numeric(unlist(strsplit(input$custom, split=",")))
        if(sum(prior) != 1 || length(prior) != lengthy){
          prior = NULL
        }
      }
      prior           
    }
  })
  
  
  getSettings = reactive({
    pamrTrain = getPamrTrain()
    pamrConfusion = getPamrConfusion()
    data = getData()
    
    if(!is.null(pamrTrain) && !is.null(pamrConfusion)){
      settingMatrix = matrix(NA, nrow = 5, ncol = 2)
      colnames(settingMatrix) = c("Questions", "Values")
      settingMatrix[1,1] = "Offset Quantile"
      settingMatrix[1,2] = names(pamrTrain$offset)
      settingMatrix[2,1] = "Offset Value"
      settingMatrix[2,2] = pamrTrain$offset[[1]]
      settingMatrix[3,1] = "RNG Seed"
      settingMatrix[3,2] = input$randomSeed
      settingMatrix[4,1] = "Contrast Sign"
      settingMatrix[4,2] = pamrTrain$sign.contrast
      settingMatrix[5,1] = "Threshold"
      settingMatrix[5,2] = round(pamrConfusion$threshold, 5)
      
      prior = as.vector(pamrTrain$prior)
      settingsPriorMatrix = matrix(NA, nrow = 2, ncol = length(prior))
      rownames(settingsPriorMatrix) = c("Class", "Prob.")
      y = data$y
      settingsPriorMatrix[1,] = levels(y)
      settingsPriorMatrix[2,] = prior

      list(settingMatrix = settingMatrix, settingsPriorMatrix = settingsPriorMatrix)
    }  
  })
  
  
  
  
  ####################
  
  
  output$settings = renderTable({
    settings = getSettings()
    if(!is.null(settings)){
      settings$settingMatrix
    }
    
  }, include.rownames = FALSE)
  
  output$settingsPrior = renderTable({
    settings = getSettings()
    if(!is.null(settings)){
      settings$settingsPriorMatrix
    }
    
  }, include.colnames = FALSE)
  
  output$fdrPlot = renderPlot({
    fdr = getFDR()
    if(!is.null(fdr)){
      pamr.plotfdr.revised(fdr)
    }
  })
  
  output$fdr = renderTable({
    fdr = getFDR()
    if(!is.null(fdr)){
      fdr$result 
    }
  })
  
  output$pamrConfusion = renderTable({
    pamrConfusion = getPamrConfusion()
    if(!is.null(pamrConfusion)){
      tt = pamrConfusion$tt
      True__Predicted = dimnames(tt)[[2]][-ncol(tt)]
      cbind(True__Predicted, pamrConfusion$tt)
    }
  }, include.rownames=FALSE)
  
  
  output$pamrConfusionTrain = renderTable({
    pamrConfusion = getPamrConfusionTrain()
    if(!is.null(pamrConfusion)){
      tt = pamrConfusion$tt
      True__Predicted = dimnames(tt)[[2]][-ncol(tt)]
      cbind(True__Predicted, pamrConfusion$tt)
    }
  }, include.rownames = FALSE)
  
  output$dat = renderDataTable({
    data = getData()
    if(!is.null(data)){
      data$dat
    }
  })
  
  output$transform = renderDataTable({
    data = getData()
    if(!is.null(data)){
      if(input$cuberoot || input$center || input$scale){
        data$transformed
      }
    }
  })
  
  output$imputedX = renderDataTable({
    data = getData()
    if(!is.null(data$imputedX)){
      data$imputedX
    } 
  })
  
  output$listgenes = renderTable({
    listgenes = getPamrListGenes()
    if(!is.null(listgenes)){
      listgenes
    }
  })
  
  output$plotcv = renderPlot({
    pamrCv = getPamrCv()
    if(!is.null(pamrCv)){
      pamr.plotcv.revised(pamrCv, 1)
    }
  })
  
  output$plotcv2 = renderPlot({
    pamrCv = getPamrCv()
    if(!is.null(pamrCv)){
      pamr.plotcv.revised(pamrCv, 2)
    }
  })
  
  output$plotcvprob = renderPlot({
    pamrCv = getPamrCv()
    if(!is.null(pamrCv) && !is.na(input$threshold)){
      data = getData()
      pamr.plotcvprob(pamrCv, data, threshold = input$threshold)  
    }
  })
  
  output$plotcen = renderPlot({
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain) && !is.na(input$threshold)){
      data = getData()
      pamr.plotcen(pamrTrain, data, threshold = input$threshold)
    }
  })

  output$predict = renderTable({
    predict = getPamrPredict()
    if(!is.null(predict)){
      predict
    }
  })
  
  output$plotpredprob = renderPlot({
    testData = getTestData()
    pamrTrain = getPamrTrain()
    if(!is.null(testData) && !is.na(input$threshold)){
      pamr.plotpredprob(pamrTrain, testData, input$threshold)
    }
  })

  output$plotTrainError = renderPlot({
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain)){
      pamr.plotTrainError(pamrTrain)      
    }
  })
    
  output$plotTestError = renderPlot({
    pamrTrain = getPamrTrain()
    testData = getTestData()
    if(!is.null(testData) && !is.null(pamrTrain)){
      pamr.plotTestError(pamrTrain, testData)
    }
    
  })
  
  #############
  
  output$originalDataText = renderText({
    data = getData()
    if(!is.null(data)){
      "Original Data"
    }
  })
  
  output$transformDataText = renderText({
    data = getData()
    if(!is.null(data)){
      if(input$cuberoot || input$center || input$scale){
        "Transformed Data"
      }
    }    
  })
  
  output$imputedDataText = renderText({
    data = getData()
    if(!is.null(data$imputedX)){
      "Imputed Data"
    }    
  })
  
  
  output$listgenesText = renderText({
    listgenes = getPamrListGenes()    
    if(!is.null(listgenes) && !is.na(input$threshold)){
      "List of Significant Genes"
    }
  })
  
  output$centroidText = renderText({
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain) && !is.na(input$threshold)){  
      "Shrunken class centroid"      
    }
  })
  
  output$confusionTrainText = renderText({
    pamrConfusion = getPamrConfusionTrain()
    if(!is.null(pamrConfusion) && !is.na(input$threshold)){
      "Training Confusion Matrix"
    }
  })
  
  output$cvConfusionMatrix <- renderText({
    pamrConfusion = getPamrConfusion()
    if(!is.null(pamrConfusion)){
      "CV Confusion Matrix"
    }
  })
  
  output$settingsText = renderText({
    settings = getSettings()
    if(!is.null(settings)){
      "Settings"
    }
  })
  
  output$fdrText = renderText({
    fdr = getFDR()
    if(!is.null(fdr)){
      "FDR results"
    }
  })
  
  output$fdrPlotText = renderText({
    fdr = getFDR()
    if(!is.null(fdr)){
      "FDR Plot"
    }
  })
  
  output$overallText = renderText({
    pamrCv = getPamrCv()
    if(!is.null(pamrCv)){
      "Overall CV Plot" 
    }
  })
  
  output$individualText = renderText({
    pamrCv = getPamrCv()
    if(!is.null(pamrCv)){
      "Individual CV Plots" 
    }
  })
  
  output$plotcvText = renderText({
    data = getData()
    pamrCv = getPamrCv()  
    if(!is.null(pamrCv) && !is.na(input$threshold)){
      "CV Probabilities Plot"
    }
  })
  
  
  output$settingsPriorText = renderText({
    settings = getSettings()
    if(!is.null(settings)){
      "Prior Distribution"
    }    
  })
  
  output$predictTableText = renderText({
    pamrPredict =getPamrPredict()
    if(!is.null(pamrPredict)){
      "Actual, Predicted Classes and Predicted Posterior Probabilities"
    }
  })
  
  output$predictPlotText = renderText({
    pamrPredict =getPamrPredict()
    if(!is.null(pamrPredict)){
      "Test Probabilities"
    }
  })
  
  output$testErrorPlotText = renderText({
    pamrPredict =getPamrPredict()
    if(!is.null(pamrPredict)){
      "Test Error"
    }
  })
  
  output$trainErrorPlotText = renderText({
    pamrTrain =getPamrTrain()
    if(!is.null(pamrTrain)){
      "Train Error"
    }
  })
  
  output$threshold <- renderUI({
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain)){
      maxvalue = round(pamrTrain$threshold[length(pamrTrain$threshold)], 5)
      numericInput("threshold", label = paste("Threshold", "(max=", maxvalue, ")", sep = ""), min = 0, max = maxvalue, value = NULL, step = 0.000001)  
    }
  })
  
  output$customvalue = renderUI({
    data = getData()
    if(!is.null(data)){
      y = data$y
      lengthy = length(levels(y))
      customnamelabel = paste(levels(y), collapse = ",")
      customnamevalue = rep(1/lengthy, lengthy)
      customnamevalue = paste(customnamevalue, collapse = ",")
      textInput("custom", label = customnamelabel, value = customnamevalue)
    }
  })
  

  ###################survival and regression problem
  
  getDataSurv = reactive({
  
  objFile = input$iFile
  survivalTimeLabel = input$survivalTimeLabel
  censoringStatusLabel = input$censoringStatusLabel
  sampleLabel = input$sampleLabel
  expressionStart = input$expressionStart
  batchLabel = input$batchLabel
  outcomeValueLabel = input$outcomeValueLabel
  
  if(!is.null(objFile) && ((!is.na(survivalTimeLabel) && !is.na(censoringStatusLabel)) || !is.na(outcomeValueLabel)) && !is.na(expressionStart)){
    
    dat = read.xlsx(objFile$datapath, 1, colNames = FALSE)   
    emptySpace = 1:(expressionStart-1)
    
    geneid = dat[-emptySpace,1]
    genenames = dat[-emptySpace,2]   
    
    if(!is.na(sampleLabel)){
      sample = as.matrix(dat[sampleLabel, c(-1,-2)])
    }
        
    if(!is.na(survivalTimeLabel)){
      y = as.numeric(dat[survivalTimeLabel, c(-1,-2)])  
    } else if(!is.na(outcomeValueLabel)){
      y = as.numeric(dat[outcomeValueLabel, c(-1,-2)]) 
    }
    
    
    x = as.matrix(dat[-emptySpace, c(-1,-2)])
    class(x) = "numeric"
    
    censoring.status = NULL
    if(!is.na(censoringStatusLabel)){
      censoring.status = as.numeric(dat[censoringStatusLabel,c(-1,-2)])  
    }
    
    imputedX = NULL
      if(sum(is.na(x)) > 0){  
        imputedummy = impute.knn(x, k = input$numberOfNeighbors)
        x = imputedummy$data
        imputedX = cbind(geneid, genenames, x)
    }
    
    transformed = NULL
    if(!is.null(input$center) || !is.null(input$scale) || !is.null(input$cuberoot)){
      x = scale(x, center = input$center, scale = input$scale)  
      if(input$cuberoot){
        x = sign(x) * abs(x)^(1/3)
      }
      transformed = cbind(geneid, genenames, x)
    }

    if(!is.na(batchLabel)){
      batchLabels = as.matrix(dat[batchLabel, c(-1,-2)])
      mydata <- list(x=x,y=y,batchlabels=factor(batchLabels))
      mydata2 <- pamr.batchadjust(mydata)      
      x = mydata2$x
    }
    
    decorrelateX = NULL
    if(!is.null(input$competingPredictorFile)){
      competing.predictors = getCompetingPredictors()
      dec = superpc.decorrelate(x,competing.predictors)
      decorrelateX = t(dec$res)
      x = decorrelateX
      decorrelateX = cbind(geneid, genenames, decorrelateX)
    }
    
    data = list(x = x, y = y, featurenames = genenames, geneid = geneid, dat = dat, sample = sample, censoring.status = censoring.status, decorrelateX = decorrelateX, transformed = transformed, imputedX = imputedX)
    data
  }
  })
  
  
  getTestDataSurv = reactive({
    
    objFile = input$testFileSurv
    survivalTimeLabel = input$survivalTimeLabel
    censoringStatusLabel = input$censoringStatusLabel
    sampleLabel = input$sampleLabel
    outcomeValueLabel = input$outcomeValueLabel
    expressionStart = input$expressionStart
    
    if(!is.null(objFile) && ((!is.na(survivalTimeLabel) && !is.na(censoringStatusLabel)) || !is.na(outcomeValueLabel))  && !is.na(expressionStart)){
      
      dat = read.xlsx(objFile$datapath, 1, colNames = FALSE)   
      emptySpace = 1:(expressionStart-1)
      
      geneid = dat[-emptySpace,1]
      genenames = dat[-emptySpace,2]   
      
      if(!is.na(sampleLabel)){
        sample = as.matrix(dat[sampleLabel, c(-1,-2)])  
      }
      
      if(!is.na(survivalTimeLabel)){
        y = as.numeric(dat[survivalTimeLabel, c(-1,-2)])  
      }
      
      else if(!is.na(outcomeValueLabel)){
        y = as.numeric(dat[outcomeValueLabel, c(-1,-2)])
      }
      
      x = as.matrix(dat[-emptySpace, c(-1,-2)])
      class(x) = "numeric"
      
      censoring.status = NULL
      if(!is.na(censoringStatusLabel)){
        censoring.status = as.numeric(dat[censoringStatusLabel,c(-1,-2)])  
      }
            
      decorrelateX = NULL
      if(!is.null(input$competingPredictorFile)){
        competing.predictors = getCompetingPredictors()
        dec = superpc.decorrelate(x,competing.predictors)
        decorrelateX = t(dec$res)
        x = decorrelateX
        decorrelateX = cbind(geneid, genenames, decorrelateX)
      }
      
      data = list(x = x, y = y, featurenames = genenames, geneid = geneid, dat = dat, sample = sample, censoring.status = censoring.status, decorrelateX = decorrelateX)
      data
    }
  })
  
  
  getCompetingPredictors = reactive({
    
    comFile = input$competingPredictorFile
    if(!is.null(comFile)){
      dat = read.xlsx(comFile$datapath, 1, colNames = TRUE)   
      n = length(dat[,1])
      competingPredictors = vector("list", n)
      names(competingPredictors) = dat[,1]
      for(i in 1:n){
        competingPredictors[[i]] = as.matrix(dat[i,c(-1,-2)])
        if(dat[i,2] == "discrete"){
          competingPredictors[[i]] = factor(competingPredictors[[i]])
        } else{
          competingPredictors[[i]] = as.numeric(competingPredictors[[i]])
        }
      }
      competingPredictors
    }
    
  })
  
  getCompetingPredictorFit = reactive({
    
    comFile = input$competingPredictorFitFile
    if(!is.null(comFile)){
      dat = read.xlsx(comFile$datapath, 1, colNames = TRUE)   
      n = length(dat[,1])
      competingPredictors = vector("list", n)
      names(competingPredictors) = dat[,1]
      for(i in 1:n){
        competingPredictors[[i]] = as.matrix(dat[i,c(-1,-2)])
        if(dat[i,2] == "discrete"){
          competingPredictors[[i]] = factor(competingPredictors[[i]])
        } else{
          competingPredictors[[i]] = as.numeric(competingPredictors[[i]])
        }
      }
      sample = as.matrix(colnames(dat)[c(-1,-2)])
      list(competing.predictors = competingPredictors, sample.labels = sample)
    }
  })  
  
  getSurvTrain = reactive({
    data = getDataSurv()
    if(!is.null(data)){
      if(input$analysisType == "Survival"){
        train.obj = superpc.train(data, type = "survival")    
      } else if (input$analysisType == "Regression"){
        train.obj = superpc.train(data, type = "regression")
      }
      print(input$analysisType)
      train.obj
    }
  })
  
  getSurvCv = reactive({
    data = getDataSurv()
    train.obj = getSurvTrain()
    
    set.seed(input$randomSeed)
    if(!is.null(data)){
      cv.obj = superpc.cv(train.obj, data)
      cv.obj      
    }
  })
  
  getlrtestObj = reactive({
    data = getDataSurv()
    train.obj = getSurvTrain()
    if(!is.null(train.obj)){
      lrtest.obj = superpc.lrtest.curv(train.obj, data, data)
      lrtest.obj
    }
  })
  
  getlrtestObjTest = reactive({
    data = getDataSurv()
    train.obj = getSurvTrain()
    testData = getTestDataSurv()
    if(!is.null(train.obj) && !is.null(testData)){
      lrtest.obj = superpc.lrtest.curv(train.obj, data, testData)
      lrtest.obj
    }
  })
  
  getPredictRed = reactive({
    data = getDataSurv()
    train.obj = getSurvTrain()
    if(!is.null(train.obj) && !is.na(input$threshold2) && !is.na(input$princomp)){
      fit.red = superpc.predict.red(train.obj, data, data, threshold = input$threshold2, n.components = input$princomp)
      fit.red
    }
  })
  
  getPredictRedFixed = reactive({
    data = getDataSurv()
    train.obj = getSurvTrain()
    if(!is.null(train.obj) && !is.na(input$threshold2)){
      fit.red = superpc.predict.red(train.obj, data, data, threshold = input$threshold2)
      fit.red
    }
  })
  
  getPredictRedCv = reactive({
    cv.obj = getSurvCv()
    data = getDataSurv()
    fit.red = getPredictRed()
    if(!is.null(fit.red) && !is.na(input$threshold2)){
      fit.redcv = superpc.predict.red.cv(fit.red, cv.obj, data, threshold = input$threshold2)
      fit.redcv
    }
  })
  
  getPredictRedTest = reactive({
    train.obj = getSurvTrain()
    data = getDataSurv()
    testData = getTestDataSurv()
    if(!is.null(train.obj) && !is.na(input$threshold2) && !is.null(testData) && !is.na(input$shrinkage)){
      fit.red = superpc.predict.red(train.obj, data, testData, threshold = input$threshold2, shrinkages = input$shrinkage, n.components = 1)
      fit.red
    }
  })

  
  getListFeatures = reactive({
    data = getDataSurv()
    train.obj = getSurvTrain()
    fit.red = getPredictRed()
    
    if(!is.null(fit.red) && !is.na(input$princomp)){
      listfeatures = superpc.listfeatures(data, train.obj, fit.red, component.number = input$princomp)
      listfeatures
    }
  })

  getPredict = reactive({
    
    train.obj = getSurvTrain()
    data = getDataSurv()
    testData = getTestDataSurv()
    
    if(!is.null(train.obj) && !is.null(testData) && !is.na(input$threshold2) && !is.na(input$princomp)){
      predict = superpc.predict(train.obj, data, testData, threshold = input$threshold2, n.components = input$princomp)  
      predict
    }
  })
  
  getPredictDiscrete = reactive({
    train.obj = getSurvTrain()
    data = getDataSurv()
    testData = getTestDataSurv()
    if(!is.null(train.obj) && !is.null(testData) && !is.na(input$threshold2)){
      predict = superpc.predict(train.obj, data, testData, threshold = input$threshold2, n.components = 1, prediction.type = "discrete")  
      predict
    }
  })
  
  getOutcome = reactive({
    testData = getTestDataSurv()
    train.obj = getSurvTrain()
    fit.groups = getPredict()
    fit = getPredictRedTest()
    if(!is.null(testData) && !is.null(fit.groups)){
      if(is.na(input$shrinkage)){
        outcome = superpc.fit.to.outcome(train.obj, testData, fit.groups$v.pred)  
      }else{
        outcome = superpc.fit.to.outcome(train.obj, testData, as.numeric(fit$v.test))
      }
      outcome
    }
  })
  

  
  ##########

  
  output$rainbowPlot = renderPlot({
    competingPredictorFit = getCompetingPredictorFit()
    data = getDataSurv()
    fit = getPredict()
    testData = getTestDataSurv()
    fit.red = getPredictRedTest()
    
    if(!is.null(data) && !is.null(testData) && !is.null(competingPredictorFit) && !is.na(input$threshold2) && !is.na(input$princomp)){
      sample.labels = competingPredictorFit$sample.labels
      competing.predictors = competingPredictorFit$competing.predictors
      if(is.na(input$shrinkage)){
        superpc.rainbowplot(data, fit$v.pred, sample.labels, competing.predictors)
      }else{
        superpc.rainbowplot(data, as.numeric(fit.red$v.test), sample.labels, competing.predictors)
      }      
    }
    
  })
  
  
  output$responsePredictionPlot = renderPlot({
    predict = getPredictDiscrete()
    testData = getTestDataSurv()
    if(!is.null(predict)){
      plot(survfit(Surv(testData$y,testData$censoring.status)~predict$v.pred), col=2:3, xlab="time", ylab="Prob survival")  
    }
  })
  
  output$predictionscore = renderTable({
    predict = getPredict()
    testData = getTestDataSurv()
    fit = getPredictRedTest()
    if(!is.null(predict)){
      if(is.na(input$shrinkage)){
        predictMatrix = matrix(predict$v.pred.1df, nrow = 1)  
      } else{
        predictMatrix = matrix(as.numeric(fit$v.test.1df), nrow = 1)
      }
      
      rownames(predictMatrix) = "Prediction Score"
      if(!is.null(testData$sample)){
        colnames(predictMatrix) = testData$sample
      }
      predictMatrix
    }
  }, digits = 5)
  
  output$teststatTable = renderTable({
    outcome = getOutcome()
    if(!is.null(outcome)){
      test = outcome$teststat.table
      rownames(test) = "Overall model fit"
      test
    }

  }, digits = 5)
  
  output$coeftable = renderTable({
    outcome= getOutcome()
    if(!is.null(outcome)){
      result = outcome$coeftable
    }
  }, digits = 5)
  
  output$plotredLrtest = renderPlot({
    fit.redcv = getPredictRedCv()
    if(!is.null(fit.redcv)){
      superpc.plotred.lrtest(fit.redcv)  
    }
  })
    
    
  output$imputedXSurv = renderDataTable({
    data = getDataSurv()
    if(!is.null(data$imputedX)){
      data$imputedX
    }
  })
  
  output$transformSurv = renderDataTable({
    data = getDataSurv()
    if(!is.null(data)){
      if(input$cuberoot || input$center || input$scale){
        data$transformed
      }
    }
  })
  
  output$plotcvsurv = renderPlot({
    cv.obj = getSurvCv()
    if(!is.null(cv.obj)){
      superpc.plotcv(cv.obj)      
    }
  })
  
  output$listfeatures = renderTable({
    listfeatures = getListFeatures()
    if(!is.null(listfeatures)){
      listfeatures
    }
  })
  
  output$threshold2 = renderUI({
    cv.obj = getSurvCv()
    if(!is.null(cv.obj)){
      maxvalue = round(cv.obj$threshold[length(cv.obj$threshold)], 5)
      numericInput("threshold2", label = paste("Threshold", "(max=", maxvalue, ")", sep = ""), min = 0, max = maxvalue, value = NULL, step = 0.000001)  
    }
  })
  
  output$shrinkage = renderUI({
    fit.red = getPredictRedFixed()
    if(!is.null(fit.red)){
      maxvalue = round(fit.red$shrinkages[length(fit.red$shrinkages)], 5)
      numericInput("shrinkage", label = paste("Shrinkage", "(max=", maxvalue, ")", sep = ""), min = 0, max = maxvalue, value = NULL, step = 0.000001)
    }
  })
  
  output$survdata = renderDataTable({
    data = getDataSurv()
    if(!is.null(data)){
      data$dat
    }
  })

  output$plotLrtest = renderPlot({
    lrtest.obj = getlrtestObj()
    if(!is.null(lrtest.obj)){
      superpc.plot.lrtest(lrtest.obj)  
      axis(3, at = lrtest.obj$threshold, labels = lrtest.obj$num.features) 
    }
  })
  
  output$plotLrtestTest = renderPlot({
    lrtest.obj = getlrtestObjTest()
    if(!is.null(lrtest.obj)){
      superpc.plot.lrtest(lrtest.obj)
      axis(3, at = lrtest.obj$threshold, labels = lrtest.obj$num.features)
    }
  })
  
  output$survPredictionPlot = renderPlot({
    train.obj = getSurvTrain()
    data = getDataSurv()
    if(!is.null(train.obj) && !is.na(input$threshold2)){
      superpc.predictionplot(train.obj, data, data, threshold = input$threshold2)      
    }
  })
  
  output$decorrelateX = renderDataTable({
    data = getDataSurv()
    if(!is.null(data)){
      data$decorrelateX      
    }
  })
  
  
  ########
  
  
  output$rainbowPlotText = renderText({
    competingPredictorFit = getCompetingPredictorFit()
    data = getDataSurv()
    testData = getTestDataSurv()
    if(!is.null(data) && !is.null(testData) && !is.null(competingPredictorFit) && !is.na(input$threshold2) && !is.na(input$princomp)){
      "Rainbow Plot"  
    }
  })
  
  output$responsePredictionPlotText = renderText({
    predict = getPredictDiscrete()
    if(!is.null(predict)){
      "Response Prediction Plot"
    }
  })
  
  output$predictionInfoText = renderText({
    predict = getPredict()
    if(!is.null(predict)){
      "Prediction info"
    }
    
  })
  
  output$originalXText = renderText({
    data = getDataSurv()
    if(!is.null(data)){
      "Original Data"      
    }
  })
  
  output$decorrelateXText = renderText({
    data = getDataSurv()
    if(!is.null(data)){
      if(!is.null(data$decorrelateX)){
        "Decorrelated X"
      }      
    }
  })

  
  output$survTrainErrorText = renderText({
    data = getDataSurv()
    if(!is.null(data)){
      "Training Error"
    }
  })
  
  output$plotredLrtestText = renderText({
    fit.redcv = getPredictRedCv()
    if(!is.null(fit.redcv)){
      "Reduced Model CV"  
    }
    
  })
  
  output$listSurvGenesText = renderText({
    listgenes = getListFeatures()
    if(!is.null(listgenes)){
      "List of Significant Genes"
    }
  })
  
  output$responsePredictionText = renderText({
    train.obj = getSurvTrain()
    if(!is.null(train.obj) && !is.na(input$threshold2)){  
      "Response Prediction Plot"
    }
  })
  
  output$transformSurvText = renderText({
    data = getDataSurv()
    if(!is.null(data)){
      if(input$cuberoot || input$center || input$scale){
        "Transformed X"
      }
    }
  })
  
  output$imputedXSurvText = renderText({
    data = getDataSurv()
    if(!is.null(data$imputedX)){
      "Imputed X"
    }
  })
  
  output$plotCvSurvText = renderText({
    cv.obj = getSurvCv()
    if(!is.null(cv.obj)){
      "CV Curves"      
    }
  })
  
  output$lrtestObjTestText = renderText({
    lrtestObjTest = getlrtestObjTest()
    if(!is.null(lrtestObjTest)){
      "Test Error"
    }
  })
  

  
})