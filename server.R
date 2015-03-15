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
      
      data = list(x = x, y = y, genenames = genenames, geneid = geneid, dat = dat, sample = sample, imputedX = imputedX)
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
      
      list(x = x, y = y, sample = sample, dat = dat)
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

  
  getPamrTrainExact = reactive({
    dat = getData()
    if(!is.null(dat)){
      data = list(x = dat$x, y = dat$y, genenames = dat$genenames, geneid = dat$geneid)    
      pamr.train(data, threshold = input$threshold, offset.percent = input$s0percentile, sign.contrast = input$sign.contrast, prior = getPrior())
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

  getPamrCvExact = reactive({
    dat = getData()
    if(!is.null(dat)){
      set.seed(input$randomSeed)
      fit = getPamrTrainExact()
      pamr.cv(fit, dat)
    }
  })
  
  
  getPamrConfusionTrain = reactive({
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain) && !is.na(input$threshold)){
      pamrTrainExact = getPamrTrainExact()
      maxvalue = pamrTrain$threshold[length(pamrTrain$threshold)]
      
      if(input$threshold <= maxvalue){
        pamrConfusion = pamr.confusion.revised(pamrTrainExact, threshold = input$threshold)
        pamrConfusion
      }
    }
  })
  
  getPamrConfusion = reactive({
    pamrCv = getPamrCv()
    if(!is.null(pamrCv) && !is.na(input$threshold)){
      
      pamrCvExact = getPamrCvExact()
      pamrTrain = getPamrTrain()
      maxvalue = pamrTrain$threshold[length(pamrTrain$threshold)]
            
      if(input$threshold <= maxvalue){
        pamrConfusion = pamr.confusion.revised(pamrCvExact, threshold = input$threshold)
        pamrConfusion
      }
    }
  })
  
  getPamrListGenes = reactive({
    
    pamrTrain = getPamrTrain()
    if(!is.null(pamrTrain) && !is.na(input$threshold)){
      pamrTrainExact = getPamrTrainExact()
      data = getData()
      pamr.listgenes(pamrTrainExact, data, threshold = input$threshold, genenames = TRUE)
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
    data = getData()
    
    if(!is.null(pamrTrain) ){
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
      settingMatrix[5,2] = if(is.na(input$threshold)){NA}else{input$threshold}
      
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
  
  output$testdat = renderDataTable({
    data = getTestData()
    if(!is.null(data)){
      data$dat
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
      pamrTrainExact = getPamrTrainExact()
      pamr.plotcen(pamrTrainExact, data, threshold = input$threshold)
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
      "Data"
    }
  })
  
  output$testDataText = renderText({
    data = getTestData()
    if(!is.null(data)){
      "Test Data"
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
    pamrTrain = getPamrTrain()
    testData = getTestData()
    if(!is.null(testData) && !is.null(pamrTrain)){
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
  
  
  #######Save Classification results
  savethis = observe({
    
  if(input$saveButton != 0){
    
  isolate({  
  dir = input$dir
  file = input$fname
      
  pamrTrain = getPamrTrain()
  pamrConfusion = getPamrConfusionTrain()
  listgenes = getPamrListGenes()
  data = getData()
  testData = getTestData()
  fdr = getFDR()
  pamrCv = getPamrCv()
  pamrConfusionCV = getPamrConfusion()
  settings = getSettings()  
  predict = getPamrPredict()
  
  titleStyle = createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD")
  wb = createWorkbook()  
  
  if(!is.null(data)){
    addWorksheet(wb, sheetName = "Data")
    writeData(wb, sheet = "Data", x = data$dat)    
  }
  
  if(!is.null(data)){
    if(!is.null(data$imputedX)){
      addWorksheet(wb, sheetName = "Imputed Data")
      writeData(wb, sheet = "Imputed Data", x = data$imputedX)
    }
  }
  
  if(!is.null(testData)){
    addWorksheet(wb, sheetName = "Test Data")
    writeData(wb, sheet = "Test Data", x = testData$dat)
  }
  
  if(!is.null(pamrTrain)){
    png(file = "TrainError.png")
    pamr.plotTrainError(pamrTrain)
    dev.off()
  }
  
  if(!is.null(pamrTrain)){      
    addWorksheet(wb, sheetName = "Training")
    addStyle(wb, sheet = "Training", titleStyle, rows = 1, cols = 1)
    writeData(wb, sheet = "Training", x = "Train Error")  
    insertImage(wb, sheet = "Training", file = "TrainError.png", width = 4, height = 4, startRow = 2)
    setColWidths(wb, "Training", cols= 1, widths = 30)    
  }
  
  TrainingRow = 23
  if(!is.null(pamrConfusion)){
    tt = pamrConfusion$tt
    True__Predicted = dimnames(tt)[[2]][-ncol(tt)]
    pamrConfusionTrain = cbind(True__Predicted, pamrConfusion$tt)
    addStyle(wb, sheet = "Training", titleStyle, rows = TrainingRow, cols = 1)
    writeData(wb, sheet = "Training", x = "Training Confusion Matrix", startRow = TrainingRow)
    writeData(wb, pamrConfusionTrain, sheet = "Training", startRow = TrainingRow + 1)
    TrainingRow = TrainingRow + nrow(pamrConfusionTrain) + 3    
  }
  
  if(!is.null(listgenes)){
    addStyle(wb, sheet = "Training", titleStyle, rows = TrainingRow, cols = 1)
    writeData(wb, sheet = "Training", x = "List of Significant Genes", startRow = TrainingRow)
    writeData(wb, sheet = "Training", x = listgenes, startRow = TrainingRow + 1)
    TrainingRow = TrainingRow + nrow(listgenes) + 3 
  }
      
  if(!is.null(pamrTrain) && !is.na(input$threshold)){
    png(file = "centroid.png")
    pamr.plotcen(pamrTrain, data, threshold = input$threshold)
    dev.off()
  }
  
  if(!is.null(pamrTrain) && !is.na(input$threshold)){
    addStyle(wb, sheet = "Training", titleStyle, rows = TrainingRow, cols = 1)
    writeData(wb, sheet = "Training", x = "Shrunken class centroid", startRow = TrainingRow)
    insertImage(wb, sheet = "Training", file = "centroid.png", width = 4, height = 4, startRow = TrainingRow + 1)
    TrainingRow = TrainingRow + 22
  }
  
  if(!is.null(fdr)){
    addStyle(wb, sheet = "Training", titleStyle, rows = TrainingRow, cols = 1)
    writeData(wb, sheet = "Training", x = "FDR results", startRow = TrainingRow)
    writeData(wb, sheet = "Training", x = fdr$result, startRow = TrainingRow + 1)
    TrainingRow = TrainingRow + nrow(fdr$results) + 3 
  }
  
  if(!is.null(fdr)){
    png(file = "fdr.png")
    pamr.plotfdr.revised(fdr)
    dev.off()
  }
  
  if(!is.null(fdr)){
    addStyle(wb, sheet = "Training", titleStyle, rows = TrainingRow, cols = 1)
    writeData(wb, sheet = "Training", x = "FDR Plot", startRow = TrainingRow)
    insertImage(wb, sheet = "Training", file = "fdr.png", width = 4, height = 4, startRow = TrainingRow + 1)
    TrainingRow = TrainingRow + 22
  }
  
  if(!is.null(pamrCv)){
    png(file = "plotcv1.png")
    pamr.plotcv.revised(pamrCv, 1)
    dev.off()
  }
  
  if(!is.null(pamrCv)){
    png(file = "plotcv2.png")
    pamr.plotcv.revised(pamrCv, 2)
    dev.off()
  }
  
  cvRow = 0
  if(!is.null(pamrCv)){
    addWorksheet(wb, sheetName = "Cross Validation")
    setColWidths(wb, "Cross Validation", cols= 1, widths = 30)    
    addStyle(wb, sheet = "Cross Validation", titleStyle, rows = 1, cols = 1)
    writeData(wb, sheet = "Cross Validation", x = "Overall CV Plot")
    insertImage(wb, sheet = "Cross Validation", file = "plotcv1.png", width = 4, height = 4, startRow = 2)
    addStyle(wb, sheet = "Cross Validation", titleStyle, rows = 24, cols = 1)
    writeData(wb, sheet = "Cross Validation", x = "Individual CV Plots", startRow = 24)
    insertImage(wb, sheet = "Cross Validation", file = "plotcv2.png", width = 4, height = 4, startRow = 25)
    cvRow = 47
  }
  
  if(!is.null(pamrCv) && !is.na(input$threshold)){
    png(file = "cvprob.png")
    pamr.plotcvprob(pamrCv, data, threshold = input$threshold)
    dev.off()
  }
  
  if(!is.null(pamrCv) && !is.na(input$threshold)){
    addStyle(wb, sheet = "Cross Validation", titleStyle, rows = cvRow, cols = 1)
    writeData(wb, sheet = "Cross Validation", x = "CV Probabilities Plot", startRow = cvRow)
    insertImage(wb, sheet = "Cross Validation", file = "cvprob.png", width = 4, height = 4, startRow = cvRow+1)
    cvRow = cvRow + 22    
    
    addStyle(wb, sheet = "Cross Validation", titleStyle, rows = cvRow, cols = 1)
    writeData(wb, sheet = "Cross Validation", x = "CV Confusion Matrix", startRow = cvRow)
    tt = pamrConfusionCV$tt
    True__Predicted = dimnames(tt)[[2]][-ncol(tt)]
    CVConf = cbind(True__Predicted, pamrConfusion$tt)
    writeData(wb, sheet = "Cross Validation", x = CVConf, startRow = cvRow+1)    
  }
  
  if(!is.null(testData) && !is.null(pamrTrain)){
    png(file = "testerror.png")
    pamr.plotTestError(pamrTrain, testData)
    dev.off()
  }
  
  
  if(!is.null(testData) && !is.null(pamrTrain)){
    addWorksheet(wb, sheetName = "Test Set Prediction")
    setColWidths(wb, "Test Set Prediction", cols= 1, widths = 30)    
    addStyle(wb, sheet = "Test Set Prediction", titleStyle, rows = 1, cols = 1)
    writeData(wb, sheet = "Test Set Prediction", x = "Test Error")
    insertImage(wb, sheet = "Test Set Prediction", file = "testerror.png", width = 4, height = 4, startRow = 2)
  }
  
  if(!is.null(testData) && !is.null(pamrTrain) && !is.na(input$threshold)){
    png(file = "testprob.png")
    pamr.plotpredprob(pamrTrain, testData, input$threshold)
    dev.off()
  }
  
  if(!is.null(testData) && !is.null(pamrTrain) && !is.na(input$threshold)){
    testRow = 23
    addStyle(wb, sheet = "Test Set Prediction", titleStyle, rows = testRow, cols = 1)
    writeData(wb, sheet = "Test Set Prediction", x = "Test Probabilities", startRow = testRow)
    insertImage(wb, sheet = "Test Set Prediction", file = "testprob.png", width = 4, height = 4, startRow = testRow + 1)
    testRow = 47
    addStyle(wb, sheet = "Test Set Prediction", titleStyle, rows = testRow, cols = 1)
    writeData(wb, sheet = "Test Set Prediction", x = "Actual, Predicted Classes and Predicted Posterior Probabilities", startRow = testRow)
    writeData(wb, sheet = "Test Set Prediction", x = predict, startRow = testRow + 1)
  }

  if(!is.null(settings)){
    addWorksheet(wb, sheetName = "Settings")
    setColWidths(wb, "Settings", cols= 1, widths = 30)    
    addStyle(wb, sheet = "Settings", titleStyle, rows = 1, cols = 1)
    writeData(wb, sheet = "Settings", x = "Settings", startRow = 1)
    writeData(wb, sheet = "Settings", x = settings$settingMatrix, startRow = 2)
    addStyle(wb, sheet = "Settings", titleStyle, rows = 9, cols = 1)
    writeData(wb, sheet = "Settings", x = "Prior Distribution", startRow = 9)
    priordummy = settings$settingsPriorMatrix[2,]
    priorMatrix = matrix(priordummy, ncol = length(priordummy))
    colnames(priorMatrix) = settings$settingsPriorMatrix[1,]
    writeData(wb, sheet = "Settings", x = priorMatrix, startRow = 10)
  }
  
  
  if(!is.null(pamrTrain)){
    fname = paste(file, "xlsx", sep = ".")
    saveWorkbook(wb, file.path(dir, fname), overwrite = TRUE)
  
    if(Sys.info()[['sysname']] == "Windows"){
      shell(file.path(dir, fname))  
    } else if(Sys.info()[['sysname']] == "Darwin"){
      system(paste("open", fname))
    }
  }
  
  })
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
    
    data = list(x = x, y = y, featurenames = genenames, geneid = geneid, dat = dat, sample = sample, censoring.status = censoring.status, decorrelateX = decorrelateX, imputedX = imputedX)
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
    if(!is.null(train.obj) && !is.na(input$threshold2) && !is.na(input$princomp) && !is.null(data$censoring.status)){
      fit.red = superpc.predict.red(train.obj, data, data, threshold = input$threshold2, n.components = input$princomp)
      fit.red
    }
  })
  
  getPredictRedFixed = reactive({
    data = getDataSurv()
    train.obj = getSurvTrain()
    if(!is.null(train.obj) && !is.na(input$threshold2) && !is.null(data$censoring.status)){
      fit.red = superpc.predict.red(train.obj, data, data, threshold = input$threshold2)
      fit.red
    }
  })
  
  getPredictRedCv = reactive({
    cv.obj = getSurvCv()
    data = getDataSurv()
    fit.red = getPredictRed()
    if(!is.null(fit.red) && !is.na(input$threshold2) && !is.null(data$censoring.status)){
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
    
    if(!is.null(data) && !is.null(testData) && !is.null(competingPredictorFit) && !is.na(input$threshold2) && !is.na(input$princomp) && !is.null(data$censoring.status)){
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
    if(!is.null(predict) && !is.null(testData$censoring.status)){
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
    data = getDataSurv()
    
    if(!is.null(fit.redcv) && !is.null(data$censoring.status)){
      superpc.plotred.lrtest(fit.redcv)  
    }
  })
    
    
  output$imputedXSurv = renderDataTable({
    data = getDataSurv()
    if(!is.null(data$imputedX)){
      data$imputedX
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

    lrtest.obj = getlrtestObj()
    if(!is.null(lrtest.obj)){
      maxvalue = round(lrtest.obj$threshold[length(lrtest.obj$threshold)], 5)
      numericInput("threshold2", label = paste("Threshold", "(max=", maxvalue, ")", sep = ""), min = 0, max = maxvalue, value = NULL, step = 0.000001)  
    }
  })
    
  output$survdata = renderDataTable({
    data = getDataSurv()
    if(!is.null(data)){
      data$dat
    }
  })
  
  output$survTestData = renderDataTable({
    testData = getTestDataSurv()
    if(!is.null(testData)){
      testData$dat
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
  
  
  ########
  
  
  output$rainbowPlotText = renderText({
    competingPredictorFit = getCompetingPredictorFit()
    data = getDataSurv()
    testData = getTestDataSurv()
    if(!is.null(data) && !is.null(testData) && !is.null(competingPredictorFit) && !is.na(input$threshold2) && !is.na(input$princomp) && !is.null(data$censoring.status)){
      "Rainbow Plot"  
    }
  })
  
  output$responsePredictionPlotText = renderText({
    predict = getPredictDiscrete()
    data = getDataSurv()
    if(!is.null(predict) && !is.null(data$censoring.status) ){
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
      "Data"      
    }
  })
  

  
  output$survTrainErrorText = renderText({
    data = getDataSurv()
    if(!is.null(data)){
      "Training Error"
    }
  })
  
  output$testSurvText = renderText({
    testData = getTestDataSurv()
    if(!is.null(testData)){
      "Test Data"
    }
  })
  
  output$plotredLrtestText = renderText({
    fit.redcv = getPredictRedCv()
    data = getDataSurv()
    if(!is.null(fit.redcv) && !is.null(data$censoring.status)){
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
  
  ####save result for survival
  
  savethis2 = observe({
    
  if(input$saveButton2 != 0){
      
  isolate({  
    dir = input$dir
    file = input$fname
        
    data = getDataSurv()
    testData = getTestDataSurv()
    train.obj = getSurvTrain()
    lrtest.obj = getlrtestObj()
    listfeatures = getListFeatures()
    cv.obj = getSurvCv()
    fit.redcv = getPredictRedCv()
    lrtest.obj.test = getlrtestObjTest()
    predict = getPredict()
    outcome = getOutcome()
    predictDiscrete = getPredictDiscrete()
    competingPredictorFit = getCompetingPredictorFit()
    fit.red = getPredictRedTest()
    fit = getPredictRedTest()   
    
    titleStyle = createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD")
    wb = createWorkbook()  
    
    if(!is.null(data)){
      addWorksheet(wb, sheetName = "Data")
      writeData(wb, sheet = "Data", x = data$dat, startRow = 1)      
    }
    
    if(!is.null(testData)){
      addWorksheet(wb, sheetName = "Test Data")
      writeData(wb, sheet = "Test Data", x = testData$dat, startRow = 1)      
    }
    
    if(!is.null(lrtest.obj)){
      png(file = "lrtesttrain.png")
      superpc.plot.lrtest(lrtest.obj)  
      axis(3, at = lrtest.obj$threshold, labels = lrtest.obj$num.features) 
      dev.off()
    }
    
    trainRow = 0
    if(!is.null(lrtest.obj)){
      addWorksheet(wb, sheetName = "Training")
      setColWidths(wb, "Training", cols= 1, widths = 30)    
      addStyle(wb, sheet = "Training", titleStyle, rows = 1, cols = 1)
      writeData(wb, sheet = "Training", x = "Training Error", startRow = 1)
      insertImage(wb, sheet = "Training", file = "lrtesttrain.png", width = 4, height = 4, startRow = 2)
      trainRow = 24
    }

    if(!is.null(listfeatures)){
      addStyle(wb, sheet = "Training", titleStyle, rows = trainRow, cols = 1)
      writeData(wb, sheet = "Training", x = "List of Significant Genes", startRow = trainRow)
      writeData(wb, sheet = "Training", x = listfeatures, startRow = trainRow + 1)
      trainRow = trainRow + nrow(listfeatures) + 3
    }
    
    if(!is.null(train.obj) && !is.na(input$threshold2)){
      png(file = "survresponseprediction.png")
      superpc.predictionplot(train.obj, data, data, threshold = input$threshold2)            
      dev.off()
    }
    
    if(!is.null(train.obj) && !is.na(input$threshold2)){ 
      addStyle(wb, sheet = "Training", titleStyle, rows = trainRow, cols = 1)
      writeData(wb, sheet = "Training", x = "Response Prediction Row", startRow = trainRow)
      insertImage(wb, sheet = "Training", file = "survresponseprediction.png", width = 4, height = 4, startRow = trainRow + 1)
    }
    
    if(!is.null(cv.obj)){
      png(file = "cvcurves.png")
      superpc.plotcv(cv.obj) 
      dev.off()      
    }
    
    cvRow = 0
    if(!is.null(cv.obj)){
      addWorksheet(wb, sheetName = "Cross Validation")
      setColWidths(wb, "Cross Validation", cols= 1, widths = 30)    
      addStyle(wb, sheet = "Cross Validation", titleStyle, rows = 1, cols = 1)
      writeData(wb, sheet = "Cross Validation", x = "CV curves", startRow = 1)
      insertImage(wb, sheet = "Cross Validation",  file = "cvcurves.png", width = 4, height = 4, startRow = 2)
      cvRow = 24
    }
    
    if(!is.null(fit.redcv) && !is.null(data$censoring.status)){
      png(file = "reducedmodelcv.png")
      superpc.plotred.lrtest(fit.redcv) 
      dev.off()
    }
    
    if(!is.null(fit.redcv) && !is.null(data$censoring.status)){
      addStyle(wb, sheet = "Cross Validation", titleStyle, rows = cvRow, cols = 1)
      writeData(wb, sheet = "Cross Validation", x = "Reduced model CV", startRow = cvRow)
      insertImage(wb, sheet = "Cross Validation",  file = "reducedmodelcv.png", width = 4, height = 4, startRow = cvRow + 1)
    }
    

    if(!is.null(lrtest.obj.test)){
      png(file = "testError.png")
      superpc.plot.lrtest(lrtest.obj.test)
      axis(3, at = lrtest.obj.test$threshold, labels = lrtest.obj.test$num.features)
      dev.off()
    }
    
    testRow = 0
    if(!is.null(lrtest.obj.test)){
      addWorksheet(wb, sheetName = "Test Set Prediction")
      setColWidths(wb, "Test Set Prediction", cols= 1, widths = 30)    
      addStyle(wb, sheet = "Test Set Prediction", titleStyle, rows = 1, cols = 1)
      writeData(wb, sheet = "Test Set Prediction", x = "Test Error", startRow = 1)
      insertImage(wb, sheet = "Test Set Prediction",  file = "testError.png", width = 4, height = 4, startRow = 2)
      testRow = 24
    }

    
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
      addStyle(wb, sheet = "Test Set Prediction", titleStyle, rows = testRow, cols = 1)
      writeData(wb, sheet = "Test Set Prediction", x = "Prediction info", startRow = testRow)
      writeData(wb, sheet = "Test Set Prediction", x = predictMatrix, startRow = testRow + 1)
      testRow = testRow + 4
      
      writeData(wb, sheet = "Test Set Prediction", x = outcome$coeftable, startRow = testRow)
      testRow = testRow + 3
      writeData(wb, sheet = "Test Set Prediction", x = outcome$teststat.table, startRow = testRow)
      testRow = testRow + 3
    }

    if(!is.null(predictDiscrete) && !is.null(data$censoring.status)){
      png(file = "testresponseprediction.png")
      plot(survfit(Surv(testData$y,testData$censoring.status)~predictDiscrete$v.pred), col=2:3, xlab="time", ylab="Prob survival")  
      dev.off()
    }
    
    if(!is.null(predictDiscrete) && !is.null(data$censoring.status)){
      addStyle(wb, sheet = "Test Set Prediction", titleStyle, rows = testRow, cols = 1)
      writeData(wb, sheet = "Test Set Prediction", x = "Response Prediction Plot", startRow = testRow)
      insertImage(wb, sheet = "Test Set Prediction",  file = "testresponseprediction.png", width = 4, height = 4, startRow = testRow + 1)
      testRow = testRow + 24
    }

    
    if(!is.null(data) && !is.null(testData) && !is.null(competingPredictorFit) && !is.na(input$threshold2) && !is.na(input$princomp) && !is.null(data$censoring.status)){
      sample.labels = competingPredictorFit$sample.labels
      competing.predictors = competingPredictorFit$competing.predictors
      png(file = "rainbowplot.png")
      if(is.na(input$shrinkage)){
        superpc.rainbowplot(data, predict$v.pred, sample.labels, competing.predictors)
      }else{
        superpc.rainbowplot(data, as.numeric(fit.red$v.test), sample.labels, competing.predictors)
      }      
      dev.off()      
    }
    
    if(!is.null(data) && !is.null(testData) && !is.null(competingPredictorFit) && !is.na(input$threshold2) && !is.na(input$princomp) && !is.null(data$censoring.status)){
      addStyle(wb, sheet = "Test Set Prediction", titleStyle, rows = testRow, cols = 1)
      writeData(wb, sheet = "Test Set Prediction", x = "Rainbow Plot", startRow = testRow)
      insertImage(wb, sheet = "Test Set Prediction",  file = "rainbowplot.png", width = 4, height = 4, startRow = testRow + 1)
    }
    
    if(!is.null(train.obj)){
      fname = paste(file, "xlsx", sep = ".")
      saveWorkbook(wb, file.path(dir, fname), overwrite = TRUE)
          
      if(Sys.info()[['sysname']] == "Windows"){
        shell(file.path(dir, fname))  
      } else if(Sys.info()[['sysname']] == "Darwin"){
        system(paste("open", fname))
      }
    }
        
  })
  }
  })  
  

  
})