library(shiny)

shinyUI(fluidPage(
  
  titlePanel("PAM - Prediction Analysis of Microarrays"),
  
  fluidRow(
    column(3,        
    wellPanel(
      
    fileInput(inputId = "iFile", label = "", accept="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"), 
    selectInput("analysisType", "", c("Classification" = "Classification", "Survival" = "Survival", "Regression" = "Regression")),
    conditionalPanel(condition = "input.analysisType == 'Classification'",
      numericInput("classLabel", label = "Class labels row", min = 1, max = 5, value = NULL, step = 1)            
    ),
    
    conditionalPanel(condition = "input.analysisType == 'Survival'",
      numericInput("survivalTimeLabel", label = "Survival Times row", min = 1, max = 5, value = NULL, step = 1),
      numericInput("censoringStatusLabel", label = "Censoring Status row", min = 1, max = 5, value = NULL, step = 1)
    ),
    
    conditionalPanel(condition = "input.analysisType == 'Regression'",
      numericInput("outcomeValueLabel", label = "Outcome Values row", min = 1, max = 5, value = NULL, step = 1)           
    ),
    
    numericInput("sampleLabel", label = "Sample labels row", min = 1, max = 5, value = NULL, step = 1),
    numericInput("batchLabel", label = "Batch labels row", min = 1, max = 5, value = NULL, step = 1),
    numericInput("expressionStart", label = "Expression data row", min = 1, max = 5, value = NULL, step = 1)
    ),
    
    wellPanel(
    conditionalPanel(condition = "input.analysisType == 'Classification'",
      uiOutput(outputId = "threshold"),
      numericInput("s0percentile", label = "Std. Dev. Factor S0 percentile (0-100)", min = 0, max = 100, value = 50, step = 0.0001),
      radioButtons("sign.contrast", "Contrast Sign", c("Both" = "both", "Positive" = "positive", "Negative" = "negative")),
      radioButtons("classPrior", "Class Prior", c("Sample Prior" = "sampleprior", "Uniform Prior" = "uniformprior", "Custom Prior" = "customprior")),
      conditionalPanel(condition = "input.classPrior == 'customprior'",  
                       uiOutput(outputId = "customvalue")
      )
    ),
    conditionalPanel(condition = "input.analysisType == 'Survival' || input.analysisType == 'Regression'",
      uiOutput(outputId = "threshold2")
    ),
    
    conditionalPanel(condition = "input.analysisType == 'Survival' || input.analysisType == 'Regression'",
      numericInput("princomp", label = "Princ Comp number for gene scores", min = 1, max = 3, value = NULL, step = 1),               
      uiOutput(outputId = "shrinkage")
    ),
   
    
    numericInput("randomSeed", label = "Random Number Generator Seed", min = 0, max = 1000000, value = 420473, step = 1),
    checkboxInput("cuberoot", label = "Transform by cube root?", value = FALSE),
    checkboxInput("center", label = "Center columns?", value = FALSE),
    checkboxInput("scale", label = "Scale columns?", value = FALSE),
    numericInput("numberOfNeighbors", "K-Nearest Neighbors Imputer: Number of Neighbors", value= 10, step=1)
    
    )
      
    ),
    
    column(9,
    
    conditionalPanel(condition = "input.analysisType == 'Classification'",       
      tabsetPanel(id = "PAM",
        tabPanel("Data", h3(textOutput("originalDataText")), dataTableOutput("dat"), h3(textOutput("imputedDataText")), dataTableOutput("imputedX"), h3(textOutput("transformDataText")), dataTableOutput("transform") ),          
        tabPanel("Training", h3(textOutput("trainErrorPlotText")), plotOutput("plotTrainError"), h3(textOutput("confusionTrainText")), tableOutput("pamrConfusionTrain"), h3(textOutput("listgenesText")), tableOutput("listgenes"), h3(textOutput("centroidText")),plotOutput("plotcen"), h3(textOutput("fdrText")), tableOutput("fdr"), h3(textOutput("fdrPlotText")), plotOutput("fdrPlot")),
        tabPanel("Cross Validation", h3(textOutput("overallText")), plotOutput("plotcv"), h3(textOutput("individualText")), plotOutput("plotcv2"), h3(textOutput("plotcvText")), plotOutput("plotcvprob"), h3(textOutput("cvConfusionMatrix")), tableOutput("pamrConfusion")),
        tabPanel("Test Set Prediction", fileInput(inputId = "testFile", label = "Test Set", accept="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"), h3(textOutput("testErrorPlotText")), plotOutput("plotTestError"), h3(textOutput("predictPlotText")), plotOutput("plotpredprob"), h3(textOutput("predictTableText")), tableOutput("predict") ),
        tabPanel("Settings",  h3(textOutput("settingsText")), tableOutput("settings"), h3(textOutput("settingsPriorText")), tableOutput("settingsPrior"))
      )
    ),
    
    conditionalPanel(condition = "input.analysisType == 'Survival' || input.analysisType == 'Regression'",
      tabsetPanel(id = "PAMSurv",
        tabPanel("Data", h3(textOutput("originalXText")), dataTableOutput("survdata"), h3(textOutput("imputedXSurvText")), dataTableOutput("imputedXSurv"), h3(textOutput("transformSurvText")), dataTableOutput("transformSurv"), h3(textOutput("decorrelateXText")), dataTableOutput("decorrelateX")),
        tabPanel("Training", fileInput(inputId = "competingPredictorFile", label = "Competing Predictors", accept="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"), h3(textOutput("survTrainErrorText")), plotOutput("plotLrtest"), h3(textOutput("listSurvGenesText")), tableOutput("listfeatures"), h3(textOutput("responsePredictionText")), plotOutput("survPredictionPlot")),
        tabPanel("Cross Validation", h3(textOutput("plotCvSurvText")), plotOutput("plotcvsurv"), h3(textOutput("plotredLrtestText")), plotOutput("plotredLrtest")),
        tabPanel("Test Set Prediction", fileInput(inputId = "testFileSurv", label = "Test Set", accept="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"), fileInput(inputId = "competingPredictorFitFile", label = "Fit with Competing Predictors", accept="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"), h3(textOutput("lrtestObjTestText")), plotOutput("plotLrtestTest"), h3(textOutput("predictionInfoText")), tableOutput("predictionscore"), tableOutput("coeftable"), tableOutput("teststatTable"), h3(textOutput("responsePredictionPlotText")), plotOutput("responsePredictionPlot"), h3(textOutput("rainbowPlotText")), plotOutput("rainbowPlot"))        
      )        
    )
     

    
    
    )
  )
))