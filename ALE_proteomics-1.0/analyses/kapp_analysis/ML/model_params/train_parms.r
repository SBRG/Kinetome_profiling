require(caret)
pipeline_seed <- 7321

tr_ctrl = trainControl(
      method = "repeatedcv" ,
      #"LOOCV"
      number = 5 ,
      # K
      repeats = 5,
      savePredictions = "final"
)
    
hyper_parms <- list(    
    rf_parms = list(ntrees = 500 , tuneLength = 10) ,
    h2o_hyper_params = list(
      activation = c("Rectifier", "Tanh", "Maxout"),
      hidden = list(c(20, 20), c(50, 50), c(30, 30, 30), c(25, 25, 25, 25)),
      input_dropout_ratio = c(0, 0.05, 0.1, 0.2, 0.5),
      l1 = c(0, 1e-4,1e-5, 1e-6),
      l2 = c(0, 1e-4,1e-5, 1e-6),
      rho = c( 0.99),
      epsilon = c(1e-09 ,1e-08 ,1e-07 ),
      rate = seq(0.0001, 0.001 , length.out = 10),
      epochs = 1000,
      stopping_rounds=2L,
      stopping_metric="MSE", ## could be "MSE","logloss","r2"
      stopping_tolerance= 0.01
    ) ,
    h2o_search_criteria = list(
      strategy = "RandomDiscrete",
      max_runtime_secs = 3600*12/(3*4),
      max_models = 5000,
      seed = 4
    )
  )
