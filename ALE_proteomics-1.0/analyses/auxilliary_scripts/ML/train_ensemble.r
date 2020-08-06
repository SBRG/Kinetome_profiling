# train ensemble:
# input:
# - train data
# - traincontrol
# - h2o grid
# - name of output variable
#
# output 
# - performance table
# - list of model objects
# - train data

train_ensemble <- function(train_dat,
                           train_ctrl, 
                           hyper_parms,
                           output_var_name,
                           use_testing,
                           test_frac ,
                           train_mice = FALSE ,
                           rm_nonopt_mods = FALSE, # can avoid memory issues 
                           shutdown_h2o_clust,
                           h2o.max_mem_size = "20G"
                           ) {
  require(caret)
  require(mice)
  require(h2o)
  
  ############################################
  # performance summary
  perform <- data.frame(model = NA,	CV.RMSE =NA,	CV.RMSE.sd =NA,	CV.R2 =NA,	CV.R2.SD =NA,	
                        CV.method =NA, test.RMSE =NA,	test.R2 =NA,		
                        n.full.dataset =NA,	Imputation =NA  )
  
  ############################################
  # remove incomplete cases , but also keep a full data set for export to the prediction function
  
  train_dat_full <- train_dat_i_full <- train_dat_i_mice_full <- train_dat
  
  # select complete cases 
  train_compl_log <- complete.cases( train_dat )
  # select complete output cases :
  output_compl_log <- ! is.na(train_dat[ , output_var_name , drop= TRUE ])
  
  print(paste("number of complete cases:",sum(train_compl_log)) )
  print(paste("number of complete output only cases:",sum(output_compl_log)) )
  
  # remove incomplete cases for the unimputed train set (woi = wo imputation)
  train_dat_woi <- train_dat[ train_compl_log , ]
  #imputed train set (i = w imputation)
  train_dat_i <- train_dat[ output_compl_log , ]
  
  ############################################
  # impute one of the sets with mice:
  
  num_cols <- sapply(train_dat_i,class)=="numeric"
  char_cols <- sapply(train_dat_i,class)=="character"
  out_col <- colnames(train_dat_i) == output_var_name
  
  if(train_mice){
        train_dat_i_mice <- train_dat_i
        
        
        # convert char to factor for mice imputation
        train_dat_i_mice[ , char_cols ] <- lapply(train_dat_i_mice[ , char_cols ] , as.factor)
        train_dat_i_mice_full[ , char_cols ] <- lapply(train_dat_i_mice_full[ , char_cols ] , as.factor)
        
        # exclude output from imputation
        predictorMatrix <- matrix( 1 , ncol(train_dat_i_mice) , ncol(train_dat_i_mice) )
        diag(predictorMatrix) <- 0
        predictorMatrix[ , out_col ] <- 0
        
        mice_imp <- mice(train_dat_i_mice , predictorMatrix = predictorMatrix)
        mice_dat <- mice::complete(mice_imp)
        
        mice_imp_full <- mice(train_dat_i_mice_full , predictorMatrix = predictorMatrix)
        mice_dat_full <- mice::complete(mice_imp_full)
        
        # convert factors back to character
        fact_cols <- sapply(mice_dat,class)=="factor"
        mice_dat[ , fact_cols ] <- lapply( mice_dat[ , fact_cols ] , as.character)
        rownames(mice_dat) <- rownames(train_dat_i_mice)
        
        fact_cols_full <- sapply(mice_dat_full,class)=="factor"
        mice_dat_full[ , fact_cols_full ] <- lapply( mice_dat_full[ , fact_cols_full ] , as.character)
        rownames(mice_dat_full) <- rownames(train_dat_i_mice_full)
        
        # use mice data 
        train_dat_i_mice <- mice_dat
        train_dat_i_mice_full <- mice_dat_full
        nrow(train_dat_i_mice)
        lapply(train_dat_i_mice,class)
  }
  
  ############################################
  # impute one of the sets with median:
  
  train_dat_i <- median_maj_impute( train_dat_i , output_var_name = output_var_name )
  stopifnot( ! any(is.na(train_dat_i[  , colnames(train_dat_i) != output_var_name ])) )
  train_dat_i_full <- median_maj_impute( train_dat_i_full , output_var_name = output_var_name )
  
  ############################################
  # create dummmys
  
  # caret dummies for non-lm models 
  dum_woi <- dummyVars(  ~ . , data = train_dat_woi )
  train_dat_woi_mmc <- predict( dum_woi , newdata = train_dat_woi )
  dum_woi_full <- dummyVars(  ~ . , data = train_dat_full )
  train_dat_mmc_full <- predict( dum_woi , newdata = train_dat_full )
  
  if(train_mice){
  # caret dummies for non-lm models 
  dum_i_mice <- dummyVars(  ~ . , data = train_dat_i_mice )
  train_dat_i_mice_mmc <- predict( dum_i_mice , newdata = train_dat_i_mice )
  dum_i_mice_full <- dummyVars(  ~ . , data = train_dat_i_mice_full )
  train_dat_i_mice_mmc_full <- predict( dum_i_mice_full , newdata = train_dat_i_mice_full )
  }
  
  # caret dummies for non-lm models 
  dum_i <- dummyVars(  ~ . , data = train_dat_i )
  train_dat_i_mmc <- predict( dum_i , newdata = train_dat_i )
  dum_i_full <- dummyVars(  ~ . , data = train_dat_i_full )
  train_dat_i_mmc_full <- predict( dum_i_full , newdata = train_dat_i_full )
  
  
  ############################################
  # create test set 
  if( use_testing ){
    train_index_woi <- sample( 1:nrow(train_dat_woi) , floor((1-test_frac)*nrow(train_dat_woi)) ) # un-stratified sampling
    
    train_dat_woi_test <- train_dat_woi[ - train_index_woi , ]
    train_dat_woi_train <- train_dat_woi[ train_index_woi , ]

    train_dat_woi_mm_train <- train_dat_woi_mm[ train_index_woi , ]
    
    train_dat_woi_mmc_test <- train_dat_woi_mmc[ - train_index_woi , ]
    train_dat_woi_mmc_train <- train_dat_woi_mmc[ train_index_woi , ]
    
    if(train_mice){
    train_index_i_mice <- sample( 1:nrow(train_dat_i_mice) , floor((1-test_frac)*nrow(train_dat_i_mice)) ) # un-stratified sampling
    
    # mice imputation casess
    train_dat_i_mice_test <- train_dat_i_mice[ - train_index_i_mice , ]
    train_dat_i_mice_train <- train_dat_i_mice[ train_index_i_mice , ]
    
    train_dat_i_mice_mmc_test <- train_dat_i_mice_mmc[ - train_index_i_mice , ]
    train_dat_i_mice_mmc_train <- train_dat_i_mice_mmc[ train_index_i_mice , ]
    }
    
    train_index_i <- sample( 1:nrow(train_dat_i) , floor((1-test_frac)*nrow(train_dat_i)) ) # un-stratified sampling
    
    # median imputation casess
    train_dat_i_test <- train_dat_i[ - train_index_i , ]
    train_dat_i_train <- train_dat_i[ train_index_i , ]

  }

  ############################################
  ############################################
  # linear model
  print(paste("lm for " , output_var_name) )
  # wo imputation
  train_lm <- train( x = train_dat_woi[  , colnames(train_dat_woi) != output_var_name ] , 
                     y = log10( train_dat_woi[  , colnames(train_dat_woi) == output_var_name , drop = TRUE ] ) , 
                     preProcess = c("center","scale"), 
                     method = "lm",
                     trControl = train_ctrl
  )
  # 
  if(use_testing){ # compute testing error 
    lm_test_RMSE <- caret::RMSE(pred = predict(train_lm , train_dat_woi_test[  , colnames(train_dat_woi_test) != output_var_name ]) , 
                                log10( train_dat_woi_test[  , colnames(train_dat_woi_test) == output_var_name , drop = TRUE ] ) )
    lm_test_R2 <- caret::R2(pred = predict(train_lm , train_dat_woi_test[  , colnames(train_dat_woi_test) != output_var_name ]) , 
                            log10( train_dat_woi_test[  , colnames(train_dat_woi_test) == output_var_name , drop = TRUE ] ))
  }else{  lm_test_RMSE <- NA;  lm_test_R2 <- NA
  }
  
  perform <-  rbind(perform, data.frame(model = train_lm$method,	CV.RMSE = train_lm$results$RMSE ,	CV.RMSE.sd = train_lm$results$RMSESD,	
                                        CV.R2 =train_lm$results$Rsquared,	CV.R2.SD =train_lm$results$RsquaredSD,	
                                        CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,	
                                        test.RMSE = lm_test_RMSE , test.R2 = lm_test_R2 ,
                                        n.full.dataset = nrow(train_lm$trainingData) ,	Imputation = "none")
  )
  
  perform<-perform[-1,]
  ############################################
  # lm  with mice imputation
  if(train_mice){
    train_lm_i_mice <- train( x = train_dat_i_mice[  , colnames(train_dat_i_mice) != output_var_name ] , 
                       y = log10( train_dat_i_mice[  , colnames(train_dat_i_mice) == output_var_name , drop = TRUE ] ) , 
                       preProcess = c("center","scale"), 
                       method = "lm",
                       trControl = train_ctrl
    )
    # 
    if(use_testing){ # compute testing error 
      lm_i_mice_test_RMSE <- caret::RMSE(pred = predict(train_lm_i_mice , train_dat_i_mice_test[  , colnames(train_dat_i_mice_test) != output_var_name ]) , 
                                  log10( train_dat_i_mice_test[  , colnames(train_dat_i_mice_test) == output_var_name , drop = TRUE ] ) )
      lm_i_mice_test_R2 <- caret::R2(pred = predict(train_lm_i_mice , train_dat_i_mice_test[  , colnames(train_dat_i_mice_test) != output_var_name ]) , 
                              log10( train_dat_i_mice_test[  , colnames(train_dat_i_mice_test) == output_var_name , drop = TRUE ] ))
    }else{  lm_i_mice_test_RMSE <- NA;  lm_i_mice_test_R2 <- NA
    }
    
    perform <-  rbind(perform, data.frame(model = train_lm_i_mice$method,	CV.RMSE = train_lm_i_mice$results$RMSE ,	CV.RMSE.sd = train_lm_i_mice$results$RMSESD,	
                                          CV.R2 =train_lm_i_mice$results$Rsquared,	CV.R2.SD =train_lm_i_mice$results$RsquaredSD,	
                                          CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,	
                                          test.RMSE = lm_i_mice_test_RMSE , test.R2 = lm_i_mice_test_R2 ,
                                          n.full.dataset = nrow(train_lm_i_mice$trainingData) ,	Imputation = "mice")
    )
  }
  ############################################
  # lm with median imputation
  train_lm_i <- train( x = train_dat_i[  , colnames(train_dat_i) != output_var_name ] , 
                       y = log10( train_dat_i[  , colnames(train_dat_i) == output_var_name , drop = TRUE ] ) , 
                       preProcess = c("center","scale"), 
                       method = "lm",
                       trControl = train_ctrl
  )
  # 
  if(use_testing){ # compute testing error 
    lm_i_test_RMSE <- caret::RMSE(pred = predict(train_lm_i , train_dat_i_test[  , colnames(train_dat_i_test) != output_var_name ]) , 
                                  log10( train_dat_i_test[  , colnames(train_dat_i_test) == output_var_name , drop = TRUE ] ) )
    lm_i_test_R2 <- caret::R2(pred = predict(train_lm_i , train_dat_i_test[  , colnames(train_dat_i_test) != output_var_name ]) , 
                              log10( train_dat_i_test[  , colnames(train_dat_i_test) == output_var_name , drop = TRUE ] ))
  }else{  lm_i_test_RMSE <- NA;  lm_i_test_R2 <- NA
  }
  
  perform <-  rbind(perform, data.frame(model = train_lm_i$method,	CV.RMSE = train_lm_i$results$RMSE ,	CV.RMSE.sd = train_lm_i$results$RMSESD,	
                                        CV.R2 =train_lm_i$results$Rsquared,	CV.R2.SD =train_lm_i$results$RsquaredSD,	
                                        CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,	
                                        test.RMSE = lm_i_test_RMSE , test.R2 = lm_i_test_R2 ,
                                        n.full.dataset = nrow(train_lm_i$trainingData) ,	Imputation = "median")
  )
  
  ############################################
  ############################################
  # enet 
  print(paste("enet for " , output_var_name) )
  # wo imputation
  train_enet <- train( x = train_dat_woi_mmc[  , colnames(train_dat_woi_mmc) != output_var_name ] , 
                     y = log10( train_dat_woi_mmc[  , colnames(train_dat_woi_mmc) == output_var_name , drop = TRUE ] ) , 
                     preProcess = c("center","scale"), 
                     tuneGrid = expand.grid( fraction = seq(0.3,0.8 , length.out = 5) , 
                                             lambda = seq(0,0.1 , length.out = 10) ),
                     method = "enet",
                     trControl = train_ctrl
  )
  # 
  if(use_testing){ # compute testing error 
    enet_test_RMSE <- caret::RMSE(pred = predict(train_enet , train_dat_woi_mmc_test[  , colnames(train_dat_woi_mmc_test) != output_var_name ]) , 
                                log10( train_dat_woi_mmc_test[  , colnames(train_dat_woi_mmc_test) == output_var_name , drop = TRUE ] ) )
    enet_test_R2 <- caret::R2(pred = predict(train_enet , train_dat_woi_mmc_test[  , colnames(train_dat_woi_mmc_test) != output_var_name ]) , 
                            log10( train_dat_woi_mmc_test[  , colnames(train_dat_woi_mmc_test) == output_var_name , drop = TRUE ] ))
  }else{  enet_test_RMSE <- NA;  enet_test_R2 <- NA
  }
  
  best.row_enet <- train_enet$results$lambda == train_enet$bestTune$lambda & train_enet$results$fraction == train_enet$bestTune$fraction
  
  perform <-  rbind(perform, data.frame(model = train_enet$method,	CV.RMSE = train_enet$results$RMSE[best.row_enet]  ,
                                        CV.RMSE.sd = train_enet$results$RMSESD[best.row_enet],	
                                        CV.R2 =train_enet$results$Rsquared[best.row_enet],	
                                        CV.R2.SD =train_enet$results$RsquaredSD[best.row_enet],	
                                        CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,
                                        test.RMSE = enet_test_RMSE , test.R2 = enet_test_R2 ,
                                        n.full.dataset = nrow(train_enet$trainingData) ,	
                                        Imputation = "none")
  )
  
  ############################################
  # enet  with mice imputation
  if(train_mice){
    train_enet_i_mice <- train( x = train_dat_i_mice_mmc[  , colnames(train_dat_i_mice_mmc) != output_var_name ] , 
                         y = log10( train_dat_i_mice_mmc[  , colnames(train_dat_i_mice_mmc) == output_var_name , drop = TRUE ] ) , 
                         preProcess = c("center","scale"), 
                         method = "enet",
                         trControl = train_ctrl
    )
    # 
    if(use_testing){ # compute testing error 
      enet_i_mice_test_RMSE <- caret::RMSE(pred = predict(train_enet_i_mice , train_dat_i_mice_mmc_test[  , colnames(train_dat_i_mice_mmc_test) != output_var_name ]) , 
                                         log10( train_dat_i_mice_mmc_test[  , colnames(train_dat_i_mice_mmc_test) == output_var_name , drop = TRUE ] ) )
      enet_i_mice_test_R2 <- caret::R2(pred = predict(train_enet_i_mice , train_dat_i_mice_mmc_test[  , colnames(train_dat_i_mice_mmc_test) != output_var_name ]) , 
                                     log10( train_dat_i_mice_mmc_test[  , colnames(train_dat_i_mice_mmc_test) == output_var_name , drop = TRUE ] ))
    }else{  enet_i_mice_test_RMSE <- NA;  enet_i_mice_test_R2 <- NA
    }
    
    best.row_enet_i_mice <- train_enet_i_mice$results$lambda == train_enet_i_mice$bestTune$lambda & train_enet_i_mice$results$fraction == train_enet_i_mice$bestTune$fraction
    
    perform <-  rbind(perform, data.frame(model = train_enet_i_mice$method,	CV.RMSE = train_enet_i_mice$results$RMSE[best.row_enet_i_mice]  ,
                                          CV.RMSE.sd = train_enet_i_mice$results$RMSESD[best.row_enet_i_mice],	
                                          CV.R2 =train_enet_i_mice$results$Rsquared[best.row_enet_i_mice],	
                                          CV.R2.SD =train_enet_i_mice$results$RsquaredSD[best.row_enet_i_mice],	
                                          CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,
                                          test.RMSE = enet_test_RMSE , test.R2 = enet_test_R2 ,
                                          n.full.dataset = nrow(train_enet_i_mice$trainingData) ,	
                                          Imputation = "mice")
    )
  }
  ############################################
  # enet with median imputation
  train_enet_i <- train( x = train_dat_i_mmc[  , colnames(train_dat_i_mmc) != output_var_name ] , 
                       y = log10( train_dat_i_mmc[  , colnames(train_dat_i_mmc) == output_var_name , drop = TRUE ] ) , 
                       preProcess = c("center","scale"), 
                       method = "enet",
                       trControl = train_ctrl
  )
  # 
  if(use_testing){ # compute testing error 
    enet_i_test_RMSE <- caret::RMSE(pred = predict(train_enet_i , train_dat_i_test[  , colnames(train_dat_i_test) != output_var_name ]) , 
                                  log10( train_dat_i_test[  , colnames(train_dat_i_test) == output_var_name , drop = TRUE ] ) )
    enet_i_test_R2 <- caret::R2(pred = predict(train_enet_i , train_dat_i_test[  , colnames(train_dat_i_test) != output_var_name ]) , 
                              log10( train_dat_i_test[  , colnames(train_dat_i_test) == output_var_name , drop = TRUE ] ))
  }else{  enet_i_test_RMSE <- NA;  enet_i_test_R2 <- NA
  }
  
  best.row_enet_i <- train_enet_i$results$lambda == train_enet_i$bestTune$lambda & train_enet_i$results$fraction == train_enet_i$bestTune$fraction
  
  perform <-  rbind(perform, data.frame(model = train_enet_i$method,	CV.RMSE = train_enet_i$results$RMSE[best.row_enet_i]  ,
                                        CV.RMSE.sd = train_enet_i$results$RMSESD[best.row_enet_i],	
                                        CV.R2 =train_enet_i$results$Rsquared[best.row_enet_i],	
                                        CV.R2.SD =train_enet_i$results$RsquaredSD[best.row_enet_i],	
                                        CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,
                                        test.RMSE = enet_test_RMSE , test.R2 = enet_test_R2 ,
                                        n.full.dataset = nrow(train_enet_i$trainingData) ,	
                                        Imputation = "median")
  )
  ############################################
  # rf 
  # wo imputation
  
  print(paste("rf for " , output_var_name) )
  
  train_rf <- train( x = train_dat_woi_mmc[  , colnames(train_dat_woi_mmc) != output_var_name ] , 
                       y = log10( train_dat_woi_mmc[  , colnames(train_dat_woi_mmc) == output_var_name , drop = TRUE ] ) , 
                       preProcess = c("center","scale"), 
                       tuneLength = hyper_parms$rf_parms$tuneLength ,
                       method = "rf",
                       ntree=hyper_parms$rf_parms$ntrees,
                       trControl = train_ctrl
  )
  # 
  if(use_testing){ # compute testing error 
    rf_test_RMSE <- caret::RMSE(pred = predict(train_rf , train_dat_woi_mmc_test[  , colnames(train_dat_woi_mmc_test) != output_var_name ]) , 
                                  log10( train_dat_woi_mmc_test[  , colnames(train_dat_woi_mmc_test) == output_var_name , drop = TRUE ] ) )
    rf_test_R2 <- caret::R2(pred = predict(train_rf , train_dat_woi_mmc_test[  , colnames(train_dat_woi_mmc_test) != output_var_name ]) , 
                              log10( train_dat_woi_mmc_test[  , colnames(train_dat_woi_mmc_test) == output_var_name , drop = TRUE ] ))
  }else{  rf_test_RMSE <- NA;  rf_test_R2 <- NA
  }
  
  best.row_rf <- train_rf$results$mtry == train_rf$bestTune$mtry
  
  perform <-  rbind(perform, data.frame(model = train_rf$method,	CV.RMSE = train_rf$results$RMSE[best.row_rf]  ,
                                        CV.RMSE.sd = train_rf$results$RMSESD[best.row_rf],	
                                        CV.R2 =train_rf$results$Rsquared[best.row_rf],	
                                        CV.R2.SD =train_rf$results$RsquaredSD[best.row_rf],	
                                        CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,
                                        test.RMSE = rf_test_RMSE , test.R2 = rf_test_R2 ,
                                        n.full.dataset = nrow(train_rf$trainingData) ,	
                                        Imputation = "none")
  )
  
  ############################################
  # rf  with mice imputation
  if(train_mice){
    train_rf_i_mice <- train( x = train_dat_i_mice_mmc[  , colnames(train_dat_i_mice_mmc) != output_var_name ] , 
                                y = log10( train_dat_i_mice_mmc[  , colnames(train_dat_i_mice_mmc) == output_var_name , drop = TRUE ] ) , 
                                preProcess = c("center","scale"), 
                                tuneLength = hyper_parms$rf_parms$tuneLength ,
                                method = "rf",
                                ntree=hyper_parms$rf_parms$ntrees,
                                trControl = train_ctrl
    )
    # 
    if(use_testing){ # compute testing error 
      rf_i_mice_test_RMSE <- caret::RMSE(pred = predict(train_rf_i_mice , train_dat_i_mice_mmc_test[  , colnames(train_dat_i_mice_mmc_test) != output_var_name ]) , 
                                           log10( train_dat_i_mice_mmc_test[  , colnames(train_dat_i_mice_mmc_test) == output_var_name , drop = TRUE ] ) )
      rf_i_mice_test_R2 <- caret::R2(pred = predict(train_rf_i_mice , train_dat_i_mice_mmc_test[  , colnames(train_dat_i_mice_mmc_test) != output_var_name ]) , 
                                       log10( train_dat_i_mice_mmc_test[  , colnames(train_dat_i_mice_mmc_test) == output_var_name , drop = TRUE ] ))
    }else{  rf_i_mice_test_RMSE <- NA;  rf_i_mice_test_R2 <- NA
    }
    
    best.row_rf_i_mice <- train_rf_i_mice$results$mtry == train_rf_i_mice$bestTune$mtry
    
    perform <-  rbind(perform, data.frame(model = train_rf_i_mice$method,	CV.RMSE = train_rf_i_mice$results$RMSE[best.row_rf_i_mice]  ,
                                          CV.RMSE.sd = train_rf_i_mice$results$RMSESD[best.row_rf_i_mice],	
                                          CV.R2 =train_rf_i_mice$results$Rsquared[best.row_rf_i_mice],	
                                          CV.R2.SD =train_rf_i_mice$results$RsquaredSD[best.row_rf_i_mice],	
                                          CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,
                                          test.RMSE = rf_test_RMSE , test.R2 = rf_test_R2 ,
                                          n.full.dataset = nrow(train_rf_i_mice$trainingData) ,	
                                          Imputation = "mice")
    )
  }
  ############################################
  # rf with median imputation
  train_rf_i <- train( x = train_dat_i_mmc[  , colnames(train_dat_i_mmc) != output_var_name ] , 
                         y = log10( train_dat_i_mmc[  , colnames(train_dat_i_mmc) == output_var_name , drop = TRUE ] ) , 
                         preProcess = c("center","scale"), 
                         tuneLength = hyper_parms$rf_parms$tuneLength ,
                         method = "rf",
                         ntree=hyper_parms$rf_parms$ntrees,
                         trControl = train_ctrl
  )
  # 
  if(use_testing){ # compute testing error 
    rf_i_test_RMSE <- caret::RMSE(pred = predict(train_rf_i , train_dat_i_test[  , colnames(train_dat_i_test) != output_var_name ]) , 
                                    log10( train_dat_i_test[  , colnames(train_dat_i_test) == output_var_name , drop = TRUE ] ) )
    rf_i_test_R2 <- caret::R2(pred = predict(train_rf_i , train_dat_i_test[  , colnames(train_dat_i_test) != output_var_name ]) , 
                                log10( train_dat_i_test[  , colnames(train_dat_i_test) == output_var_name , drop = TRUE ] ))
  }else{  rf_i_test_RMSE <- NA;  rf_i_test_R2 <- NA
  }
  
  best.row_rf_i <- train_rf_i$results$mtry == train_rf_i$bestTune$mtry
  
  perform <-  rbind(perform, data.frame(model = train_rf_i$method,	CV.RMSE = train_rf_i$results$RMSE[best.row_rf_i]  ,
                                        CV.RMSE.sd = train_rf_i$results$RMSESD[best.row_rf_i],	
                                        CV.R2 =train_rf_i$results$Rsquared[best.row_rf_i],	
                                        CV.R2.SD =train_rf_i$results$RsquaredSD[best.row_rf_i],	
                                        CV.method =paste0(  train_ctrl$repeats ,"x",train_ctrl$number,"fold" ) ,
                                        test.RMSE = rf_test_RMSE , test.R2 = rf_test_R2 ,
                                        n.full.dataset = nrow(train_rf_i$trainingData) ,	
                                        Imputation = "median")
  )
  
  ############################################
  ############################################
  # deep  learning wo imputation
  print(paste("dl for " , output_var_name))
  
  h2o.init(max_mem_size = h2o.max_mem_size)
  
  grid_id <-paste0("dl_grid_",output_var_name)
  
  #we need to log10 the output before creating a h2o object
  train_dat_woi_mmc_df_h2o <- train_dat_woi_mmc
  train_dat_woi_mmc_df_h2o[  , colnames(train_dat_woi_mmc_df_h2o) == output_var_name ] <-  log10( train_dat_woi_mmc_df_h2o[  , colnames(train_dat_woi_mmc_df_h2o) == output_var_name ] )
  
  
  train_dat_woi_mmc_h2o <-  as.h2o( train_dat_woi_mmc_df_h2o[ sample( 1:nrow(train_dat_woi_mmc_df_h2o) , nrow(train_dat_woi_mmc_df_h2o) , replace = FALSE ) ,  ] )
  
  train_dl <- h2o.grid( x =  setdiff(names(train_dat_woi_mmc_h2o), output_var_name) ,
                         y =  output_var_name ,
                         nfolds = 5, 
                         training_frame = train_dat_woi_mmc_h2o,
                         algorithm = "deeplearning" ,
                         grid_id = grid_id,
                         hyper_params = hyper_parms$h2o_hyper_params ,
                         search_criteria = hyper_parms$h2o_search_criteria,
                         keep_cross_validation_predictions = FALSE,
                        keep_cross_validation_models = FALSE,
                        keep_cross_validation_fold_assignment = FALSE
  )
  
  
  
  grid <- h2o.getGrid( grid_id ,sort_by="RMSE",decreasing=FALSE)
  
  best_model_dl <- h2o.getModel(grid@model_ids[[1]]) ## model with lowest logloss
  
  # save the best model, 
  h2o.saveModel(best_model_dl, path = paste0( "analyses/kapp_analysis/ML/saved_models/h2o/" , grid_id , "_best_model") , force = TRUE)
  # rm non-optimal models to save memory and envoke gc
  if(rm_nonopt_mods){
    h2o.rm( unlist(grid@model_ids)[-1] )
  }
  gc()
  
  ##########################################
  
  h2o.cv.summary <- best_model_dl@model$cross_validation_metrics_summary
  
  if(use_testing){ # compute testing error 
    dl_test_RMSE <- caret::RMSE(pred = as.data.frame(predict(best_model_dl , as.h2o(train_dat_woi_mmc_test) ))$predict , 
                                obs = train_dat_woi_mmc_test[ , colnames(train_dat_woi_mmc_test) == output_var_name ] )
    dl_test_R2 <- caret::R2(pred = as.data.frame(predict(best_model_dl , as.h2o(train_dat_woi_mmc_test) ))$predict , 
                            obs = train_dat_woi_mmc_test[ , colnames(train_dat_woi_mmc_test) == output_var_name ] )
  }else{  dl_test_RMSE <- NA;  dl_test_R2 <- NA
  }
  
  
  perform <-  rbind(perform, data.frame(model = best_model_dl@algorithm,	
                                        CV.RMSE = h2o.cv.summary["rmse","mean"]  ,
                                        CV.RMSE.sd = h2o.cv.summary["rmse","sd"],	
                                        CV.R2 =h2o.cv.summary["r2","mean"],	
                                        CV.R2.SD =h2o.cv.summary["r2","sd"],	
                                        CV.method =paste0(  1 ,"x",best_model_dl@parameters$nfolds,"fold" ) ,
                                        test.RMSE = dl_test_RMSE , 
                                        test.R2 = dl_test_R2, 
                                        n.full.dataset = nrow(train_dat_woi_mmc_df_h2o) ,	
                                        Imputation = "none")
  )
  ############################################
  # dl with mice imputation
  if(train_mice){
    grid_id <-paste0("dl_grid_",output_var_name,"_mice_impu")
    
    #we need to log10 the output before creating a h2o object
    train_dat_i_mice_mmc_df_h2o <- train_dat_i_mice_mmc
    train_dat_i_mice_mmc_df_h2o[  , colnames(train_dat_i_mice_mmc_df_h2o) == output_var_name ] <-  log10( train_dat_i_mice_mmc_df_h2o[  , colnames(train_dat_i_mice_mmc_df_h2o) == output_var_name ] )
    
    
    train_dat_i_mice_mmc_h2o <-  as.h2o( train_dat_i_mice_mmc_df_h2o[ sample( 1:nrow(train_dat_i_mice_mmc_df_h2o) , nrow(train_dat_i_mice_mmc_df_h2o) , replace = FALSE ) ,  ] )
    
    train_dl <- h2o.grid( x =  setdiff(names(train_dat_i_mice_mmc_h2o), output_var_name) ,
                          y =  output_var_name ,
                          nfolds = 5, 
                          training_frame = train_dat_i_mice_mmc_h2o,
                          algorithm = "deeplearning" ,
                          grid_id = grid_id,
                          hyper_params = hyper_parms$h2o_hyper_params,
                          search_criteria = hyper_parms$h2o_search_criteria,
                          keep_cross_validation_predictions = FALSE,
                          keep_cross_validation_models = FALSE,
                          keep_cross_validation_fold_assignment = FALSE
    )
    
    
    
    grid <- h2o.getGrid( grid_id ,sort_by="RMSE",decreasing=FALSE)
    
    best_model_dl_i_mice <- h2o.getModel(grid@model_ids[[1]]) ## model with lowest logloss
    
    # save the best model, 
    h2o.saveModel(best_model_dl_i_mice, path = paste0( "analyses/kapp_analysis/ML/saved_models/h2o/" , grid_id , "_best_model")  , force = TRUE)
    # rm non-optimal models to save memory and envoke gc
    if(rm_nonopt_mods){
      h2o.rm( unlist(grid@model_ids)[-1] )
    }
    gc()
  
    ##########################################
  
    h2o.cv.summary <- best_model_dl_i_mice@model$cross_validation_metrics_summary
  
    
    if(use_testing){ # compute testing error 
      dl_test_RMSE <- caret::RMSE(pred = as.data.frame(predict(best_model_dl_i_mice , as.h2o(train_dat_i_mice_mmc_test) ))$predict , 
                                  obs = train_dat_i_mice_mmc_test[ , colnames(train_dat_i_mice_mmc_test) == output_var_name ] )
      dl_test_R2 <- caret::R2(pred = as.data.frame(predict(best_model_dl_i_mice , as.h2o(train_dat_i_mice_mmc_test) ))$predict , 
                              obs = train_dat_i_mice_mmc_test[ , colnames(train_dat_i_mice_mmc_test) == output_var_name ] )
    }else{  dl_test_RMSE <- NA;  dl_test_R2 <- NA
    }
    
    
    perform <-  rbind(perform, data.frame(model = best_model_dl_i_mice@algorithm,	
                                          CV.RMSE = h2o.cv.summary["rmse","mean"]  ,
                                          CV.RMSE.sd = h2o.cv.summary["rmse","sd"],	
                                          CV.R2 =h2o.cv.summary["r2","mean"],	
                                          CV.R2.SD =h2o.cv.summary["r2","sd"],	
                                          CV.method =paste0(  1 ,"x",best_model_dl_i_mice@parameters$nfolds,"fold" ) ,
                                          test.RMSE = dl_test_RMSE , 
                                          test.R2 = dl_test_R2, 
                                          n.full.dataset = nrow(train_dat_i_mice_mmc_df_h2o) ,	
                                          Imputation = "mice")
    )
  }
  ############################################
  # dl with median imputation
  
  grid_id <-paste0("dl_grid_",output_var_name,"_median_impu")
  
  #we need to log10 the output before creating a h2o object
  train_dat_i_mmc_df_h2o <- train_dat_i_mmc
  train_dat_i_mmc_df_h2o[  , colnames(train_dat_i_mmc_df_h2o) == output_var_name ] <-  log10( train_dat_i_mmc_df_h2o[  , colnames(train_dat_i_mmc_df_h2o) == output_var_name ] )
  
  
  train_dat_i_mmc_h2o <-  as.h2o( train_dat_i_mmc_df_h2o[ sample( 1:nrow(train_dat_i_mmc_df_h2o) , nrow(train_dat_i_mmc_df_h2o) , replace = FALSE ) ,  ] )
  
  train_dl <- h2o.grid( x =  setdiff(names(train_dat_i_mmc_h2o), output_var_name) ,
                        y =  output_var_name ,
                        nfolds = 5, 
                        training_frame = train_dat_i_mmc_h2o,
                        algorithm = "deeplearning" ,
                        grid_id = grid_id,
                        hyper_params = hyper_parms$h2o_hyper_params,
                        search_criteria = hyper_parms$h2o_search_criteria,
                        keep_cross_validation_predictions = FALSE,
                        keep_cross_validation_models = FALSE,
                        keep_cross_validation_fold_assignment = FALSE
  )
  
  
  grid <- h2o.getGrid( grid_id ,sort_by="RMSE",decreasing=FALSE)
  
  best_model_dl_i <- h2o.getModel(grid@model_ids[[1]]) ## model with lowest logloss
  
  # save the best model, 
  h2o.saveModel(best_model_dl_i, path = paste0( "analyses/kapp_analysis/ML/saved_models/h2o/" , grid_id , "_best_model") , force = TRUE)
  # rm non-optimal models to save memory and envoke gc
  if(rm_nonopt_mods){
    h2o.rm( unlist(grid@model_ids)[-1] )
  }
  gc()
  
  ##########################################
  
  h2o.cv.summary <- best_model_dl_i@model$cross_validation_metrics_summary
  
  if(use_testing){ # compute testing error 
    dl_test_RMSE <- caret::RMSE(pred = as.data.frame(predict(best_model_dl_i , as.h2o(train_dat_i_mmc_test) ))$predict , 
                                obs = train_dat_i_mmc_test[ , colnames(train_dat_i_mmc_test) == output_var_name ] )
    dl_test_R2 <- caret::R2(pred = as.data.frame(predict(best_model_dl_i , as.h2o(train_dat_i_mmc_test) ))$predict , 
                            obs = train_dat_i_mmc_test[ , colnames(train_dat_i_mmc_test) == output_var_name ] )
  }else{  dl_test_RMSE <- NA;  dl_test_R2 <- NA
  }
  
  
  perform <-  rbind(perform, data.frame(model = best_model_dl_i@algorithm,	
                                        CV.RMSE = h2o.cv.summary["rmse","mean"]  ,
                                        CV.RMSE.sd = h2o.cv.summary["rmse","sd"],	
                                        CV.R2 =h2o.cv.summary["r2","mean"],	
                                        CV.R2.SD =h2o.cv.summary["r2","sd"],	
                                        CV.method =paste0(  1 ,"x",best_model_dl_i@parameters$nfolds,"fold" ) ,
                                        test.RMSE = dl_test_RMSE , 
                                        test.R2 = dl_test_R2, 
                                        n.full.dataset = nrow(train_dat_i_mmc_df_h2o) ,	
                                        Imputation = "median")
  )
  
  
  ##########################################
  
  if(shutdown_h2o_clust){h2o.shutdown()}
  
  ##########################################
  ##########################################
  # output structures
  if(train_mice){ 
  out <- list(
        perform_stats = perform ,
        models = list( train_lm = train_lm , train_lm_i_mice = train_lm_i_mice , train_lm_i = train_lm_i ,
                       train_enet = train_enet , train_enet_i_mice = train_enet_i_mice , train_enet_i = train_enet_i , 
                       train_rf = train_rf , train_rf_i_mice = train_rf_i_mice , train_rf_i = train_rf_i ,
                       best_model_dl = best_model_dl , best_model_dl_i_mice = best_model_dl_i_mice , best_model_dl_i = best_model_dl_i
        ),
        train_data = list( train_dat_woi = train_dat_woi , train_dat_i_mice = train_dat_i_mice , train_dat_i = train_dat_i,
                           train_dat_woi_mmc = train_dat_woi_mmc , train_dat_i_mice_mmc = train_dat_i_mice_mmc , train_dat_i_mmc = train_dat_i_mmc
        ),
        full_data = list( train_dat_full = train_dat_full,
                          train_dat_mmc_full = train_dat_mmc_full,
                          train_dat_i_full = train_dat_i_full,
                          train_dat_i_mmc_full = train_dat_i_mmc_full,
                          train_dat_i_mice_full = train_dat_i_mice_full,
                          train_dat_i_mice_mmc_full = train_dat_i_mice_mmc_full
                          )
    )
  }else{
    out <- list(
      perform_stats = perform ,
      models = list( train_lm = train_lm ,train_lm_i = train_lm_i ,
                     train_enet = train_enet ,train_enet_i = train_enet_i , 
                     train_rf = train_rf ,  train_rf_i = train_rf_i ,
                     best_model_dl = best_model_dl , best_model_dl_i = best_model_dl_i
      ),
      train_data = list( train_dat_woi = train_dat_woi , train_dat_i = train_dat_i,
                         train_dat_woi_mmc = train_dat_woi_mmc , train_dat_i_mmc = train_dat_i_mmc
      ),
      full_data = list( train_dat_full = train_dat_full,
                        train_dat_mmc_full = train_dat_mmc_full,
                        train_dat_i_full = train_dat_i_full,
                        train_dat_i_mmc_full = train_dat_i_mmc_full
      )
    )
  }
  return( out )
}

# given a datframe X and the name of the output column (which will not be imputed) impute with median for neumnerics and 
# the majority for chars (factors are not touched)
median_maj_impute <- function(X , output_var_name){
  stopifnot(output_var_name %in% colnames(X))
  
  to_imp <- apply(X , 2, function(x) any(is.na(x)) )
  to_imp[ colnames(X) == output_var_name  ] <- FALSE
  
  
  num_cols <- sapply(X,is.numeric)
  char_cols <- sapply(X,is.character)
  out_col <- colnames(X) == output_var_name
  
  for( i in 1:ncol(X)){
    if( to_imp[i] & num_cols[i] ){  
      X[,i][is.na(X[,i])] <- median(X[,i] , na.rm = TRUE)
    }
    if( to_imp[i] & char_cols[i] ){  
      X[,i][is.na(X[,i])] <- names(sort(table(X[,i]),decreasing=TRUE))[1]
    }
  }
  return(X)
}

