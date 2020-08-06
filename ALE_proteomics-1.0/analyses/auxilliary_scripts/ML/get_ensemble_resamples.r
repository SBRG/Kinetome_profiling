get_ensemble_resamples <- function( ml_res , 
                         h2o_mod_log ,
                         mod_imputed_log #,
                         #epochs
                         
) {
  train_res <- ml_res$train_res
  
  resamp_df <- data.frame( RMSE = numeric(0) , 
                           Rsquared = numeric(0) , 
                           MAE = numeric(0) , 
                           Resample = character(0) , 
                           ML_model = character(0) , 
                           keff_source = character(0) ) 
  
  library(h2o)
  h2o.init()
  
  # go over outputs:
  for ( i in seq_along(train_res) ){
    # go over ML models
    for( j in seq_along( train_res[[i]]$model ) ){
      if(h2o_mod_log[j]){# h2o model? need to rerun CV in this case
          
          model_name <- paste0( names(train_res)[i] , ifelse(mod_imputed_log[j] , "_median_impu" , "") , "_best_model" )
        
          #train_res[[i]]$model[j]$best_model_dl@model_id
          mod <- h2o.loadModel( paste0( "analyses/kapp_analysis/ML/saved_models/h2o/dl_grid_",model_name,"/",train_res[[i]]$model[j]$best_model_dl@model_id ) )
          
          output_var_name <- names(train_res)[i]
          
          # unimputed
          if(! mod_imputed_log[j] ){
              #we need to log10 the output before creating a h2o object
              train_dat_mmc_df_h2o <- train_res[[i]]$train_data$train_dat_woi_mmc
          }else{
              train_dat_mmc_df_h2o <- train_res[[i]]$train_data$train_dat_i_mmc
          }
          
          train_dat_mmc_df_h2o[  , colnames(train_dat_mmc_df_h2o) == output_var_name ] <-  log10( train_dat_mmc_df_h2o[  , colnames(train_dat_mmc_df_h2o) == output_var_name ] )
          train_dat_mmc_df_h2o <-  as.h2o( train_dat_mmc_df_h2o[ sample( 1:nrow(train_dat_mmc_df_h2o) , nrow(train_dat_mmc_df_h2o) , replace = FALSE ) ,  ] )
          
          all_parms <- train_res[[i]]$model[j]$best_model_dl@allparameters
          
          all_parms$keep_cross_validation_predictions <- TRUE
          all_parms$training_frame <- train_dat_mmc_df_h2o
          
          #all_parms$epochs <- epochs
          
          all_parms_dl <- all_parms[ names(all_parms) %in% formalArgs(h2o.deeplearning) ] 
          
          for ( k in 1:5){ # repeat 5 times 
              refit <- do.call(  h2o.deeplearning , all_parms_dl ) 
              
              dl_resamp <- refit@model$cross_validation_metrics_summary
              
              dl_resamp_R2 <- as.numeric( unlist( dl_resamp[ "r2", grep("^cv" , colnames(dl_resamp) ) , drop = TRUE ] ) )
              dl_resamp_RMSE <- as.numeric( unlist( dl_resamp[ "rmse", grep("^cv" , colnames(dl_resamp) ), drop = TRUE  ] ) )
              
              resamp_df <- rbind( resamp_df , 
                                  data.frame( RMSE = dl_resamp_RMSE , Rsquared = dl_resamp_R2 , MAE = NA , 
                                      Resample = paste0(names(dl_resamp_R2) ,"_rep_" , k ) , ML_model = names(train_res[[i]]$model)[j] ,
                                      keff_source = names(train_res)[i] , stringsAsFactors = FALSE)
                                  )
          }
      }else{ # caret model: get resamples from train object
          resamp_j <- train_res[[i]]$model[[j]]$resample
          resamp_j$ML_model <-  names(train_res[[i]]$model)[j]
          resamp_j$keff_source <- names(train_res)[i]
          resamp_df <- rbind( resamp_df, resamp_j ) 
      }
    }
  }
  return( resamp_df ) 
}

ml_res <- readRDS("analyses/kapp_analysis/ML/saved_predictions/kappmax_preds_and_models_5000_dl_mods_3600_s_max_runtime_4x_cutoff.rds")


set.seed( 6324 )
resamps <- get_ensemble_resamples (ml_res,
                                   h2o_mod_log = grepl("_dl", names(ml_res$train_res[[1]]$models )),
                                   mod_imputed_log = grepl("_i$", names(ml_res$train_res[[1]]$models ))#,
                                  # epochs = 1e3
) 

library(dplyr)
resamps %>% group_by(ML_model , keff_source) %>% 
  summarise( m.RMSE = mean(RMSE), sd.RMSE = sd(RMSE),  m.Rsquared = mean(Rsquared) , m.sd = sd(Rsquared)) %>%  print(n=Inf)

resamps %>% filter(keff_source == "kappmax_davidi_per_pp_per_s") %>% group_by(ML_model , keff_source) %>% 
  summarise( m.RMSE = mean(RMSE), sd.RMSE = sd(RMSE),  m.Rsquared = mean(Rsquared) , m.sd = sd(Rsquared)) %>% print(n=Inf)

resamps %>% filter(ML_model == "best_model_dl" , keff_source == "kcat_iv_ML_per_AS_per_s")

save( resamps , file = "analyses/kapp_analysis/ML/saved_resamples/kappmax_preds_and_models_5000_dl_mods_3600_s_max_runtime_4x_cutoff_resamps_org_epochs.rdata" )
