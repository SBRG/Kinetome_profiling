# predict from ensemble
# input: 
# - train_res: list of model outputs containins lists of model objects and train data
# - mods_in_ensemble names of models as they are found in train_res
# - feature_data_per_mod, # named char vector
# reacts_to_pred, # reaction ids of reactions to predict
# h2o_mod_log
# output:
# - long format df with kappmax predictions with direction-specific bigg IDs for all outputs and models

# the current version expects a h2o cluster to be running

predict_from_ensemble <- function(train_res,
                                  mods_in_ensemble, #char
                                  feature_data_per_mod, # named char vector
                                  reacts_to_pred, # reaction ids of reactions to predict
                                  h2o_mod_log # logical indicating which models are based on h2o
                                  ) {
  pred_li <- vector(  mode = "list" , length = length(train_res) )
  names(pred_li) <- names(train_res)
  # init output df w results in long format
  indiv_pred_df <- data.frame(keff_source = character() , 
                        model = character() , 
                        react_id  = character() , 
                        kappmax_per_pp_per_s  = numeric() , 
                        stringsAsFactors = FALSE)
  
  for (i in  seq_along(pred_li) ) {
    stopifnot(all(mods_in_ensemble %in% names(train_res[[i]]$models) ) )
    for(j in seq_along(mods_in_ensemble) ){
        mod_j <- mods_in_ensemble[j]
        dataset_j <- feature_data_per_mod[ names(feature_data_per_mod) == mod_j ]
        newdata <- train_res[[i]]$full_data[[ dataset_j ]]
        # remove cases where features contain NAs
        obs_w_compl_feats <- apply( newdata[ , colnames(newdata)!=names(pred_li)[i] ] , 1 , function(x) ! any(is.na(x)) )
        react_id_miss_obs <- rownames(newdata)[ ! obs_w_compl_feats ]
        newdata_compl_feats <- newdata[ obs_w_compl_feats , ]
        if(h2o_mod_log[j]){newdata_compl_feats_rownames <- rownames(newdata_compl_feats) # need to save those because h2o is not saving them
                           newdata_compl_feats <- as.h2o(newdata_compl_feats) }
        pred_j <- predict(  object = train_res[[i]]$models[[ mod_j ]] ,
                            newdata = newdata_compl_feats
        )
        if(h2o_mod_log[j]){pred_j <- as.data.frame(pred_j)$predict
                           names(pred_j) <- newdata_compl_feats_rownames 
        }
        # add NAs for reactions that did now get a prediction because obs were missing (do this to allow median replacement lateron)
        pred_j[ react_id_miss_obs ] <- NA
        indiv_pred_df <- rbind(indiv_pred_df , data.frame( 
                                    keff_source = names(pred_li)[i],
                                    model = mod_j , 
                                    react_id  = names(pred_j) , 
                                    kappmax_per_pp_per_s = 10^pred_j ,
                                    stringsAsFactors = FALSE ) 
        )
    }
  }
  library(dplyr)
  
  ens_pred_df <- indiv_pred_df %>% 
                      group_by( keff_source, model ) %>% # add a median-imputed version
                      mutate( kappmax_per_pp_per_s_med_imp_temp = ifelse(is.na(kappmax_per_pp_per_s) , 
                                                                    median(kappmax_per_pp_per_s , na.rm=TRUE), 
                                                                    kappmax_per_pp_per_s ) ) %>%  
                      group_by( keff_source , react_id ) %>% 
                      summarise( kappmax_per_pp_per_s_med_imp = mean(kappmax_per_pp_per_s_med_imp_temp) ,
                                 kappmax_per_pp_per_s_med_imp_ls = 10^ mean( log10(kappmax_per_pp_per_s_med_imp_temp) ) , 
                                 kappmax_per_pp_per_s = mean(kappmax_per_pp_per_s,na.rm=TRUE) ) %>% 
                      as.data.frame(  stringsAsFactors = FALSE)
  
  return(list( indiv_pred_df = indiv_pred_df , 
               ens_pred_df = ens_pred_df
               ) )
}


