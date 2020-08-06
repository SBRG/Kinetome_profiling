
load_env <- new.env()
source("analyses/kapp_analysis/ML/model_params/train_parms.r" , local = load_env)

# parameters

# created in "prepare_kappmax_labels_for_ML.R"
y_files <- c("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_Y.rds",
             "output_data/kappmax/Davidi_kappmax/davidi_kappmax_Y.rds",
             "output_data/kappmax/KO_MFA_Davidi_kappmax/KO_mfa_davidi_kappmax_Y.rds",
             "output_data/kcat_in_vitro/kcat_in_vitro_ML_Y.rds"
             )

names(y_files) <- c("kappmax_KO_ALE_per_pp_per_s",
                    "kappmax_davidi_per_pp_per_s",
                    "kappmax_KO_ALE_davidi_per_pp_per_s",
                    "kcat_iv_ML_per_AS_per_s"
                    ) # names should correspond to y variable name
train_ctrl <- load_env$tr_ctrl
hyper_parms <- load_env$hyper_parms
set.seed(load_env$pipeline_seed)

# for predictions
mods_in_ensemble <-  c( "train_enet",
                        "train_enet_i","train_rf",  
                        "train_rf_i","best_model_dl",   
                        "best_model_dl_i" )
feature_data_per_mod <- rep( "train_dat_mmc_full" , length (mods_in_ensemble) )
names(feature_data_per_mod) <- mods_in_ensemble
feature_data_per_mod[ grepl("_i$",mods_in_ensemble) ] <- "train_dat_i_mmc_full"
h2o_mod_log <- grepl("_dl",mods_in_ensemble)

#####################################################
    
# prepare train set:
# input: 
# -  destination of file with features
# - model obect rdata file that is used to remove spontaneous reactions
# -  destination of kappmax df with reaction-specific bigg IDs
#
# output:
# - train data (still contains reactions with NA labels, but no spontaneous or exchange reactions) 
source("analyses/auxilliary_scripts/ML/prep_train_data.r")
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
source("analyses/auxilliary_scripts/ML/train_ensemble.r")
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
source("analyses/auxilliary_scripts/ML/predict_from_ensemble.r")

#####################################################

# output structures
train_dat <- vector(mode = "list", length = length(y_files) )
train_res <- vector(mode = "list", length = length(y_files) )
preds <- vector(mode = "list", length = length(y_files) )


names(train_dat) <- names(train_res) <- names(preds) <- names(y_files)

#####################################################
# i <- 1
for ( i in seq_along(y_files) ) { 
    print(paste("processing input file ", i) )

    train_dat[[i]] <- prep_train_data( 
                     feature_file = "data/ML_inputs/ML_features.rdata",
                     model_file = "data/M-models/iML1515/iML1515.rdata" ,
                     output_file = y_files[i]
                     )

    train_res[[i]] <- train_ensemble( 
                    train_dat = train_dat[[i]] , 
                    train_ctrl  = train_ctrl, 
                    hyper_parms = hyper_parms,
                    output_var_name = names(y_files)[i] ,
                    train_mice = FALSE ,
                    rm_nonopt_mods = FALSE,
                    use_testing = FALSE,
                    test_frac = 0.2,
                    shutdown_h2o_clust = FALSE
    )
}    


preds <- predict_from_ensemble(train_res = train_res,
                               mods_in_ensemble = mods_in_ensemble, #char
                               feature_data_per_mod = feature_data_per_mod, # named char vector
                               reacts_to_pred = rownames(train_res$kappmax_KO_ALE$full_data$train_dat_full) , 
                               h2o_mod_log = h2o_mod_log
)


# there shouldn't be any NAs in the predictions:
stopifnot( ! any(is.na(preds)) )

head(preds$ens_pred_df)
# check preformances
lapply( train_res , function(x){x$perform_stats} )

saveRDS( list(preds = preds , 
              train_res = train_res ,
              pars = list (hyper_parms = hyper_parms,
                           train_ctrl=train_ctrl,
                           pipeline_seed = load_env$pipeline_seed,
                           y_files = y_files,
                           mods_in_ensemble = mods_in_ensemble
                           )
          ) , 
          file = paste0("analyses/kapp_analysis/ML/saved_predictions/kappmax_preds_and_models_",
                        hyper_parms$h2o_search_criteria$max_models,"_dl_mods_",
                        hyper_parms$h2o_search_criteria$max_runtime_secs,"_s_max_runtime"
                        ,"_4x_cutoff.rds" )
)

# continue in "prepare_kappmax_ML_preds_for_mech_models.R"
        