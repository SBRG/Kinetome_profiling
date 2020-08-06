# prepare ensemble predictions for MOMENT and ME modelling
library(reshape2)
# parameters:
pred_to_use <- "kappmax_per_pp_per_s_med_imp"

######################################################################
# result of "kappmax_ML_pipeline.r"
ml_res <- readRDS("analyses/kapp_analysis/ML/saved_predictions/kappmax_preds_and_models_5000_dl_mods_3600_s_max_runtime_4x_cutoff.rds")
ens_preds <- as.data.frame(ml_res$preds$ens_pred_df ,  stringsAsFactors = FALSE)
# select how unpredicted (because of missing features) observations are handled
ens_preds <- ens_preds[ , c("keff_source"  ,   "react_id" , pred_to_use ) ]
sum(is.na(ens_preds))

######################################################################
# convert to  wide
ml_preds_w <- reshape2::dcast( ens_preds , 
                     react_id ~ keff_source  , 
                     value.var = pred_to_use)
# check NAs
any(is.na(ml_preds_w))
ml_preds_w[ apply(ml_preds_w , 1, function(x){any(is.na(x))} ) ,  ] 
numcols <- sapply(ml_preds_w , is.numeric)
head(ml_preds_w)
ml_preds_w[ , numcols  ] <-  apply( ml_preds_w[ , numcols  ] , 2, function(x){x[is.na(x)]<-median(x,na.rm = TRUE) ; return(x)} )
head(ml_preds_w)
any(is.na(ml_preds_w))

numcols_names <- names(ml_preds_w)[numcols]
######################################################################
# merge "measured" kappmaxes and kcat in vitro
for ( i in seq_along(numcols_names)){
  output_name <- numcols_names[i]
  train_dat <- ml_res$train_res[[output_name]]$train_data$train_dat_i
  print(sum(!is.na(train_dat[ , output_name ])))
  output <- train_dat[ , output_name , drop = FALSE]
  output$react_id <- rownames(output)
  ml_preds_w <- merge( ml_preds_w , output , by = "react_id" , all.x = TRUE , suffixes = c("_pred","_meas") )
}
  
head(ml_preds_w)
sapply(ml_preds_w , function(x){sum(!is.na(x))} )

ml_preds_w_repl <- ml_preds_w
######################################################################
# replace predicted values with corresponding measurements

for ( i in seq_along(numcols_names)){
  ml_preds_w_repl[ ,  paste0(numcols_names[i],"_repl") ] <- ifelse(  ! is.na(ml_preds_w[ , paste0(numcols_names[i],"_meas") ]) ,
                                                                  ml_preds_w[ , paste0(numcols_names[i],"_meas") ],
                                                                  ml_preds_w[ , paste0(numcols_names[i],"_pred") ] )
  ml_preds_w_repl[ ,  paste0(numcols_names[i],"_NULL_med") ] <- ifelse(  ! is.na(ml_preds_w[ , paste0(numcols_names[i],"_meas") ]) ,
                                                                     ml_preds_w[ , paste0(numcols_names[i],"_meas") ],
                                                                     median( ml_preds_w[ , paste0(numcols_names[i],"_meas") ] , na.rm = TRUE ) )
}

head(ml_preds_w_repl)
######################################################################
# sanity check
head(ml_preds_w_repl[!is.na(ml_preds_w_repl$kappmax_davidi_per_pp_per_s_meas),])

######################################################################
#  select measured vals imputed by ML for export
ml_preds_export <- ml_preds_w_repl[ , grepl("_repl$|react_id|_NULL_med$",colnames(ml_preds_w_repl)) ]

###########
# check for reactions that should be excluded (prep_train_data.r should have taken care of Exch and spont, but double-check)

library(sybil)
load("data/M-models/iML1515/iML1515.rdata")
iML1515_irr <- mod2irrev(iML1515)
#exch reactions? 
grep("EX_" , ml_preds_export$react_id)
grep("acgam" , ml_preds_export$react_id)
#spont reactions? 
iML_spont_reacts <- iML1515_irr@react_id[iML1515_irr@gpr=="s0001"] 
any( ml_preds_export$react_id %in% iML_spont_reacts )
# ATPM?
grep("ATPM" , iML1515_irr@react_id , value = TRUE)
grep("ATPM" , ml_preds_export$react_id , value = TRUE)
# remove ATPM
ml_preds_export <- ml_preds_export[ ml_preds_export$react_id != "ATPM" , ]

# biomass
grep("BIOMASS" , iML1515_irr@react_id , value = TRUE)
grep("BIOMASS" , ml_preds_export$react_id , value = TRUE)
# remove biomass reaction:
ml_preds_export <- ml_preds_export[ !grepl("BIOMASS" , ml_preds_export$react_id ) , ]

######################################################################
# additional checks, exports should have no NAs and be larger 0
stopifnot(! any( is.na(ml_preds_export) ))
stopifnot( all( ml_preds_export>0 ) )

######################################################################
# correlation among exported vectors
cor(log10(ml_preds_export[,2:5]))
######################################################################
# export wide format 

head(ml_preds_export)

write.csv(ml_preds_export , file= "output_data/kappmax/KO_MFA_Davidi_kappmax_w_ML/KO_MFA_Davidi_kappmax_w_ML_4x_cutoff.csv" , row.names = FALSE)
