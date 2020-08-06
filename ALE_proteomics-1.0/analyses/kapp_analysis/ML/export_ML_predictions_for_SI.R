library(dplyr)

ml_preds_export <- read.csv( file= "output_data/kappmax/KO_MFA_Davidi_kappmax_w_ML/KO_MFA_Davidi_kappmax_w_ML_4x_cutoff.csv" , stringsAsFactors = FALSE)

# rename variables
names(ml_preds_export) <- sub( "_repl" , "_ensemble_model" , names(ml_preds_export))
names(ml_preds_export) <- sub( "_NULL_med" , "_median_imputed" , names(ml_preds_export))

ml_preds_export <- ml_preds_export %>%  rename( "kcat_in_vitro_per_s_ensemble_model" = "kcat_iv_ML_per_AS_per_s_ensemble_model",
                                                "kcat_in_vitro_per_s_median_imputed" = "kcat_iv_ML_per_AS_per_s_median_imputed"
                                                )

# exclude transporters adn outer memebrane porins
transport_react_df_irr <- read.csv("data/M-models/iML1515/iML1515_list_of_transporters.csv" , stringsAsFactors = FALSE)
transport_react_df_irr$bigg_id <-  sub("_f$" , "" , transport_react_df_irr$bigg_id)
sum(ml_preds_export$react_id %in%  transport_react_df_irr$bigg_id)

setdiff( transport_react_df_irr$bigg_id , ml_preds_export$react_id)

ml_preds_export_wo_transp <- ml_preds_export[ ! ml_preds_export$react_id %in%  transport_react_df_irr$bigg_id  , ] 


# filter exchanges
exchnge_rxns_log <- grepl( "ex$|ex_f$|ex_b$|ATPM|BIOMASS" , ml_preds_export_wo_transp$react_id) 
sum(exchnge_rxns_log)
ml_preds_export_wo_transp <- ml_preds_export_wo_transp[ !exchnge_rxns_log , ]


summary(ml_preds_export_wo_transp)


# filter spontaneous reactions 
library(sybil)
load("data/M-models/iML1515/iML1515.rdata")
spont_log <- grepl("s0001" , iML1515@gpr )
sum(spont_log)
spont_sybil_ids <- iML1515@react_id[ spont_log ]
spont_df_log <- gsub( "_f$|_b$" , "", ml_preds_export_wo_transp$react_id ) %in% spont_sybil_ids
sum(spont_df_log)
ml_preds_export_wo_transp_wo_spont <-  ml_preds_export_wo_transp[ ! spont_df_log ,  ]

head(ml_preds_export_wo_transp_wo_spont)

# export for SI
write.csv( ml_preds_export_wo_transp_wo_spont , file = "output_data/SM/ML_predictions/turnover_number_data_w_ML_predictions.csv" , row.names = FALSE )
