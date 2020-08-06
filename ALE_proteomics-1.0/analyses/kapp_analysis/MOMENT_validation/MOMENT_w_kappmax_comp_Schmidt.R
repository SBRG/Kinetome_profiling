# run MOMENT and compare to Schmidt data, also export predictions
library(dplyr)

#script parameters
write_outputs <- TRUE
exclude_transp <- TRUE # do not use predictions for transporters 

# parameters and M model for moment
glob.parms <- new.env()
source("analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_global_parms.R" , local = glob.parms )
source("analyses/auxilliary_scripts/MOMENT/model_proteome_MOMENT.r"  )

#######################################
# load df of predictions

kcat.ml.preds <- read.csv(  "output_data/kappmax/KO_MFA_Davidi_kappmax_w_ML/KO_MFA_Davidi_kappmax_w_ML_4x_cutoff.csv" , stringsAsFactors = FALSE)

if(  exclude_transp ){
  transport_react_df_irr <- read.csv("data/M-models/iML1515//iML1515_list_of_transporters.csv" , stringsAsFactors = FALSE)
  transport_react_df_irr$bigg_id <- sub( "_f$" , "" , transport_react_df_irr$bigg_id ) 
  sum(kcat.ml.preds$react_id %in%  transport_react_df_irr$bigg_id)
  kcat.ml.preds[ kcat.ml.preds$react_id %in%  transport_react_df_irr$bigg_id  , sapply(kcat.ml.preds,is.numeric) ] <-  65
}

pred.names <- names(kcat.ml.preds)[-1]

###############################
# load schmidt data
schmidt <- read.csv("data/proteomics/Schmidt_2016/Schmidt_2016_copy_numbers.csv" , stringsAsFactors = FALSE , fileEncoding = "latin1" )
schmidt <- schmidt[,1:34] # delete empty cols
# conditions as they occur in schmidt column names:
conditions <- c( "Acetate" ,"Glycerol","Fumarate", "Succinate" , "Pyruvate", "Glucose" , 
                 "Glucosamine", "Xylose" , "Mannose" , "Galactose" , "Fructose", "LB" )

condition_models <- c( "Acetate" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Acetate_model.R",
                       "Glycerol" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Glycerol_model.R",
                       "Fumarate" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Fumarate_model.R", 
                       "Succinate" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Succinate_model.R" , 
                       "Pyruvate" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Pyruvate_model.R", 
                       "Glucose" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Glucose_model.R",
                       "Glucosamine" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Glucosamine_model.R", 
                       "Xylose"  = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Xylose_model.R", 
                       "Mannose"  = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Mannose_model.R", 
                       "Galactose"  = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Galactose_model.R" , 
                       "Fructose" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_Fructose_model.R" ,
                       "LB" = "analyses/kapp_analysis/MOMENT_validation/Schmidt_models/iML/Schmidt_Moment_LB_model.R" 
)

# convert schmidt copy numbers per cell to fg per cell (n/mol = 6.0221e23 , n /cell * mol/n * g/mol * 1e15 fg/g-> fg/cell)
schmidt[ , conditions ] <- lapply( schmidt[ , conditions ] , function(x) { x * 6.0221e-23 * schmidt$Molecular.weight..Da. * 1e15 } )
schmidt_s <- subset(schmidt , select = c(conditions , "Bnumber" ) )

##################################################################

MOMENT_schmidt_out <- model_proteome_MOMENT( conditions = conditions, # char vector 
                           condition_mods = condition_models, # char of R files to source
                           prot_abs = schmidt_s, # df of protein abundacnes as weight per cell with column "Bnumber" and conditions
                           kcat_df = kcat.ml.preds # wide df of kcat vectors with column bigg.id
   )


##################################################################

MOMENT_schmidt_out$out_df %>% group_by( source ) %>% summarise(mean(R2) , mean(RMSE))
out_df <- MOMENT_schmidt_out$out_df

if(write_outputs){
  save(MOMENT_schmidt_out, exclude_transp, file = paste0("analyses/kapp_analysis/MOMENT_validation/MOMENT_saved_preds/Moment_Schmidt_data_exclude_transp_",exclude_transp,"_4x_cutoff.rdata") )
  write.csv(MOMENT_schmidt_out$out_df , file = paste0("analyses/kapp_analysis/MOMENT_validation/MOMENT_saved_preds/Moment_Schmidt_perfor_stats_exclude_transp_",exclude_transp,"_4x_cutoff.csv") )
}

