
write_outputs <- TRUE

library(tidyverse)

source("analyses/auxilliary_scripts/MOMENT/model_proteome_MOMENT.r"  )

#############################
# load and prepare MOMENT data

# created by "MOMENT_w_kappmax_comp_Schmidt.R"
load("analyses/kapp_analysis/MOMENT_validation/MOMENT_saved_preds/Moment_Schmidt_data_exclude_transp_TRUE_4x_cutoff.rdata")

mom_ab <- MOMENT_schmidt_out$out_det

# note that the schmidt data is in fg/cell, while MOMENT is in weight/weightDw
mom_ab$mech_model<- "MOMENT"
mom_ab$mech_model[mom_ab$source=="Schmidt_data"]<- "experiment"


#############################
# load and prepare ME data 

# load long format abudnace data 
me_ab <- read.csv( "data/ME_model_predictions/mes_sim_abus_KOALE_all_conditions_4x_cutoff_040619.csv" , stringsAsFactors = FALSE )

names(me_ab)
# get bnumbers that are associated with transport and Membrane in ME 
me_tranp_Mem_bnums <- unique( me_ab$gene_bnum[ me_ab$annotation %in% c("Membrane","Transport" ) ] )
length(me_tranp_Mem_bnums)
me_ab <- me_ab[ ! me_ab$gene_bnum %in% me_tranp_Mem_bnums ,]


# remove index col
me_ab <- me_ab[ , names(me_ab)!="X"]
# select weight and rename to match MOMENT
me_ab <- me_ab %>% 
  select(- c(abundance..mmol. , sum.metabolic.reactions , annotation) ) %>% 
  rename( abundance = abundance.mg.) %>% 
  filter( growth_condition != "Glycerol + AA") %>% # conditioons not run for MOMENT
  select( source ,growth_condition ,gene_bnum, abundance) 
me_ab$mech_model <- "ME"


##########################################
# combine MOMENT and ME 

 ab <- rbind( mom_ab , me_ab )

##################################


# create character of genes that are found in every condition and source

list_of_bnums_MOM <- split( ab$gene_bnum[ab$mech_model=="MOMENT"] ,  paste0( ab$source[ab$mech_model=="MOMENT"] , "_" , ab$growth_condition[ab$mech_model=="MOMENT"] ) , drop = TRUE )
list_of_bnums <- split( ab$gene_bnum ,  paste0( ab$source , "_" , ab$growth_condition) , drop = TRUE )

bnum_overlap <- Reduce( intersect , list_of_bnums )
length(bnum_overlap)

#####################################
# compute performance on intersected 

intersected_perf <- intersect_mechMod_preds(ab , data_source_name = "Schmidt_data")

intersected_perf$out_df

intersected_perf$out_df %>% filter(mech_model == "ME" & ! is.na(RMSE) )


if(write_outputs){
  save(intersected_perf, file = paste0("analyses/kapp_analysis/MOMENT_validation/MOMENT_saved_preds/Moment_Schmidt_data_exclude_transp_",exclude_transp,"_intersect_4x_cutoff.rdata") )
  write.csv(intersected_perf$out_df , file = paste0("analyses/kapp_analysis/MOMENT_validation/MOMENT_saved_preds/Moment_Schmidt_perfor_stats_exclude_transp_",exclude_transp,"_intersect_4x_cutoff.csv") )
}


