library(sybil)
library(dplyr)

# check reaction id overalp of Douglas' flux maps with M Model
mfa <- read.csv("data/McCloskey KO strain flux maps/Table_S16_fluxesDataSampled.csv" , stringsAsFactors = FALSE )
# add converted fluxes that were 'unlumped'
mfa_unlumped <- read.csv("output_data/mfa_ids/Table_S16_fluxesDataSampled_unmapped_rxnids_unlumped.csv" , stringsAsFactors = FALSE )
# update some sampling bounds:
err_rows <- mfa_unlumped$sampling_lb > mfa_unlumped$sampling_ub
sum(err_rows)
mfa_unlumped$sampling_lb_org <-mfa_unlumped$sampling_lb
mfa_unlumped$sampling_ub_org <-mfa_unlumped$sampling_ub
mfa_unlumped$sampling_lb[err_rows] <- mfa_unlumped$sampling_ub_org[err_rows]
mfa_unlumped$sampling_ub[err_rows] <- mfa_unlumped$sampling_lb_org[err_rows]
mfa_unlumped <- subset(mfa_unlumped , select = -c(sampling_lb_org,sampling_ub_org) )
sum( mfa_unlumped$sampling_lb > mfa_unlumped$sampling_ub)


head(mfa)
head(mfa_unlumped)

# select shared columns and remove lumped reactions that were converted
mfa_s_f <- mfa %>% select( intersect(colnames(mfa) , colnames(mfa_unlumped) ) ) %>% filter( ! rxn_id %in% mfa_unlumped$lumped_rxn_id )
head(mfa_s_f)
# note taht rbind matches columns by name
mfa_j <- rbind( mfa_s_f , subset(mfa_unlumped , select = - lumped_rxn_id) )

head(mfa_j)
tail(mfa_j)

# convert to sybil reac_ids
mfa_j$rxn_id_sybil <- sub("_LPAREN_e_RPAREN_" , "(e)" , mfa_j$rxn_id )
mfa_j$rxn_id_sybil <- sub("_DASH" , "_" , mfa_j$rxn_id_sybil ) 

#######################################
# manual translation
mfa_j$rxn_id_sybil[mfa_j$rxn_id_sybil == "Ec_biomass_iJO1366_WT_53p95M" ] <- "BIOMASS_Ec_iJO1366_WT_53p95M"
mfa_j$rxn_id_sybil[mfa_j$rxn_id_sybil == "D__LACtex" ] <- "D_LACtex"
mfa_j$rxn_id_sybil[mfa_j$rxn_id_sybil == "EX_glc(e)" ] <- "EX_glc__D(e)"



summary(mfa_j[ grep( "EX_co2\\(e\\)_unlabeled" , mfa_j$rxn_id_sybil ) ,  ]$sampling_ave )
summary(mfa_j[ grep( "EX_glc\\(e\\)" , mfa_j$rxn_id_sybil ) ,  ]$sampling_ave )
summary(mfa_j$sampling_ave)
mfa_j[ order(mfa_j$sampling_ave, decreasing = T) ,  ] 

######
# also map to strain IDs 

write.csv(mfa_j , "output_data/mfa_ids/Table_S16_fluxesDataSampled_sybil_iJO.csv" , row.names = FALSE)
