# combine KO proteomics with un-lumped MFA maps

#script parameters:
n_samples <- 500 # number of samples when simulating kappmax variability
MC_seed <- 6346
fmol_per_gDW_calib_obs_cutoff <- 5e4 # comparable to 10 polypeptides per fl cytoplasm, as was used in Davidi et al.
mfa_flux_avg_cutoff <- 1e-10 # in mmol gDW-1 hr-1
mfa_flux_factor_ci_cutoff <- 4 # would be 4 in "Observable fluxes were those where  the estimated flux value was at least four times larger than the 95% confidence interval and did not include the value zero.Standard deviations of observable" (10.1021/acs.analchem.5b04914))
remove_outliers <- FALSE 
outliers_to_rm <- NA

library(tidyverse)
options(tibble.width = Inf) # for printing tbls
library(reshape2)

##############################################################
# load and modify data 

#load top3_preds, created in "analyses/data_prep/compute_abundance_proxies.R"
load("data/proteomics/ALE_Proteomics_1_26_2019_top3_preds_wo_ups.rdata")
# subset top3_preds to relevant cols
top3_preds <- top3_preds %>% select("sample_id","bio_rep","Master.Protein.Accessions" , starts_with("fmol_per_gDW_") ) 

# load KO MFA data created in "translate_lumped_mfa_to_sybil_iJO.r":
mfa <- read.csv("output_data/mfa_ids/Table_S16_fluxesDataSampled_sybil_iJO.csv" , stringsAsFactors = FALSE)
# add "mfa" prefix to sampling cols
colnames(mfa)[grep("sampling",colnames(mfa))] <- paste0( "mfa_" , colnames(mfa)[grep("sampling",colnames(mfa))] )

# sanity check bounds
all(mfa$mfa_sampling_lb < mfa$mfa_sampling_ub)
mfa_sample_id_map <- read.csv("data/McCloskey KO strain flux maps/mfa_id_sample_id_map.csv" , stringsAsFactors = FALSE )
mfa <- merge(mfa, mfa_sample_id_map , by = "simulation_id" )


# convert to direction-specific ids and numbers
mfa <- mfa %>% mutate(mfa_sampling_ub_temp = mfa_sampling_ub,
                        mfa_sampling_lb_temp = mfa_sampling_lb,
                        rxn_id_sybil_rev = rxn_id_sybil,
                        rxn_id_sybil = if_else( mfa_sampling_ave < 0 , paste0(rxn_id_sybil,"_b") , rxn_id_sybil  ) ,
                        rxn_id = if_else( mfa_sampling_ave < 0 , paste0(rxn_id,"_b") , rxn_id  ) ,
                        mfa_sampling_ub = if_else( mfa_sampling_ave < 0 , - mfa_sampling_lb_temp , mfa_sampling_ub  ) ,
                        mfa_sampling_lb = if_else( mfa_sampling_ave < 0 , - mfa_sampling_ub_temp , mfa_sampling_lb  ) ,
                        mfa_sampling_ave = if_else( mfa_sampling_ave < 0 , -mfa_sampling_ave, mfa_sampling_ave) ) %>%
                        select(-c(mfa_sampling_ub_temp,mfa_sampling_lb_temp))

# check 
stopifnot(all( mfa$mfa_sampling_lb <= mfa$mfa_sampling_ave & mfa$mfa_sampling_ave <= mfa$mfa_sampling_ub ))

# mark observable fluxes (as defined in Mccloskey 2016 (10.1021/acs.analchem.5b04914) , but where the factor 4  is set by mfa_flux_factor_ci_cutoff)
mfa <- mfa %>% mutate(flux_observable =  mfa_flux_factor_ci_cutoff * abs(mfa_sampling_ub - mfa_sampling_lb) < abs(mfa_sampling_ave) &  
                        ( ( mfa_sampling_ub>0 & mfa_sampling_lb>0) | (mfa_sampling_ub<0 & mfa_sampling_lb<0) ) &
                        mfa_sampling_ave > mfa_flux_avg_cutoff
)

#g-R maps for iJO1366
# note that these only contain reactions catalyzed by homomers, a fact that is used in the merge
# also note that even though all reactions are catalyzed by homomers, theres are reactions where the catalyzing gene is catalyzing as a heteromer as well
g_R_one_to_many <- read.csv( "data/M-models/iJO1366/iJO1366_gene-reaction-map_one-to-many_irr.csv", stringsAsFactors = FALSE)
g_R_one_to_one <- read.csv( "data/M-models/iJO1366/iJO1366_gene-reaction-map_one-to-one_irr.csv", stringsAsFactors = FALSE)

###############################################################

# inner join of t3 abs with g R map (removing non-matched prot abs)
#t3_1t1 <- merge( top3_preds , subset( g_R_one_to_one ,select =  - c(X,organism) ) , by.x = "Master.Protein.Accessions", by.y = "name" )
t3_1tn <- merge( top3_preds ,  
                 subset( g_R_one_to_many ,select =  - c(X,organism) ) , 
                 by.x = "Master.Protein.Accessions", by.y = "name"  )

# count mapped cases (note that n reaction will be smaller for the 1tn map)
t3_1tn %>% filter(one_to_one) %>% group_by(sample_id,bio_rep) %>% summarise(n = n())
t3_1tn %>% group_by(sample_id,bio_rep) %>% summarise(n = n())

length(intersect(t3_1tn$react_id , mfa$rxn_id_sybil))
sort(intersect( t3_1tn$sample_id , mfa$sample_id ) )
# number of reactions we have protein for in some condition, 
# but the reaction does not appear in the MFA map
length(setdiff( t3_1tn$react_id , mfa$rxn_id_sybil ) )

# add mfa data
# left join to keep all proteins, even if no mfa data available 
t3_1tn_mfa <- merge( t3_1tn , mfa , by.x =c( "react_id" , "sample_id" ) , 
                     by.y = c( "rxn_id_sybil" , "sample_id" )  , all.x = TRUE )
# remove non-KO samples that we dont have MFA data for that are left from the left join
t3_1tn_mfa <- t3_1tn_mfa %>% filter(sample_id %in% mfa$sample_id)

# confirm that the mfa map never contains flux in more than one direction (use this fact when deciding if all fluxes per gene were measured)
mfa %>% mutate(rxn_id_sybil_rev = sub("(.+)_b$","\\1",rxn_id_sybil)) %>%  group_by(simulation_id , rxn_id_sybil_rev) %>% 
  summarise(n = n()) %>% ungroup() %>%  summarise(all(n==1))

# check if all iJO reactions associated with gene i are occuring in the MFA map.
# we do this to decide whether to calc kapp in 1tn cases
t3_1tn_mfa <- t3_1tn_mfa %>% group_by( Master.Protein.Accessions , sample_id , bio_rep ) %>% 
  mutate( all_flx_per_gene_meas_and_homo_r = ifelse(all(!reversible) , # logical indicating whether all fluxes per gene were measured, and if all catalyzed reactions are catalyzed by a homomeric enzyme
                                                    !any(is.na(mfa_sampling_ave)) & all_r_per_g_homo  ,
                                                    sum(is.na(mfa_sampling_ave))==sum(direction=="bwd") & all_r_per_g_homo  # mfa map never contains flux in more than one direction 
                                                    ) , 
          n_flux_per_gene_meas = sum(!is.na(mfa_sampling_ave)) , 
          mfa_sampling_ave_sum_per_gene = if_else( all_flx_per_gene_meas_and_homo_r , # sum of flux of all reactions measured for gene 
                                                   sum(mfa_sampling_ave , na.rm = TRUE), # removal of NAs is still necessary in the case of irreversible reactions, where NAs are found for one direction despite all_flx_per_gene_meas_and_homo_r can be TRUE
                                                   as.numeric(NA) )
          )

head(as.data.frame(t3_1tn_mfa[t3_1tn_mfa$all_flx_per_gene_meas_and_homo_r , ]),50)
as.data.frame(t3_1tn_mfa[t3_1tn_mfa$all_flx_per_gene_meas_and_homo_r & t3_1tn_mfa$reversible , ]) %>% 
  filter(Master.Protein.Accessions =="folD" & sample_id == "ptsEvo2" & bio_rep == "B1") %>%  arrange(Master.Protein.Accessions, sample_id) #%>% slice(1000:1050)
                
head(as.data.frame(t3_1tn_mfa[ t3_1tn_mfa$reversible , ] %>% arrange(Master.Protein.Accessions, sample_id)),50)

head(as.data.frame(t3_1tn_mfa),20)

# inspect reversibel reacts 
rev_reacts <- unique(grep("_b$",t3_1tn_mfa$react_id,value= TRUE))
rev_reacts <- sort( c(rev_reacts , sub("(.+)_b$" ,"\\1",rev_reacts)) )
t3_1tn_mfa %>% filter( react_id %in% rev_reacts & sample_id =="EVO1")

head(as.data.frame(t3_1tn_mfa[t3_1tn_mfa$Master.Protein.Accessions == "mtn" & t3_1tn_mfa$sample_id=="EVO1",]),20)

#####
# some counts 
# first make sure tht we have one observation for each sample, biorep , reaction combination

# count 1tn reactions for which we have all fluxes measured
t3_1tn_mfa %>% filter( (! one_to_one) & fmol_per_gDW_calib_avg > fmol_per_gDW_calib_obs_cutoff & flux_observable ) %>% group_by(sample_id , react_id) %>% 
  summarise(all_flx_per_gene_meas_and_homo_r = first(all_flx_per_gene_meas_and_homo_r) ) %>%
    group_by(all_flx_per_gene_meas_and_homo_r, sample_id) %>% summarise( n_reacts = n() ) %>% group_by(all_flx_per_gene_meas_and_homo_r) %>%
  summarise(avg_n_reacts = mean(n_reacts) ) # average number 1tn reactions w and wo all flx measured
# => the majority of 1tn reactions does not have all fluxes measured

# count observed 1-t-1 pairs
# (group bio replicates)
t3_1tn_mfa %>%   group_by( sample_id ,react_id) %>% 
  summarise( onetone_meas = any( one_to_one & fmol_per_gDW_calib_avg > fmol_per_gDW_calib_obs_cutoff & flux_observable) ) %>% 
  group_by( sample_id ) %>% summarise(n_reacts = sum(onetone_meas , na.rm = TRUE) ) %>% 
  summarise(avg_n_reacts = mean(n_reacts) ) # average number 1t1 reactions with flux and protein measured


########################################################

# compute kapps
t3_1tn_mfa <- t3_1tn_mfa %>% ungroup() %>% 
  mutate( kapp_per_s = if_else( fmol_per_gDW_calib_avg > fmol_per_gDW_calib_obs_cutoff & flux_observable , #& one_to_one, # case: 1t1
                              ( mfa_sampling_ave / 3600 )  / ( fmol_per_gDW_calib_avg * 1e-12 ) , 
                              as.numeric(NA)

                      )
          )

t3_1tn_mfa <- t3_1tn_mfa %>% group_by( react_id , bio_rep ) %>% # kappmax per bio rep
  mutate( kappmax_per_s_br = ifelse( ! all( is.na(kapp_per_s) ) , max(kapp_per_s , na.rm = TRUE) , as.numeric(NA) ) ,
          kappmax_per_s_br_n = sum(!is.na(kapp_per_s) ),
          kappmax_per_s_br_sample_id_of_max = ifelse( ! all(is.na(kapp_per_s) ) , sample_id[ which.max(kapp_per_s)] , NA )
          ) %>% 
  group_by(react_id) %>%
  # note the hard-coded bio rep grouping
  mutate( kappmax_per_s_br_avg = mean( c( first(kappmax_per_s_br[bio_rep == "B1"]),first(kappmax_per_s_br[bio_rep == "B2"]) ) , na.rm = TRUE),
          kappmax_per_s_br_sd = sd(c( first(kappmax_per_s_br[bio_rep == "B1"]),first(kappmax_per_s_br[bio_rep == "B2"]) ) ,  na.rm = TRUE),
          kappmax_per_s_br_avg_n = sum(!is.na(c( first(kappmax_per_s_br[bio_rep == "B1"]),first(kappmax_per_s_br[bio_rep == "B2"]) )))
          )
      
  
#######################################################
#######################################################
# MC simulation of kappmax variability (sample for 1t1)

# First, summarise biological replicates of protein abundance and estimate standard deviation
# select relevant columns for this 
t3_1tn_mfa_s <- t3_1tn_mfa %>% 
  select(react_id , sample_id , Master.Protein.Accessions , bio_rep ,
         fmol_per_gDW_calib_avg , mfa_sampling_ave , mfa_sampling_ave_sum_per_gene ,mfa_sampling_var , one_to_one , flux_observable , all_flx_per_gene_meas_and_homo_r ) %>% 
  group_by(sample_id , Master.Protein.Accessions , react_id ) %>% 
  summarise( #react_id = first(react_id),
             fmol_per_gDW_calib_tavg_bavg = mean(fmol_per_gDW_calib_avg, na.rm = TRUE) , 
             fmol_per_gDW_calib_tavg_bsd = sd(fmol_per_gDW_calib_avg , na.rm = TRUE) , 
             one_to_one = first(one_to_one) , 
             flux_observable = first(flux_observable),
             mfa_sampling_ave = first(mfa_sampling_ave) , 
             mfa_sampling_ave_sum_per_gene = first(mfa_sampling_ave_sum_per_gene) ,
             mfa_sampling_var = first(mfa_sampling_var),
             all_flx_per_gene_meas_and_homo_r = first(all_flx_per_gene_meas_and_homo_r)
             )


# The actual simulation:

# explore signal to noise ratio in protein abundances (use it to impute missing sds )
hist(log10(t3_1tn_mfa_s$fmol_per_gDW_calib_tavg_bavg / t3_1tn_mfa_s$fmol_per_gDW_calib_tavg_bsd))
plot( log10(t3_1tn_mfa_s$fmol_per_gDW_calib_tavg_bavg) , log10(t3_1tn_mfa_s$fmol_per_gDW_calib_tavg_bsd ) )
# fit lm for signal to noise, use later to impute sds
mu_sd_lm <- lm( log10(fmol_per_gDW_calib_tavg_bsd ) ~ log10(fmol_per_gDW_calib_tavg_bavg) , data = t3_1tn_mfa_s )
summary(mu_sd_lm)

set.seed(MC_seed)

t3_1tn_mfa_MCsamp <- t3_1tn_mfa_s %>% rowwise() %>% 
  # only  observable flux
  filter(  fmol_per_gDW_calib_tavg_bavg > fmol_per_gDW_calib_obs_cutoff & flux_observable & (!is.na(mfa_sampling_ave) | ! is.na(mfa_sampling_ave) ) ) %>% 
  # impute NAs in P abundances with signal/ noise ratio
  mutate( fmol_per_gDW_calib_tavg_bsd_imp = replace( fmol_per_gDW_calib_tavg_bsd , 
                                                 is.na(fmol_per_gDW_calib_tavg_bsd) , 
                                                 10^predict( mu_sd_lm , data.frame(fmol_per_gDW_calib_tavg_bavg  = fmol_per_gDW_calib_tavg_bavg ) ) )) %>% 
  #draw the kapp samples for each gene at each condition, NAs are introduced when no sd for protein is available
  mutate( kapp_samp = #ifelse( one_to_one , # 1t1  case 
                              list( rnorm( n_samples , mean = mfa_sampling_ave / 3600 , sd = sqrt(mfa_sampling_var) / 3600 ) / 
                                    rnorm(n_samples , mean = fmol_per_gDW_calib_tavg_bavg * 1e-12 , sd = fmol_per_gDW_calib_tavg_bsd_imp * 1e-12 )  ),
          kapp_samp_avg = mean(kapp_samp, na.rm = TRUE ),
          kapp_samp_sd = sd(kapp_samp, na.rm = TRUE ),
          kapp_samp_95CI_lb = quantile(kapp_samp,c(0.025) , na.rm = TRUE),
          kapp_samp_95CI_ub = quantile(kapp_samp,c(0.975), na.rm = TRUE ),
          kapp_samp_n = n_samples
          ) %>% 
  unnest() %>%
  # index MC-samples within each reaction within each condition
  group_by( sample_id , react_id ) %>% 
  mutate(MC_samp_index = row_number()) %>% 
  # calc kappmax within a reaction for each MC-sample index (resulting in n_samples unique kappmaxes per reaction)
  group_by(MC_samp_index , react_id) %>% 
  mutate(kappmax_samp = ifelse( ! all(is.na(kapp_samp) ) , max(kapp_samp , na.rm = TRUE), as.numeric(NA) ),
         sample_id_of_max = ifelse( ! all(is.na(kapp_samp) ) , sample_id[ which.max(kapp_samp) ] , as.character(NA) ) # in which condition was this max found?
         ) %>% 
  group_by(react_id , sample_id_of_max) %>% 
  mutate( count_sample_id_of_max = n() )

# check, (works best for small n_samples)
t3_1tn_mfa_MCsamp %>% filter(react_id == "ACGS" | react_id == "ADSS"  ) %>% as.data.frame() %>% head(20)

# summarise simulated kapps
t3_1tn_mfa_MCsamp_summ_kapp <- t3_1tn_mfa_MCsamp %>% group_by(sample_id, Master.Protein.Accessions, react_id) %>% 
  summarise_at( vars(kapp_samp_avg,
                     kapp_samp_sd,
                     kapp_samp_95CI_lb,
                     kapp_samp_95CI_ub,
                     kapp_samp_n),
                first)

# summarise simulated kappmaxes
t3_1tn_mfa_MCsamp_summ_kappmax <- t3_1tn_mfa_MCsamp %>%
  # calc summary stats across all MC samples for each reaction
  group_by(react_id) %>% 
  summarise( kappmax_samp_avg = mean(kappmax_samp, na.rm = TRUE) ,
             kappmax_samp_med = median(kappmax_samp, na.rm = TRUE) ,
             kappmax_samp_sd = sd(kappmax_samp , na.rm = TRUE),
             kappmax_samp_95CI_lb = quantile(kappmax_samp,c(0.025) , na.rm = TRUE),
             kappmax_samp_95CI_ub = quantile(kappmax_samp,c(0.975), na.rm = TRUE ),
             kappmax_samp_maj_sample_id = sample_id_of_max[which.max(count_sample_id_of_max)] # sample most frequently contributing the highest kapp
             )

t3_1tn_mfa_MCsamp_summ_kappmax  %>% as.data.frame() %>% head()

# check object size
format(object.size(t3_1tn_mfa_MCsamp),"Mb")

#######################################################
# add MC sampling results to original kapp,max data

# left join sampling results to main data

t3_1tn_mfa <- merge( t3_1tn_mfa , t3_1tn_mfa_MCsamp_summ_kapp , by=c("sample_id","Master.Protein.Accessions","react_id" ) , all.x = TRUE)
t3_1tn_mfa <- merge( t3_1tn_mfa , t3_1tn_mfa_MCsamp_summ_kappmax , by=c("react_id" ) , all.x = TRUE)
head( t3_1tn_mfa[ ! is.na(t3_1tn_mfa$kappmax_samp_avg) , ] )

########################################################
# summarise to reaction level

t3_1tn_mfa_kappmax <- t3_1tn_mfa %>% group_by(react_id, bio_rep) %>% 
  #summarise biorep-specific stats
  mutate ( kappmax_per_s_n_avg_brs_used  = mean(kappmax_per_s_br_n , na.rm = TRUE) ,
           kappmax_per_s_br_avg_sample_id_of_max  = ifelse (! all( is.na(kappmax_per_s_br)) , kappmax_per_s_br_sample_id_of_max[ which.max(kappmax_per_s_br) ] , NA ) 
           ) %>% 
  group_by(react_id) %>% 
  summarise_at( vars( 
            one_to_one,
            all_r_per_g_homo,
            Master.Protein.Accessions,
            kappmax_per_s_br_avg ,
            kappmax_per_s_br_sd,
            kappmax_per_s_br_avg_n,
            kappmax_per_s_n_avg_brs_used,
            kappmax_per_s_br_avg_sample_id_of_max,
            kappmax_samp_avg,
            kappmax_samp_med,
            kappmax_samp_sd,
            kappmax_samp_95CI_lb,
            kappmax_samp_95CI_ub,
            kappmax_samp_maj_sample_id ) , 
            first
            )

t3_1tn_mfa_kappmax[ ! is.na(t3_1tn_mfa_kappmax$kappmax_per_s_br_avg) , ]


########################################################
# write outputs

write.csv(t3_1tn_mfa_kappmax , file = "output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_summarised.csv")
write.csv(t3_1tn_mfa , file = "output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_full.csv")


# save R objects including parameters
res_li <- list( t3_1tn_mfa_kappmax = t3_1tn_mfa_kappmax , 
                t3_1tn_mfa = t3_1tn_mfa,  
                MC_seed = MC_seed , 
                fmol_per_gDW_calib_obs_cutoff = fmol_per_gDW_calib_obs_cutoff , 
                mfa_flux_avg_cutoff = mfa_flux_avg_cutoff , 
                remove_outliers = remove_outliers , 
                outliers_to_rm = outliers_to_rm)

save( res_li ,
      file = "output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax.rdata")

########################################################

