library(dplyr)

# load mutation annotation
mut_annot <- read.csv("data/ALE_db_genome_mutation_annotation/CCK_genetic_targets.csv", stringsAsFactors = FALSE)
mut_annot <- mut_annot %>% group_by( exp ) %>% 
  filter(flask != 0) %>% # exclude non-evolved strains
  mutate(rank = dense_rank(ale)) %>% as.data.frame()
mut_annot$shared_id <- with(mut_annot, paste0(exp,"_",rank,"_I",isolate) )

head(mut_annot[mut_annot$exp == "CCK_pgi",c(1:5,31)])

# sanity check, compare pgi1 clone 2 to the ALEdb entry:
mut_annot %>% filter( shared_id == "CCK_pgi_1_I2")

#################################
# filter relevant ALEs:

# load list of the 21 strains used in the study
samplenames <- read.csv("data/ALE_proteomics_sample_names.csv", stringsAsFactors = FALSE )
head(samplenames)
# convert to match notation
samplenames$isolate <- sub( ".+EP I(\\d*) .+" , "\\1" , samplenames$Strain.Alias )
samplenames$exp_suffix <- tolower( substr( samplenames$Strain.Alias, 5,7 ) )
samplenames$exp <- paste0("CCK_",samplenames$exp_suffix)
samplenames$ale <- sub( ".+Evo0(\\d*)EP .+" , "\\1" , samplenames$Strain.Alias )
samplenames$shared_id <- with(samplenames,paste0(exp,"_",ale,"_I",isolate))

# inner join 
ko_ALE_muts<- merge(samplenames, mut_annot , by = "shared_id")

# make sure we catch all 21 samples
length(table(ko_ALE_muts$shared_id))

table(ko_ALE_muts$shared_id)

# remove all muts that are found in all strains:
all_strain_mut_ids <- names(table(ko_ALE_muts$Mut.ID)[ table(ko_ALE_muts$Mut.ID) == 21 ])
ko_ALE_muts[ko_ALE_muts$Mut.ID %in% all_strain_mut_ids , c(1:8,15:25)] %>% arrange(Mut.ID)

ko_ALE_muts <- ko_ALE_muts[ ! ko_ALE_muts$Mut.ID %in% all_strain_mut_ids , ]


################################
library(stringr)
# extract bnumbers
ko_ALE_muts$bnum1 <- sub(".+bnum\': \'(b\\d{4}).+" , "\\1" , ko_ALE_muts$genetic.features)
ko_ALE_muts$bnum2 <- sub(".+bnum\': \'b\\d{4}/(b\\d{4}).+" , "\\1" , ko_ALE_muts$genetic.features)

bnum_matches <- str_match_all( ko_ALE_muts$genetic.features,"(b\\d{4})" )

# first bnumber listed
ko_ALE_muts$bnum1 <- vapply(bnum_matches , function(x){x[1,2]} , FUN.VALUE = "a")
# second bnumber listed
ko_ALE_muts$bnum2 <- vapply(bnum_matches , function(x){ if(nrow(x)>1){
                                                              x[2,2]
                                                        }else{
                                                              as.character(NA)
                                                        }
                                                      } , 
                            FUN.VALUE = c("a"))

# count coding & metabolic 

# count metabolic cases by maching bnumbers to iML 1515:
library(sybil)

load("data/M-models/iML1515/iML1515.rdata")
iML_genes <- iML1515@allGenes # bnumbers

ko_ALE_muts$in_iML1515 <- vapply(bnum_matches , function(x){any(x[,2] %in% iML_genes)} , 
                                   FUN.VALUE = TRUE )

table( coding = ko_ALE_muts$coding , iML = ko_ALE_muts$in_iML1515 )

table( coding = ko_ALE_muts$coding )

ko_ALE_muts %>% filter(in_iML1515 & as.logical(coding) ) %>%
  select(mutation.target.annotation) %>% 
  table() %>% sort(decreasing = TRUE)

################################################

# export 
write.csv(ko_ALE_muts , "output_data/annotated_CCK_mutations_kcatALE_samples.csv",row.names = FALSE)

##################################################
# inspect icd kapps:
t3_1tn_mfa <- read.csv("output_data/kappmax/KO_MFA_kappmax/KO_MFA_kappmax_full.csv" , stringsAsFactors = FALSE )

# average kapps over bio reps
t3_1tn_mfa_m <- t3_1tn_mfa %>% group_by(react_id , sample_id) %>% summarise( kapp_per_s_m = mean(kapp_per_s,na.rm=TRUE) , 
                                                                             Master.Protein.Accessions = first(Master.Protein.Accessions) ,
                                                                             kappmax_per_s_br = first(kappmax_per_s_br),
                                                                             kapp_samp_95CI_lb = first(kapp_samp_95CI_lb),
                                                                             kapp_samp_95CI_ub = first(kapp_samp_95CI_ub)
                                                                             ) %>% 
  ungroup()

# add kapps to mutations 
ko_ALE_muts_w_kapp <- merge(ko_ALE_muts , t3_1tn_mfa_m , by.x = c("sample_id","mutation.target.annotation") ,by.y = c("sample_id","Master.Protein.Accessions") )

ko_ALE_muts_w_kapp %>% select(sample_id, mutation.target.annotation , kappmax_per_s_br , kapp_samp_95CI_lb , kapp_samp_95CI_ub)

############################

t3_1tn_mfa_m %>% filter(Master.Protein.Accessions=="zwf") %>% 
  select(sample_id, Master.Protein.Accessions , kapp_per_s_m , kapp_samp_95CI_lb , kapp_samp_95CI_ub) %>% arrange(sample_id)

ko_ALE_muts_w_kapp %>% filter(mutation.target.annotation=="zwf") %>% select(sample_id, mutation.target.annotation , kappmax_per_s_br , kapp_samp_95CI_lb , kapp_samp_95CI_ub)

# left join mutation info to kapp table
t3_1tn_mfa_m_w_muts <- merge( t3_1tn_mfa_m , ko_ALE_muts,  ,by.x = c("sample_id","Master.Protein.Accessions") , by.y = c("sample_id","mutation.target.annotation"), all.x = TRUE)

#################################################################################
### inspect cases where we have a metabolic enzyme (or complex members) mutated more than once :

######
# icd
icd_kapps <- t3_1tn_mfa_m_w_muts %>% filter(Master.Protein.Accessions=="icd") %>% 
  select(sample_id,  Master.Protein.Accessions , kapp_per_s_m , kapp_samp_95CI_lb , kapp_samp_95CI_ub , Details) %>% arrange(sample_id)
# the one mutated observation has a slightly lower  kapp than the wt, 

#################################
#################################
# export for SI:

ko_ALE_muts_SI <- ko_ALE_muts %>% select(sample_id , Reference.Seq, Position, Mutation.Type , Sequence.Change, 
                                         Details , mutation.target.annotation , coding , range , in_iML1515 ) %>% 
  rename( sample.id = sample_id , in.iML1515 = in_iML1515)

# remove "evo"
ko_ALE_muts_SI$sample.id <- sub("^EVO(1|2)$", "WT\\1" , ko_ALE_muts_SI$sample.id)
ko_ALE_muts_SI$sample.id <- sub("EVO|Evo", "" , ko_ALE_muts_SI$sample.id)
head(ko_ALE_muts_SI)
# sample.id : name of starting strain in ALE
# Reference.Seq: Reference Genome
# Position : Position of mutation
# Mutation.Type : Type of mutation
# Sequence.Change : Sequence change                  
# Details : Details on sequence change
# mutation.target.annotation 
# coding : Did the mutation affect a coding region?              
# range : Positional range of mutation
# in.iML1515 : Are any of the affected genes included in the iML1515 model?

write.csv(ko_ALE_muts_SI , file = "output_data/SM/mutation_table/mutation_table.csv" , row.names = FALSE)
