########
# WARNING run on fresh session, plyr/ dplyr conflicts can occur

####
# parameters
# use top 75 percent most abundant prots
#abundance_cutoff <- 0.25 # percentile/100 for under which smaples are discarded
abundance_cutoff <- 0 # percentile/100 for under which smaples are discarded

library(dplyr)
library(ggplot2)
library(reshape2)


# created in "compute_abundance_proxies.R"
load("data/proteomics/ALE_Proteomics_1_26_2019_ab_proxies.rdata")

#########################################################
# get iML genes to look at ametabolic genes

# gene name map necessary to translate model bnumbers to gene names 
gene_name_map <- read.csv("data/data_maps/ecoli_K12_gene_name_bnumber_map_uniprot_270219.csv" , sep = "\t" , header = TRUE , stringsAsFactors = FALSE)
# resolve ambiguity by using first entry (mostly a problem for the gene names), in the description "In this field, we provide the "primary" gene name (the first one used in 
# the GN line of the relevant Swiss-Prot entry) and alternative gene name(s)."
gene_name_map$first_bnum <- sapply( strsplit( gene_name_map$Ordered_locus_name , split = ";" ) , function(x){x[1]} )
gene_name_map$first_gene_name <- sapply( strsplit( gene_name_map$Gene_name , split = ";" ) , function(x){x[1]} )


library(sybil)

load("data/M-models/iML1515/iML1515.rdata")
iML_genes <- iML1515@allGenes # bnumbers
head(iML_genes)
# two options: 
iML_genes_names_first <- gene_name_map$first_gene_name[ gene_name_map$first_bnum %in% iML_genes]
iML_genes_names_w_syns <- gene_name_map$Gene_name[ gene_name_map$first_bnum %in% iML_genes]
iML_genes_names_w_syns <- unique( unlist(strsplit( iML_genes_names_w_syns , ";")) )
length(iML_genes_names_w_syns);length(iML_genes)

##########################################################
# sample correlation in top3

# use top 75 percent most abundant prots ()
top3_quant_cutoff <- quantile( ab_proxies$top3 , abundance_cutoff )
log10(top3_quant_cutoff)

# remove UPS for this part of the analysis
top3 <- ab_proxies %>% #select( - ibaq ) %>% 
  filter( !grepl("ups$", Master.Protein.Accessions ) )

# filter low-abundance proteins
top3_f <- ab_proxies %>% filter( top3 > top3_quant_cutoff ) %>% #select( - ibaq ) %>% 
  filter( !grepl("ups$", Master.Protein.Accessions ) )
#filter non-metabolic proteins
#top3_f <- ab_proxies %>% filter( Master.Protein.Accessions %in% iML_genes_names_w_syns ) %>% select( - ibaq ) %>% filter( !grepl("ups$", Master.Protein.Accessions ) )
top3_f_repmean <- top3_f %>% group_by(sample_id , Master.Protein.Accessions) %>% summarise(top3_repmean = mean(top3,na.rm = TRUE)) 

top3_wide <- acast( top3 , formula =  Master.Protein.Accessions ~ sample_id + bio_rep , value.var = "top3")
top3_f_wide <- acast( top3_f , formula =  Master.Protein.Accessions ~ sample_id + bio_rep , value.var = "top3")
top3_f_repmean_wide <- acast( top3_f_repmean , formula =  Master.Protein.Accessions ~ sample_id , value.var = "top3_repmean")
#top3_wide <- acast( top3 , formula =  Master.Protein.Accessions ~ sample_id , value.var = "top3")

# only use complete obs
top3_f_wide_compl_obs <- top3_f_wide[ apply(top3_f_wide,1,function(x){!any(is.na(x))}), ]

# save set of genes observed in all samples
fully_cov_genes <- rownames( top3_f_wide )[ apply( top3_f_wide , 1 , function(x) !any(is.na(x)) ) ]
fully_cov_genes_sel <- rownames( top3_f_wide )[ apply( top3_f_wide , 1 , function(x) !any(is.na(x)) ) ]

########################################################
# correlation between bio reps
top3_f_selfj <- merge( top3_f ,top3_f , by = c( "sample_id" , "Master.Protein.Accessions") )
top3_f_selfj <- top3_f_selfj[ top3_f_selfj$bio_rep.x == "B1" & top3_f_selfj$bio_rep.y == "B2" , ]
replicate_R2s<- top3_f_selfj %>% group_by(sample_id) %>% 
  filter( ! is.na(top3.x) & ! is.na(top3.y) ) %>% 
  summarise( rep_R2 = cor( log10(top3.x+1) , log10(top3.y+1)  , use = "pairwise.complete")^2 ,
             n_for_R2 = n()
             ) %>% as.data.frame()

#a <- top3_f_selfj %>% filter( Master.Protein.Accessions %in% fully_cov_genes_sel )
replicate_R2s_full_obs <- top3_f_selfj %>% filter( Master.Protein.Accessions %in% fully_cov_genes_sel ) %>% 
  group_by(sample_id) %>% summarise( rep_R2 = cor( log10(top3.x+1) , log10(top3.y+1)  , use = "pairwise.complete")^2 ) %>% as.data.frame()
replicate_R2s_full_obs %>% filter(rep_R2<0.85)

summary(replicate_R2s$rep_R2)
replicate_R2s %>% ggplot(aes(x=rep_R2)) + geom_histogram() + xlim(NA,1) + xlab(expression("R"^2*" between biological replicates (top3)") )

#######################################################
# coverage:

# number of genes detected in at least one of the two bio reps
samp_sizes <- ab_proxies %>% 
  filter(top3>0) %>% 
  group_by(sample_id ,Master.Protein.Accessions) %>% 
  summarise(f = first(Master.Protein.Accessions)) %>% 
  group_by(sample_id) %>% 
  summarise( n_genes_observed = n() )

########################################################
# merge coverage and rep correlations

out_dat  <- merge(replicate_R2s , samp_sizes , by = "sample_id")

########################################################
# translate ids and filter to samples used in publication:

source("data/proteomics/sample_id_name_map.R")

# only keep samples used in the kapp study
out_dat <- out_dat %>%  filter(sample_id %in% names(sample_id_map) )
out_dat$sample_name <- sample_id_map[match(out_dat$sample_id , names(sample_id_map))]
out_dat <- out_dat %>% select(sample_name , rep_R2 , n_genes_observed)
names(out_dat) <- c("strain" , "R2 between biological replicates" , "number of proteins detected")

###########################################################
# output data:

write.csv( out_dat, file = "output_data/SI/proteomics_rep_R2_and_coverage/proteomics_rep_R2_and_coverage.csv" , row.names = FALSE )

summary(out_dat$`R2 between biological replicates`)
summary(out_dat$`number of proteins detected`)



