#040219
# read PSM_files and combine them into one df 
# rename some sample ids 
# apply Quality filters

# save to rdata

ups2 <- read.csv("data/proteomics/UPS2.csv",stringsAsFactors = FALSE)
ups2$UniProt_Accession_ext <- paste0( ups2$UniProt_Accession , "ups" )



# vector of PSM files provided by Anaamika
files <- list.files("data/proteomics/ALE_Proteomics_1_26_2019/")

# split at _ to get sample ids
sample_ids <- sapply( strsplit(x = files , split = "_") , "[" , 1 )

files <- paste0( "data/proteomics/ALE_Proteomics_1_26_2019/" , files )

names(files) <- sample_ids

# define column classes manually
col_classes <- c( "character" ,"character","character","character","character","character",NA,NA,
    "character","character",NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,"character","character",NA,NA,
    NA,NA,"character","character",NA,NA,NA,NA,
    NA )

assign("last.warning", NULL, envir = baseenv())

options( warn = 1 )

for( i in seq_along(files)){
  print(paste(i,names(files)[i]))
  
  file_i_df <- read.delim( files[i] , stringsAsFactors = FALSE , quote = "\""  ) 
  
  # add column for identifying sample
  file_name_i <- names(files)[i]
  file_i_df[ "sample_id"  ] <- gsub("(.+)B1|B2$","\\1",file_name_i)
  file_i_df[ "bio_rep"  ] <- gsub(".+(B1|B2)$","\\1",file_name_i)
  
  if(i==1){
    full_df <- file_i_df
  }else{
    full_df <- rbind(full_df , file_i_df)  
  }
}

head(full_df)
lapply(full_df, class)

# re-name some samples
table(full_df$sample_id,full_df$bio_rep)

full_df$sample_id[ full_df$sample_id=="57" ] <- "C13_A3"
full_df$sample_id[ full_df$sample_id=="133" ] <- "C13_A1"
full_df$sample_id[ full_df$sample_id=="ALE4" ] <- "GLU_A4"
full_df$sample_id[ full_df$sample_id=="F380" ] <- "GLU_A8"
full_df$sample_id[ full_df$sample_id=="F406" ] <- "GLU_A6"
full_df$sample_id[ full_df$sample_id=="F418" ] <- "GLU_A10"

table(full_df$sample_id,full_df$bio_rep)


########################################################
# Anaamikas QC
full_df_s <- full_df[ full_df$Confidence == "High",  ]
full_df_s <- full_df_s[ full_df_s$PSM.Ambiguity != "Rejected",  ]
# use Protein.Accessions to replace missing Master protein accessions:
full_df_s$Master.Protein.Accessions[ full_df_s$Master.Protein.Accessions == "" ] <- full_df_s$Protein.Accessions[ full_df_s$Master.Protein.Accessions == "" ] 
# remove peptides that match multiple master proteins
full_df_s <- full_df_s[ ! grepl( ";" , full_df_s$Master.Protein.Accessions )  ,  ]

########################################################

# ambiguous masters after QC: 
length( grep(";",full_df_s$Master.Protein.Accessions,value = TRUE) )
sum( full_df_s$Master.Protein.Accessions == "" )

unique(grep("ups" , full_df_s$Master.Protein.Accessions, value = T ))

# check covergae of the standard
cbind(ups2 , ups2$UniProt_Accession_ext %in% full_df_s$Master.Protein.Accessions )
# note that only half of the eight 0.5 fmol standard proteins were deteccted in at elast one sample, 
# and only five of the 5fmol proteins

# looking at the ups coverage in a single sample:
cbind(ups2 , abs = ups2$UniProt_Accession_ext %in% full_df_s$Master.Protein.Accessions[full_df_s$sample_id =="pgiEVO8"] )
# lowest abundance not found, only one protein for the 2nd lowest

#######################################################

save( full_df_s , file = "data/proteomics/ALE_Proteomics_1_26_2019.rdata" , compress = TRUE)
