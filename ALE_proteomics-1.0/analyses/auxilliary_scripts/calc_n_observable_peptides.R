library("seqinr")

ind_f001_faa <- read.fasta( file = "data/ALEdb_genomes/variant_files/IND_MG1655_F001/IND_MG1655_F001_translated.fasta" , 
                            seqtype = "AA" , as.string = TRUE ) 

concat_faa <- c(ind_f001_faa)

head(ind_f001_faa)
tail(ind_f001_faa)

ups_log <- (grepl( "HUMAN_UPS", names(concat_faa) ) )

names(concat_faa)[ ups_log ] <- gsub( "(.+ups)\\|.+HUMAN_UPS" , "\\1" , names(concat_faa)[ ups_log ] )

out_df <- data.frame( prot_id = names(concat_faa) , n_observable_peptides = NA , n_Aminoacids = NA )

for( i in 1:nrow(out_df) ){
    entry_i <- concat_faa[[ out_df$prot_id[i] ]] 
    # add ~ at trypsin cuts
    entry_i <- gsub("(K|R)","\\1~",entry_i)
    # split at ~ 
    peps <- strsplit( entry_i , split = "~" )[[1]]
    # count peps longer 5 as (6 as is the shortest ones in Anaamikas first 12 samples) 
    out_df$n_observable_peptides[i] <- length( nchar(peps) >=6 )
    out_df$n_Aminoacids[i] <- nchar(entry_i)
}

out_df[ duplicated(out_df$prot_id) ,  ]

write.csv( out_df , row.names = FALSE,  file = "data/ALEdb_genomes/variant_files/IND_MG1655_F001/n_observable_peptides.csv" )
