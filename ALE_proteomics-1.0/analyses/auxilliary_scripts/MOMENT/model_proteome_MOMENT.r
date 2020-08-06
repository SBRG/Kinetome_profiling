
source("analyses/auxilliary_scripts/MOMENT/MOMENT_Main_functions_1-1.R")

model_proteome_MOMENT <- function( conditions, # char vector 
                                   condition_mods, # char of R files to source
                                   prot_abs, # df of protein abundacnes as weight per cell with column "Bnumber" and conditions
                                   data_source_name = "Schmidt_data", # name of the data used in outputs
                                   kcat_df # wide df of kcat vectors with column react_id
      ){
      
      stopifnot( all(conditions %in% names(condition_mods) ))
      stopifnot( "Bnumber" %in% colnames(prot_abs) )
      stopifnot( "react_id" %in% colnames(kcat_df) )
  
      ######
      # initalize output df list
      out_li <- vector( mode = "list" , length = length(conditions) )
      names(out_li) <- conditions
      out_li <- lapply(out_li , function(x)  data.frame( "source" = pred.names , "R2" = NA , "RMSE" = NA , "n" = NA , "obj" = NA ) )
      for(i in seq_along( out_li ) ){
        out_li[[i]]$condition <- names(out_li)[i]
      }
      
      # record genes used in model-data comparison
      genes_used_li <- vector( mode = "list" , length = length(conditions) )
      names(genes_used_li) <- conditions
      for(i in seq_along( genes_used_li ) ){
        genes_used_li[[i]] <- vector( mode = "list" , length = length(pred.names) )
        names(genes_used_li[[i]]) <- pred.names
      }
      
      # record fluxes model-data comparison
      fluxes_li <- vector( mode = "list" , length = length(conditions) )
      names(fluxes_li) <- conditions
      for(i in seq_along( fluxes_li ) ){
        fluxes_li[[i]] <- vector( mode = "list" , length = length(pred.names) )
        names(fluxes_li[[i]]) <- pred.names
      }
      
      # out details
      out_det <-data.frame( source = character() , 
                            growth_condition = character() , 
                            gene_bnum = character() , 
                            abundance = numeric(), # note that this will hav mixed units #proetin weight per dw (weight fraction) for MOMENT pred nad fg / cell for data
                            stringsAsFactors = FALSE
      )

      for( k in seq_along(conditions) ){
        condition_k <- conditions[k]
        # load condition-specific M model 
        cond.model <- new.env()
        source( condition_models[ condition_k ] ,
                local = cond.model )
        
        prot_abs_k <- subset(prot_abs , select = c("Bnumber"  , condition_k ) ) 
        
        out_det <- rbind( out_det , 
                          data.frame( source = data_source_name , 
                                      growth_condition = condition_k ,
                                      gene_bnum = prot_abs_k$Bnumber ,
                                      abundance = prot_abs_k[,condition_k] , #fg /cell
                                      stringsAsFactors = FALSE
                          )
        )
        
        for (i in seq_along(pred.names) ){
          kcat <- convert_to_moment_kcat( subset( kcat_df , select = c( "react_id" , pred.names[i] ) ) )
          
          try_ret<- try(
            init_model <- build_initial_model(model = cond.model$model.org , mod2 = cond.model$mod2 , 
                                              Kappa = kcat , # use modelled kapmax
                                              mw = glob.parms$mw , # mw and RHS were  not saved in earlier versions of build_initial_model() and might htus not be found in mcmc_res
                                              RHS =  glob.parms$RHS ,
                                              moment_method = glob.parms$moment_method , solver = glob.parms$solver , medval =  NULL , 
                                              allow_kcat_impu = TRUE
            )
          )
          
          fluxes_li[[k]][[i]] <- data.frame( react_id = init_model$irrev_model@react_id , 
                                             init_model$sol$fluxes[ seq_along(init_model$irrev_model@react_id) ]  , stringsAsFactors = FALSE)
          
          if(class(try_ret) !=  "try-error"){
            gene_states <- get_gene_states(init_model = init_model )
            
            out_det <- rbind(out_det , 
                             data.frame( source = pred.names[i],
                                         growth_condition = condition_k ,
                                         gene_bnum =  gene_states$gene.bnums , #prot_abs_k$Bnumber ,
                                         abundance = gene_states$gene.dw.frac ), #proetin weight per dw (weight fraction)
                             stringsAsFactors = FALSE
            )
            
            
            schmidt.pred <- merge( gene_states , prot_abs_k  , by.x = "gene.bnums", by.y = "Bnumber"  )
            
            # normalize to total metabolic sum of proteins after merge 
            # gene.dw.frac is the weight fraction 
            schmidt.pred$prot.frac <- schmidt.pred$gene.dw.frac/sum(schmidt.pred$gene.dw.frac)
            # normalize proteomics data
            schmidt.pred[condition_k] <- schmidt.pred[condition_k] / sum(schmidt.pred[condition_k], na.rm = TRUE)
            
            # remove NAs as well taht originate from non-matched genes
            pred.zero_log <- schmidt.pred$prot.frac<=0 | is.na( schmidt.pred$prot.frac )
            obs.zero_log <- schmidt.pred[condition_k]<=0 | is.na( schmidt.pred[condition_k] )
            
            sum(!pred.zero_log & !obs.zero_log)
            
            pred <- log10(schmidt.pred$prot.frac[!pred.zero_log & !obs.zero_log])
            obs <- log10(schmidt.pred[condition_k][!pred.zero_log & !obs.zero_log, ])
            
            out_li[[k]]$R2[i] <- cor( pred , obs , use = "pairwise.complete.obs") ^2
            out_li[[k]]$RMSE[i] <- sqrt( mean( ( pred - obs )^2  , na.rm = TRUE ) )
            out_li[[k]]$n[i] <- sum( !is.na(pred) & !is.na(obs) )
            out_li[[k]]$obj[i] <- init_model$sol$obj
            
            genes_used_li[[k]][[i]] <- as.character( schmidt.pred$gene.bnums[ !pred.zero_log & !obs.zero_log ] )
          }else{
            out_li[[k]]$R2[i] <- NA
            out_li[[k]]$RMSE[i] <- NA
            out_li[[k]]$n[i] <- NA
            out_li[[k]]$obj[i] <- NA
            
            genes_used_li[k][i] <- NA
            genes_used_li[[k]][[i]] <- NA
          }
          
        }
      }
      
      out_df <- do.call( rbind , out_li )
      
      return(list( out_df = out_df , 
                   genes_used_li = genes_used_li, 
                   out_det = out_det , 
                   fluxes_li = fluxes_li
                   )
      )
}

convert_to_moment_kcat <- function(df){
  # add direction column
  df$dirxn <- ifelse( grepl("_b$",df$react_id) , -1, 1 )
  # remove suffixes from ids
  df$react_id <- sub( "_b$|_f$" , "" , df$react_id )
  names(df) <- c( "rxn_id" , "val" , "dirxn" )
  return(df)
}

# ensure that all keff vectors use the same set of genes in each condition
intersect_mechMod_preds <- function(ab_df, # abundace dataframe contianing columns "source" , "growth_condition" , "gene_bnum" , "abundance" , "mech_model"
                                    data_source_name = "Schmidt_data", # abundace dataframe contianing columns "source" , "growth_condition" , "gene_bnum" , "abundance" , "mech_model" for experimental data
                                    custom_gene_set_li # list of genes (with names that match conditions) that are used for comparison with experimental data
){
  
  # filter zer oabundance and NAs
  ab_df <- ab_df %>% filter( abundance>0 & !is.na(abundance) )
  
  # extract data
  ab_data_df <- ab_df %>% 
    filter( source %in% c( data_source_name  )
    ) 
  # remove data from main df
  ab_df <- ab_df %>% 
    filter( ! source %in% c( data_source_name  )
    )
  
  # investigate duplicates in predictions
  dup_test <- ab_df %>% filter( source != data_source_name)  %>% 
    group_by(source , growth_condition, mech_model ) %>% 
    summarise(n_dups = sum(duplicated(gene_bnum))) %>% as.data.frame()
  # no occurences of dups
  stopifnot(all(dup_test$n_dups == 0) )
  
  # now get condition-specific overlap
  conds <- unique(ab_df$growth_condition)
  
  if(! missing(custom_gene_set_li)){stopifnot(all(conds %in% names(custom_gene_set_li) )) }
  
  source <- unique(ab_df$source)
  
  mech_models <- unique(ab_df$mech_model)
  
  cond_spec_overlap <- vector(mode = "list" , length = length(conds) )
  names(cond_spec_overlap) <- conds
  
  cond_spec_overlap_ME <- vector(mode = "list" , length = length(conds) )
  names(cond_spec_overlap_ME) <- conds
  
  # create output structures
  out_li <- vector( mode = "list" , length = length(conds) )
  names(out_li) <- conds
  out_li <- lapply(out_li , function(x)  data.frame( "source" = source , "R2" = NA , "RMSE" = NA , "n" = NA , "obj" = NA ) )
  for(i in seq_along( out_li ) ){
    out_li[[i]]$condition <- names(out_li)[i]
  }
  
  
  out_li_ME <- vector( mode = "list" , length = length(conds) )
  names(out_li_ME) <- conds
  out_li_ME <- lapply(out_li_ME , function(x)  data.frame( "source" = source , "R2" = NA , "RMSE" = NA , "n" = NA , "obj" = NA ) )
  for(i in seq_along( out_li_ME ) ){
    out_li_ME[[i]]$condition <- names(out_li_ME)[i]
  }
  
  # save observations usedd for performance comparison for scatterplots
  ab_df_overlapped <- data.frame( source =character() ,
                               growth_condition =character() ,
                               gene_bnum    =character() ,
                               weight_frac_predicted = numeric() ,
                               weight_frac_observed = numeric() ,
                               mech_model =character(), 
                               stringsAsFactors = FALSE) 
  
  for(k in seq_along(mech_models)){
    mech_model_k <- mech_models[k]
    
    ab_df_k <- ab_df %>% filter(mech_model == mech_model_k)
    
    for ( i in seq_along(conds) ){
      if( conds[i] %in% ab_df_k$growth_condition){
        
        ab_df_k_i <- ab_df_k %>% filter(growth_condition == conds[i])
        ab_data_df_i <- ab_data_df %>% filter(growth_condition == conds[i]) # also filter data or growth condition
        
        # add data to mech_model- and condition-specific predicitons
        ab_df_k_i <- rbind( ab_df_k_i , ab_data_df_i )
        
        list_of_bnums_i <- split( ab_df_k_i$gene_bnum ,  ab_df_k_i$source , drop = TRUE )
        
        if( missing(custom_gene_set_li) ){
          bnum_overlap_i <- Reduce( intersect , list_of_bnums_i )
        }else{
          bnum_overlap_i <- custom_gene_set_li[[ conds[i] ]]
        }
        
        if(mech_model_k=="MOMENT"){
          cond_spec_overlap[[i]] <- bnum_overlap_i
        }else if (mech_model_k=="ME"){
          cond_spec_overlap_ME[[i]] <- bnum_overlap_i
        }else{
          stop("ab_df$mech_model needs to be \"MOMENT\" or \"ME\"")
        }
        
        for( j in seq_along(source) ){
          source_j <- source[j]
          
          # select correct sources 
          obs_df <- ab_data_df_i[, c("gene_bnum", "abundance")]          
          pred_df <- ab_df_k_i[ab_df_k_i$source == source_j  , c("gene_bnum", "abundance")]          
          
          # data has duplicated bnumbers 
          anyDuplicated( obs_df$gene_bnum )
          anyDuplicated( pred_df$gene_bnum )
          
          # remove dups and missing bnums
          obs_df <- obs_df [ !duplicated( obs_df$gene_bnum ) & obs_df$gene_bnum != "" ,  ]
          
          obs_pred_df <- merge( obs_df , pred_df , by.x =  "gene_bnum", by.y = "gene_bnum" , suffixes = c("_obs", "_pred") )
          
          # normalize by sum before focussing on condition-overlap only
          
          #metabolic proteome weight fractions
          obs_pred_df$met_frac_obs <- obs_pred_df$abundance_obs / sum(obs_pred_df$abundance_obs)
          obs_pred_df$met_frac_pred <- obs_pred_df$abundance_pred / sum(obs_pred_df$abundance_pred)
          
          # filter for overlpping bnums only
          obs_pred_df_ol <- obs_pred_df [ obs_pred_df$gene_bnum %in% bnum_overlap_i , ] 
          
          # save for scatterplots
          ab_df_overlapped <- rbind( ab_df_overlapped , data.frame( source = source_j , 
                                                              growth_condition = conds[i] , 
                                                              gene_bnum = obs_pred_df_ol$gene_bnum ,
                                                              weight_frac_predicted = obs_pred_df_ol$met_frac_pred ,
                                                              weight_frac_observed = obs_pred_df_ol$met_frac_obs  ,
                                                              mech_model = mech_model_k  ) )
          
          # pearson's R
          if(mech_model_k=="MOMENT"){
            out_li[[i]]$R2[j] <- cor( log10(obs_pred_df_ol$met_frac_pred) , log10(obs_pred_df_ol$met_frac_obs) , use = "pairwise.complete.obs") ^2
            out_li[[i]]$RMSE[j] <- sqrt( mean( ( log10(obs_pred_df_ol$met_frac_pred) - log10(obs_pred_df_ol$met_frac_obs) )^2  , na.rm = TRUE ) )
            out_li[[i]]$n[j] <- sum( !is.na(obs_pred_df_ol$met_frac_pred) & !is.na(obs_pred_df_ol$met_frac_obs) )
            
          }else if(mech_model_k=="ME"){
            out_li_ME[[i]]$R2[j] <- cor( log10(obs_pred_df_ol$met_frac_pred) , log10(obs_pred_df_ol$met_frac_obs) , use = "pairwise.complete.obs") ^2
            out_li_ME[[i]]$RMSE[j] <- sqrt( mean( ( log10(obs_pred_df_ol$met_frac_pred) - log10(obs_pred_df_ol$met_frac_obs) )^2  , na.rm = TRUE ) )
            out_li_ME[[i]]$n[j] <- sum( !is.na(obs_pred_df_ol$met_frac_pred) & !is.na(obs_pred_df_ol$met_frac_obs) )
          }
        }
      }
    }
  }
  
  # collect perfromance stats to df
  out_df_MOMENT <- Reduce(  rbind , out_li )
  out_df_MOMENT$mech_model <- "MOMENT"
  
  out_df_ME <- Reduce(  rbind , out_li_ME )
  out_df_ME$mech_model <- "ME"
  
  out_df <- rbind(out_df_MOMENT , out_df_ME)
  
  return(list(out_df=out_df , cond_spec_overlap = cond_spec_overlap , 
              cond_spec_overlap_ME = cond_spec_overlap_ME , 
              ab_df_overlapped = ab_df_overlapped ) )
}

