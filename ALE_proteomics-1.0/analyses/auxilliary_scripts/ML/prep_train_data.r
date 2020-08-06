# prepare train set:
# input: 
# -  destination of file with features 
# - model obect rdata file that is used to remove spontaneous reactions
# -  destination of kappmax df with reaction-specific bigg IDs
#
# output:
# - train data (still contains reactions with NA labels, but no spontaneous or exchange reactions)

prep_train_data <- function(feature_file = "../kcat_ML/data/loc_and_glob_struct_and_ECs_use_bind_sit_FALSE_MFA-FBA_bm_norm_FALSE_use_met_dG_Km_TRUE_iML1515.rdata",
                            model_file = "data/M-models/iML1515/iML1515.rdata" , 
                            output_file) {
  require(sybil)
  
  ld_env <- new.env()
  load( feature_file , envir = ld_env )
  X <- ld_env$glob
  
  mod_env <- new.env()
  load( model_file , envir = mod_env )
  stopifnot( length(ls(mod_env)) == 1 ) # we make use of the fact  that the model file only contains one object
  irr_mod <- mod2irrev(mod_env$iML1515)
  spont_reacts <- irr_mod@react_id[irr_mod@gpr=="s0001"]
  spont_reacts <- sub( "_f$", "" ,spont_reacts)
  
  # replace _f in forward reactions
  X$bigg_id <- sub( "_f$", "" ,X$bigg_id)
  
  # remove exchange reacts
  X <- X[ ! X$exchange_log , ]
  # remove spont reacts
  X <- X[ ! X$bigg_id %in% spont_reacts , ]
  
  ######################
  # rename features
  nam <- names(X)
  nam[ nam == "seq_hydrophobicity.kd_AS_avg" ] <- "AS_hydrophobicity"
  nam[ nam == "struct_CA_DEPTH.msms_AS_avg" ] <- "AS_alpha_C_depth"
  nam[ nam == "struct_ASA.dssp_AS_avg" ] <- "AS_solvent_access"
  nam[ nam == "seq_SS.sspro_AS_maj" ] <- "AS_maj_sec_struct"
  nam[ nam == "seq_RSA.accpro_AS_maj" ] <- "AS_exposed"
  names(X) <- nam
  
  
  # impute additional features
  X$m_efficiency[ is.na(X$m_efficiency) ] <- median(X$m_efficiency , na.rm = TRUE)
  X$log10_m_subst_conc[ is.na(X$log10_m_subst_conc) ] <- median(X$log10_m_subst_conc , na.rm = TRUE)
  X$log10_m_prod_conc[ is.na(X$log10_m_prod_conc) ] <- median(X$log10_m_prod_conc , na.rm = TRUE)
  X$log10_Km_m[ is.na(X$log10_Km_m) ] <- median(X$log10_Km_m , na.rm = TRUE)
  
  ################
  # add labels (output data)
  Y <- readRDS(file = output_file)
  stopifnot( "bigg_id" %in% names(Y) )
  
  print(paste("n intersecting bigg ids :" , length(intersect(Y$bigg_id , X$bigg_id))))
  print(paste("n label rows :" , nrow(Y) ))
  
  ##################
  # merge X and Y
  
  XY <- merge( X, Y , by = "bigg_id" , all = TRUE )
  
  rownames(XY)  <- XY$bigg_id
  
  # select 
  XY <- subset( XY , select = - c( bigg_id ,  # pH and temp not releavant for in vivo kcats
                                 exchange_log , react_cyto_log,
                                 sasa,alpha_frac,beta_frac
    )
  )
  
  return(XY)
}