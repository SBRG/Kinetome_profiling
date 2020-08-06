##############################################################################
# returns a list that extends the output of original cfba_moment*() functions
library(sybil)
library(sybilccFBA)


build_initial_model <- function( model, # of class "modelorg"
                                 mod2 , # e.g. mod2irrev(model)
                     Kappa, # char matrix with 3 cols, rxn_id, va, dirxn (as used in sybilccFBA)
                     mw, # df with 2 columns
                     RHS = 0.32, # the proteome budget. [g/gDW]
                     # Adadi cites that .56 g/gDW of the Ecoli cell are in Proteins, and fitted that 48% of that are used in metabolism, 0.56*0.48~= 0.27 = RHS
                     # Proteomaps shows about 57% of the proteome to be used in metabolism, this would result in: 0.56*0.57~= 0.32 = RHS
                     moment_method = "moment", # or "moment_mr". the flavor of moment to be used in all downstream analyses
                     solver = "glpkAPI",
                     cplex_parms_li =NULL , #named (integers as found in ?cplexConstants ) list of cplex parameters
                     medval = 3600*22.6 , 
                     verbose = 0,
                     allow_kcat_impu = FALSE,
                     gene.lb = 0 # default lower bound used for gene variables
){
  ###
  
  if(moment_method == "moment"){
    # source modified version of ccfbamoment. currently the only difference is the return of a vector to keep track of the row of gs (in the problem matrix) that are associated with reactions 
    source("analyses/auxilliary_scripts/MOMENT/cfba_moment_1-1.R")
    mom_ret <- cfba_moment_dh_1.1( model = model , mod2 = mod2 , Kcat = Kappa , MW = mw , verbose=verbose , RHS=RHS , solver=solver , cplex_parms_li =  cplex_parms_li , medval= medval, allow_kcat_impu = allow_kcat_impu , gene.lb = gene.lb) 
  }else if (moment_method == "moment_mr"){
    stop("_mr not yet implemented")
  }else{ stop("moment_method should be either \"moment\" or \"moment_mr\"") }
  
  if(mom_ret$sol$obj==0){warning("Problem in build_initial_model(): model does not support growth.")}
  
  return( c( mom_ret , list( "model" = model , 
                             "irrev_model" = mod2, 
                             "moment_method" = moment_method ,
                             "Kappa" = Kappa , 
                             "solver" = solver , 
                             "cplex_parms_li" = cplex_parms_li,
                             "allow_kcat_impu" = allow_kcat_impu ,
                             "mw" = mw , 
                             "RHS" = RHS, 
                             "medval" =  medval
  )
  )
  )
}

get_gene_states <- function(init_model){
      gpr_i <- gpr(init_model$model)
      gene.bnums <- unique(unlist(regmatches(gpr_i,gregexpr("b\\d{4}",gpr_i)) ))
      gene.mws <- init_model$mw[ match( gene.bnums , init_model$mw$Synonym ) , ]$mw  #[g/mmol]
      
      gene.cols <- as.numeric( init_model$geneCol[ match( gene.bnums , init_model$geneCol[,"gene"] ) , "Col" ]  )
      gene.state <- init_model$sol$fluxes[gene.cols] # [mmol/gDW] states of the gene variables associated with the current reaction
      
      gene.dw.frac_li <- as.list(gene.state*gene.mws)
      names(gene.dw.frac_li) <- gene.bnums
      
      ret_df <- data.frame(gene.bnums=gene.bnums , 
                           gene.mws=gene.mws , 
                           gene.state=gene.state , 
                           gene.dw.frac = unlist(gene.dw.frac_li)   )
      
      return( ret_df )
}


