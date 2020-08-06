load("data/Moment_data/mw_iML1515.rdata")
mw_iML1515$mw <- mw_iML1515$mw/1000
mw <- mw_iML1515

# initial Kappa:
init.kcat <-  23.8 
# Adadi cites that .56 g/gDW of the Ecoli cell are in Proteins, and fitted that 48% of that are used in metabolism, 0.56*0.48~= 0.27 = RHS. Proteomaps shows about 57% of the proteome to be used in metabolism, this would result in: 0.56*0.57~= 0.32 = RHS
RHS <- 0.32 
moment_method <- "moment"
solver <- "cplexAPI" 
