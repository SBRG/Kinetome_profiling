library(sybil);library(sybilccFBA)


# set Parameters for initializing the model:
load( file = "data/M-models/iML1515/iML1515.rdata" )
model.org <- iML1515

# read JSON for LB medium
library(rjson)
lb <- rjson::fromJSON( file = "data/M-models/iML1515/Media/LB_Media.json")

names(lb)[ names(lb) == "EX_dad__2_e" ] <- "EX_dad_2_e"


class(lb)
for( i in seq_along(lb) ){
  reac_id_bigg <- names(lb)[i]
  reac_id_sybil <- sub( "_e$" , "(e)" , reac_id_bigg )
  print(reac_id_sybil)
  stopifnot( reac_id_sybil %in% react_id(model.org) )
  model.org@react_rev[ react_id(model.org)==reac_id_sybil ] <- TRUE
  lowbnd(model.org)[react_id(model.org)==reac_id_sybil] = -1000#lb[[i]]
}

  

mod2=mod2irrev(model.org)

# to match aerobic glucose growth
uppbnd(mod2)[react_id(mod2)=="EX_glc__D(e)_b"]=0
lowbnd(mod2)[react_id(mod2)=="EX_cbl1(e)"] = -1000 # From Orth  et al 2011: "Because only a very small amount of B12 is required for growth, the lower bound on cob(I)alamin uptake is arbitrary and never actually constraining in practice."
#uppbnd(mod2)[react_id(mod2)=="EX_glyc(e)"] = 0
#uppbnd(mod2)[react_id(mod2)=="EX_ac(e)"] = 0
uppbnd(mod2)[react_id(mod2)=="EX_o2(e)_b"] = 1000
lowbnd(mod2)[react_id(mod2)=="ATPM"] = 0



b <- findExchReact(mod2)
cbind(b@react_id[b@uptake],b@uppbnd[b@uptake])
