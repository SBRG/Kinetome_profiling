# aparmeters for simulating iJO with 

library(sybil);library(sybilccFBA) # use for mw dat


# set Parameters for initializing the model:
load( file = "data/M-models/iML1515/iML1515.rdata" )
model.org <- iML1515

model.org@react_rev[ react_id(model.org)=="EX_ac(e)" ] <- TRUE
lowbnd(model.org)[react_id(model.org)=="EX_ac(e)"] = -1000

model.org@react_id[model.org@lowbnd < 0 ] 
a <- findExchReact(model.org)
a@react_id[a@uptake]


mod2=mod2irrev(model.org)

# to match aerobic glucose growth as we did it in AF1260:
uppbnd(mod2)[react_id(mod2)=="EX_glc__D(e)_b"]=0
lowbnd(mod2)[react_id(mod2)=="EX_cbl1(e)"] = -1000 # From Orth  et al 2011: "Because only a very small amount of B12 is required for growth, the lower bound on cob(I)alamin uptake is arbitrary and never actually constraining in practice."
#uppbnd(mod2)[react_id(mod2)=="EX_glyc(e)"] = 0
#uppbnd(mod2)[react_id(mod2)=="EX_ac(e)_b"] = 1000
uppbnd(mod2)[react_id(mod2)=="EX_o2(e)_b"] = 1000
lowbnd(mod2)[react_id(mod2)=="ATPM"] = 0

b <- findExchReact(mod2)
b@react_id[b@uptake]


