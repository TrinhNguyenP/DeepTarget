## this is used to check w the original data provided by Sanju
## common samples is obtain inside the functions. there is difference between getting common celllines and selecing unique drug.
rm(list=ls());
##
library ( DeepTarget)
##data ("OntargetM")
onTarget.O=readRDS('../../../12_07_22/drug.OnTarget-master/Data/onTarget_v3.0.RDS')
S.Drug <- c('K70301465','K09951645')
## check to see whether it has duplicated drug.
## prepare data
KO.GES <- onTarget.O$avana
sec.prism <- onTarget.O$secondary_prism
d.mt= onTarget.O$mutations_matrix
d.expr <- onTarget.O$expression_20Q4
metadata <- onTarget.O$drugCategory[,c("broad_id_trimmed", 'name','target', 'drug_category','moa') ]
sec.prism.f <- sec.prism[ which ( row.names(sec.prism)  %in% S.Drug),]
## because we use mapping to calculate the correlation so we have to feed in the righ one.
## let's do the one with least NA, and most NA as we selected unique drug.
rowSums(is.na(sec.prism.f[1:2,]))
rowSums(is.na(sec.prism.f))
## 1.3 has higher NA. 2.4 has less NA
#sec.prism.f.NA.more <- sec.prism.f[c(1,3),]
sec.prism.f.NA.less <- sec.prism.f[c(2,4),]

## Compute a correlation between the every gene crispr KO vs each drug response

# List.sim=mclapply(S.Drug,function(x) GetSim(x,DRS=sec.prism.f.NA.more, GES=KO.GES),mc.cores = 2)

List.sim=mclapply(S.Drug,function(x) GetSim(x,DRS=sec.prism.f.NA.less, GES=KO.GES),mc.cores = 2)

names(List.sim)=S.Drug
#saveRDS(List.sim,
#       file = paste('drugVScrispr.Similarity.List.RDS', Sys.Date(), '.RDS', sep=''))
## Map to get the stats from the list of simalirty to obtain the best and max corelation.

DrugTargetSim <- PredTarget(Sim.GES.DRS=List.sim, D.M = metadata)
DrugGeneMaxSim <- PredMaxSim(Sim.GES.DRS=List.sim, D.M = metadata)

## note: this is mapped w the overlaped samples in the files associated wt the function.( Same way w Sanju)
MutantInteract <- DoInteractMutant (Predtargets=DrugTargetSim,Mutant=d.mt,DRS=sec.prism.f.NA.less,GES=KO.GES)
## assign the estimate as strength, and P val from the interaction model.
TargetMutSpecificity <- data.frame(MaxTgt_Inter_Mut_strength=sapply(MutantInteract, function(x) x[1]), MaxTgt_Inter_Mut_Pval=sapply(MutantInteract, function(x) x[2]))
row.names(TargetMutSpecificity) <- row.names(DrugTargetSim)
## note: this is mapped w the overlaped samples in the files associated wt the function.( Same way w Sanju)
ExpInteract <- DoInteractExp (DrugTargetSim,d.expr,sec.prism.f.NA.less,KO.GES,CutOff = 2)
TargetExpSpecificity <- data.frame(MaxTgt_Inter_Exp_strength=sapply(ExpInteract, function(x) x[1]), MaxTgt_Inter_Exp_Pval=sapply(ExpInteract, function(x) x[2]))
row.names(TargetExpSpecificity) <- row.names(DrugTargetSim)
identical ( DrugTargetSim[,1],DrugGeneMaxSim[,1])
all(sapply(list(row.names(DrugTargetSim),
                row.names(TargetMutSpecificity)), FUN = identical, row.names(TargetExpSpecificity)))
Pred.d <- cbind ( DrugTargetSim,DrugGeneMaxSim,TargetMutSpecificity,TargetExpSpecificity)

## whether interaction is true or false based on cut-off. estimate and p val from lm model
Pred.d$Whether_interaction_Ex_based=sapply(ExpInteract, function(x) x[1]<0 & x[2]<0.2 )
## mutation interaction with P <0.1
Pred.d$predicted_resistance_mutation = Pred.d$MaxTgt_Inter_Mut_Pval<0.1
head(Pred.d )
write.csv (Pred.d, "../../../Testing/Checking //Pred.lessNA.csv" )
## for less NA
####
## the figures.
## Drug mutant binding.
DOI = 'dabrafenib'
GOI ='BRAF'
DMB (DN=DOI,GN=GOI,Pred=Pred.d,Mutant=d.mt,DRS= sec.prism.f.NA.less,GES= KO.GES)
## figure 3.a predction of secondary targets of drugs that mediate the eresponse.
DOI = 'ibrutinib'
GOI ='BTK'
DTR ( DN=DOI,GN=GOI,Pred=Pred.d, Exp=d.expr,DRS= sec.prism.f.NA.less,GES= KO.GES,CutOff= 2)
DOI = 'ibrutinib'
GOI ='ERBB2'
DTR ( DN=DOI,GN=GOI,Pred=Pred.d, Exp=d.expr,DRS= sec.prism.f.NA.less,GES= GES,CutOff= 2)
############ second case.

## the more NA.
## this is used to check w the orginal data provided by Sanju
rm(list=ls());
library ( DeepTarget)
##data ("OntargetM")
onTarget.O=readRDS('../../../12_07_22/drug.OnTarget-master/Data/onTarget_v3.0.RDS')
S.Drug <- c('K70301465','K09951645')
## check to see whether it has duplicated drug.

## prepare data
KO.GES <- onTarget.O$avana
sec.prism <- onTarget.O$secondary_prism
d.mt= onTarget.O$mutations_matrix
d.expr <- onTarget.O$expression_20Q4
metadata <- onTarget.O$drugCategory[,c("broad_id_trimmed", 'name','target', 'drug_category','moa') ]
sec.prism.f <- sec.prism[ which ( row.names(sec.prism)  %in% S.Drug),]
## because we use mapping to calculate the correlation so we have to feed in the righ one.
## let's do the one with least NA, and most NA as we selected unique drug.
rowSums(is.na(sec.prism.f[1:2,]))
rowSums(is.na(sec.prism.f))
## 1.3 has higher NA. 12.4 has less NA
sec.prism.f.NA.more <- sec.prism.f[c(1,3),]


## Compute a correlation between the every gene crispr KO vs each drug response

# List.sim=mclapply(S.Drug,function(x) GetSim(x,DRS=sec.prism.f.NA.more, GES=KO.GES),mc.cores = 2)

List.sim=mclapply(S.Drug,function(x) GetSim(x,DRS=sec.prism.f.NA.more , GES=KO.GES),mc.cores = 2)

names(List.sim)=S.Drug
#saveRDS(List.sim,
#       file = paste('drugVScrispr.Similarity.List.RDS', Sys.Date(), '.RDS', sep=''))
## Map to get the stats from the list of simalirty to obtain the best and max corelation.

DrugTargetSim <- PredTarget(Sim.GES.DRS=List.sim, D.M = metadata)
DrugGeneMaxSim <- PredMaxSim(Sim.GES.DRS=List.sim, D.M = metadata)

## note: this is mapped w the overlaped samples in the files associated wt the function.( Same way w Sanju)
MutantInteract <- DoInteractMutant (Predtargets=DrugTargetSim,Mutant=d.mt,DRS=sec.prism.f.NA.more ,GES=KO.GES)
## assign the estimate as strength, and P val from the interaction model.
TargetMutSpecificity <- data.frame(MaxTgt_Inter_Mut_strength=sapply(MutantInteract, function(x) x[1]), MaxTgt_Inter_Mut_Pval=sapply(MutantInteract, function(x) x[2]))
row.names(TargetMutSpecificity) <- row.names(DrugTargetSim)
## note: this is mapped w the overlaped samples in the files associated wt the function.( Same way w Sanju)
ExpInteract <- DoInteractExp (DrugTargetSim,d.expr,sec.prism.f.NA.more ,KO.GES,CutOff = 2)
TargetExpSpecificity <- data.frame(MaxTgt_Inter_Exp_strength=sapply(ExpInteract, function(x) x[1]), MaxTgt_Inter_Exp_Pval=sapply(ExpInteract, function(x) x[2]))
row.names(TargetExpSpecificity) <- row.names(DrugTargetSim)
identical ( DrugTargetSim[,1],DrugGeneMaxSim[,1])
all(sapply(list(row.names(DrugTargetSim),
                row.names(TargetMutSpecificity)), FUN = identical, row.names(TargetExpSpecificity)))
Pred.d <- cbind ( DrugTargetSim,DrugGeneMaxSim,TargetMutSpecificity,TargetExpSpecificity)

## whether interaction is true or false based on cut-off. estimate and p val from lm model
Pred.d$Whether_interaction_Ex_based=sapply(ExpInteract, function(x) x[1]<0 & x[2]<0.2 )
## mutation interaction with P <0.1
Pred.d$predicted_resistance_mutation = Pred.d$MaxTgt_Inter_Mut_Pval<0.1
head(Pred.d )
write.csv (Pred.d, "../../../Testing/Checking //Pred.more.NA.csv" )
## the figures.
## Drug mutant binding.
DOI = 'dabrafenib'
GOI ='BRAF'
DMB (DN=DOI,GN=GOI,Pred=Pred.d,Mutant=d.mt,DRS= sec.prism.f.NA.more,GES= KO.GES)
## figure 3.a predction of secondary targets of drugs that mediate the eresponse.
DOI = 'ibrutinib'
GOI ='BTK'
DTR ( DN=DOI,GN=GOI,Pred=Pred.d, Exp=d.expr,DRS= sec.prism.f.NA.more,KO.GES,CutOff= 2)
DOI = 'ibrutinib'
GOI ='EGFR'
DTR ( DN=DOI,GN=GOI,Pred=Pred.d, Exp=d.expr,DRS= sec.prism.f.NA.more,KO.GES,CutOff= 2)






