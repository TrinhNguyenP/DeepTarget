rm(list=ls());
library ( DeepTarget)
data ("OntargetM")

S.Drug <- c('K70301465','K09951645')

## gene effect score from KO
KO.GES <- OntargetM$avana_CRISPR
dim(KO.GES)
## Response scores from treatment
sec.prism <- OntargetM$secondary_prism

## Compute a correlation between the every gene crispr KO vs each drug response
List.sim=mclapply(S.Drug,function(x) GetSim(x,DRS=sec.prism, GES=KO.GES),mc.cores = 2)
names(List.sim)=S.Drug
#saveRDS(List.sim,
 #       file = paste('drugVScrispr.Similarity.List.RDS', Sys.Date(), '.RDS', sep=''))
## Map to get the stats from the list of simalirty to obtain the best and max corelation.
metadata <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(Sim.GES.DRS=List.sim, D.M = metadata)
DrugGeneMaxSim <- PredMaxSim(Sim.GES.DRS=List.sim, D.M = metadata)
d.mt <- OntargetM$mutations_mat
MutantInteract <- DoInteractMutant (Predtargets=DrugTargetSim,Mutant=d.mt,DRS=sec.prism,GES=KO.GES)
## assign the estimate as strength, and P val from the interaction model.
TargetMutSpecificity <- data.frame(MaxTgt_Inter_Mut_strength=sapply(MutantInteract, function(x) x[1]), MaxTgt_Inter_Mut_Pval=sapply(MutantInteract, function(x) x[2]))
row.names(TargetMutSpecificity) <- row.names(DrugTargetSim)
d.expr <- OntargetM$expression_20Q4
ExpInteract <- DoInteractExp (DrugTargetSim,d.expr,sec.prism,KO.GES,CutOff = 2)
TargetExpSpecificity <- data.frame(MaxTgt_Inter_Exp_strength=sapply(ExpInteract, function(x) x[1]), MaxTgt_Inter_Exp_Pval=sapply(ExpInteract, function(x) x[2]))
row.names(TargetExpSpecificity) <- row.names(DrugTargetSim)

identical ( DrugTargetSim[,1],DrugGeneMaxSim[,1])
all(sapply(list(row.names(DrugTargetSim),
                row.names(TargetMutSpecificity)), FUN = identical, row.names(TargetExpSpecificity)))
Pred.d <- cbind ( DrugTargetSim,DrugGeneMaxSim,TargetMutSpecificity,TargetExpSpecificity)
head(Pred.d )
## whether interaction is true or false based on cut-off. estimate and p val from lm model
Pred.d$Whether_interaction_Ex_based=sapply(ExpInteract, function(x) x[1]<0 & x[2]<0.2 )
## mutation interaction with P <0.1
Pred.d$predicted_resistance_mutation = Pred.d$MaxTgt_Inter_Mut_Pval<0.1
Pred.d
write.csv (Pred.d, "../../../Testing/Checking //Pred.cleanData.csv" )
## Drug mutant binding.
DOI = 'dabrafenib'
GOI ='BRAF'
DMB (DN=DOI,GN=GOI,Pred=Pred.d,Mutant=d.mt,DRS= sec.prism,GES= KO.GES)
## figure 3.a predction of secondary targets of drugs that mediate the eresponse.
DOI = 'ibrutinib'
GOI ='BTK'
DTR ( DN=DOI,GN=GOI,Pred=Pred.d, Exp=d.expr,DRS= sec.prism,KO.GES,CutOff= 2)
DOI = 'ibrutinib'
GOI ='EGFR'
DTR ( DN=DOI,GN=GOI,Pred=Pred.d, Exp=d.expr,DRS= sec.prism,KO.GES,CutOff= 2)



