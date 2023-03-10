# DeepTarget
## V2 Of GDeepTarget For Bioconductor
## install the devtools prior install the package.
install.packages("devtools")
## call the library
library(devtools)
## intall Deeptarget
install_github("TrinhNguyenP/DeepTarget")
## call the deeptarget
library (DeepTarget)
## data 
data (OntargetM)

## use two drug as an example.
S.Drugs <- c('K70301465','K09951645')
## Gene effect scores from KO method
KO.GES <- OntargetM$avana_CRISPR
### viablity scores from drug treament
sec.prism <- OntargetM$secondary_prism
## ## Obtain the similarity between the viablity scores from drug target treatment vs Gene effect score from KO method
sim <- mclapply(S.Drugs,function(x) GetSim(x,DRS=sec.prism, GES=KO.GES),mc.cores = 2)
##
names(sim) <- S.Drugs
## meta data from the drug treatment
Meta.data <- OntargetM$DrugMetadata
## obtain the max similarity scores for the known targeted genes from the drugs
DrugTargetSim <- PredTarget(sim,Meta.data)
## expression
d.expr <- OntargetM$expression_20Q4
## finding the interaction of the viablity scores vs Gene effect score in term of high and low expression
ExpInteract <- DoInteractExp (DrugTargetSim,d.expr,sec.prism,KO.GES,CutOff = 2)
###
TargetExpSpecificity <- data.frame(MaxTgt_Inter_Exp_strength=sapply(ExpInteract, function(x) x[1]), MaxTgt_Inter_Exp_Pval=sapply(ExpInteract, function(x) x[2]))
###
row.names(Target_Exp_specificity) <- row.names(DrugTargetSim)
# mutant 
d.mt <- OntargetM$mutations_mat
## ## finding the interaction of the viablity scores vs Gene effect score in term of mutant vs non-mutant
MutantInteract <- DoInteractMutant (Predtargets=DrugTargetSim,Mutant=d.mt,DRS=sec.prism,GES=KO.GES)
## assign the estimate as strength, and P val from the interaction model.
TargetMutSpecificity <- data.frame(MaxTgt_Inter_Mut_strength=sapply(MutantInteract, function(x) x[1]), MaxTgt_Inter_Mut_Pval=sapply(MutantInteract, function(x) x[2]))
#####
row.names(Target_Mut_specificity) <- row.names(DrugTargetSim)
## obtain the best similarity scores for the known targeted genes from the drugs
Drug.Gene.max.sim <- PredMaxSim(Sim.GES.DRS=sim, D.M = Meta.data)
#########
identical ( DrugTargetSim[,1],Drug.Gene.max.sim[,1] )
all(sapply(list(row.names(DrugTargetSim),
                row.names(TargetMutSpecificity)), FUN = identical, row.names(TargetExpSpecificity)))
 #####
 Pred.d <- cbind ( DrugTargetSim,DrugGeneMaxSim,TargetMutSpecificity,TargetExpSpecificity)
## whether interaction is true or false based on cut-off. estimate and p val from lm model
Pred.d$Whether_interaction_Ex_based=sapply(ExpInteract, function(x) x[1]<0 & x[2]<0.2 )
## mutation interaction with P <0.1
Pred.d$predicted_resistance_mutation = Pred.d$MaxTgt_Inter_Mut_Pval<0.1
####### obtain the best similarity scores f
Low.Exp = sapply(Pred.d[,3],function(x)errHandle(sum(d.expr[x,] < 2)) )
### Use cutoff<2
Pred.d$Low.Exp.Group <- Low.Exp
####### Plot the interaction bt viablilty and resonse to the drug for WT and mutant form.
DOI = 'dabrafenib'
GOI = 'BRAF'
DMB (DN=DOI,GN=GOI,Pred=Pred.d,Mutant=d.mt,DRS= sec.prism,GES= KO.GES)
###### Plot the interaction bt viablilty and resonse to the drug for low and high expression
DOI = 'ibrutinib'
GOI ='BTK'
DTR ( DN=DOI,GN=GOI,Pred=Pred.d, Exp=d.expr,DRS= sec.prism,GES= KO.GES,CutOff= 2)
## perfom the pathway analyis based on the drug meta data ranking by similarity between the viablity scores from drug target treatment vs Gene effect score from KO method 
Pwy.Enr <- DoPWY(Sim.GES.DRS=sim,D.M = Meta.data)
### 
## ## Obtain the similarity between the viablity scores from drug target treatment vs Gene effect score from KO method in term of low expression group
sim.lowExpr <- mclapply(S.Drugs,function(x) GetSim(x,DRS=sec.prism, GES=KO.GES),mc.cores = 2)
##

