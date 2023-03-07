# DeepTarget
## V2 Of GDeepTarget For Bioconductor

library (DeepTarget)

data (OntargetM)
### 5 random drugs
set.seed (12345)
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
S.Drugs <- sample(All.Drugs, 5)

## Gene effect scores from KO method
KO.GES <- OntargetM$avana_CRISPR
### viablity scores 
sec.prism <- OntargetM$secondary_prism
## ## Obtain the similarity between the viablity scores from drug target treatment vs Gene effect score from KO method
sim <- mclapply(S.Drugs,function(x) GetSim(x,DRS=sec.prism, GES=KO.GES),mc.cores = 2)
names(sim) <- S.Drugs
## meta data from the drug treatment
Meta.data <- OntargetM$DrugMetadata
## obtain the max similarity scores for the known targeted genes from the drugs
DrugTargetSim <- PredTarget(sim,Meta.data)
## expression
d.expr <- OntargetM$expression_20Q4
## finding the interaction of the viablity scores vs Gene effect score in term of high and low expression
ExpInteract <- DoInteractExp (DrugTargetSim,d.expr,sec.prism,KO.GES,CutOff = 2)
# mutant 
d.mt <- OntargetM$mutations_mat
## ## finding the interaction of the viablity scores vs Gene effect score in term of mutant vs non-mutant
MutantInteract <- DoInteractMutant (Predtargets=DrugTargetSim,Mutant=d.mt,DRS=sec.prism,GES=KO.GES)
## obtain the best similarity scores for the known targeted genes from the drugs
Drug.Gene.max.sim <- PredMaxSim(Sim.GES.DRS=sim, D.M = Meta.data)
## perfom the pathway analyis based on the drug meta data ranking by similarity between the viablity scores from drug target treatment vs Gene effect score from KO method 
Meta.data <- OntargetM$DrugMetadata
Pwy.Enr <- DoPWY(Sim.GES.DRS=sim,D.M = Meta.data)

