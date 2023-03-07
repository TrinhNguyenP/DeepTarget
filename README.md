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
## set seed
set.seed (12345)
### Drug IDs
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
## 5 random drugs
S.Drugs <- sample(All.Drugs, 5)

## Gene effect scores from KO method
KO.GES <- OntargetM$avana_CRISPR
### viablity scores from drug treament
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

