
GetSim <- function(drugName=DrugName, DRS= drugResponseScore,GES= GeneEffectScores ){
  ## make sure that users already use the common cell lines in both datasets >> the result will be precise.
  common.c=intersect(colnames(DRS),
                     colnames(GES))

  DRS.c=DRS[drugName,common.c ]
  GES.c=GES[,common.c]
  ## calculate corelation.
  out.Sim.GES.DRS  =sapply(1:nrow(GES.c),function(x)
    unlist(cor.test_trimmed_v0.default(DRS.c, GES.c[x,])))
  out.Sim.GES.DRS =t(out.Sim.GES.DRS)
  row.names(out.Sim.GES.DRS  ) <- row.names(GES.c)
  FDR <- fdrCor(out.Sim.GES.DRS[,1])
  out.Sim.GES.DRS <- cbind(out.Sim.GES.DRS, FDR=FDR)
 # head(out.Sim.GES.DRS)
  out.Sim.GES.DRS
}
