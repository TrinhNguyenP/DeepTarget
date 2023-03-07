
## perform the interaction based on cut-off

DoInteractMutant <- function(Predtargets=Drug.targeted.sim,Mutant=Mutant,DRS=DRS, GES=GES )
  {
  interactFeatures=lapply(1:nrow(Predtargets), function(x)
  {
    ## drug response based on the drug ID
    DRS.f =errHandle(DRS[Predtargets[x,1],])
    ## Gene effect scores based on the targeted gene
    GES.f=errHandle(GES[Predtargets[x,3],])
    ## groups based on the cut-off. Lower and higher based on the targed gene
    Mutant.Group =errHandle(Mutant[Predtargets[x,3],] )
    ## True is low, False is high.
    errHandle(summary(lm(DRS.f ~
                             GES.f * Mutant.Group))$coefficients[4,c(1,4)])

    }
)

  names(interactFeatures)=Predtargets$drugName
  interactFeatures
}


