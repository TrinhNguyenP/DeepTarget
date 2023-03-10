## this is used to pull out the targeted gene by the drug that has the max corelation.
## need to check this.

PredTarget <- function(Sim.GES.DRS=sim, D.M = Drug.Metadata){
    ## need to test this function based on 1 drug, 2 drugs.
    ## Make sure the user use the same meta data for the list of Sim.GES.DRS
    Drug.Id.i <- match(names(Sim.GES.DRS), D.M[,1])
    D.M.f <- D.M[Drug.Id.i,]
    ## based on the name of the drug.
    Splt.targets = str_split(as.character(D.M.f[,3]), ", ")
    ## this is extract P val, and cor values and turn to the matrix where genes are rows and drugs are column.
    corrMatPval=sapply(Sim.GES.DRS, function(x) x[,1])
    corrMat=sapply(Sim.GES.DRS, function(x) x[,2])
    corrMatFDR=sapply(Sim.GES.DRS, function(x) x[,3])
    # corrMatFDR=apply(corrMat_P, 2, function(x) fdrCor(x))
    ## map to get the name of drug with corelation values based on the genes.
    ## errHandle if can't find from the gene effect score list will return NA.
    ## x is the column of the drug.
    ## if there is only one drug. it din't work well.

    Stat.Drug=sapply(1:length(Splt.targets),
                                       function(x) {
                                           ret_cor=corrMat[
                                               match(Splt.targets[[x]],
                                                     rownames(corrMat)), x]
                                           names(ret_cor)=errHandle(Splt.targets[[x]])
                                           ret_cor
                                       } )

   ## if there is one drug.
     if ( nrow(D.M.f)==1){
         ## add errHandle function to avoid the eror of returing -Inf
        Target.Max.cor = errHandle(max(Stat.Drug,na.rm=T ))
        Target.Max.Name=names(Stat.Drug)[which.max(Stat.Drug)]
        ## drug Id will be the same as the metadata.
        Drug.id <- names(Sim.GES.DRS)
        names(Target.Max.cor) <- as.character(D.M.f[,2])
    }else{
        names(Stat.Drug) <- as.character(D.M.f[,2])
        Target.Max.cor=sapply(Stat.Drug, function(x) max(x, na.rm=T))
        Target.Max.Name=sapply(Stat.Drug,function(x) names(which.max(x)))
        Drug.id = D.M.f[match(names(Target.Max.cor), D.M.f$name),1]
    }
    #which ( row.names(corrMat)=="CACNA1C")

    ###
    ### Dataframe of the knowntarget_prediction ( get the annotation from broad.)

    ##  Pull the names of the targeted gene.

    ## if no name found from mapping, add NA.
    Target.Max.Name[sapply(Target.Max.Name, length)==0]=NA
    ### obtain drug ID ( make sure that we map the correct one)
### these assigned as known target due to the drug information.
    Target.Pred=data.frame(
        DrugID = Drug.id,
        drugName=names(Target.Max.cor),
        MaxTargetName=unlist(Target.Max.Name),
        Maxcorr=Target.Max.cor)
   #
    # Known Target Significance
    ## errror if there is one row having NA.
    ## pull out p val and FDR.
    Target.Pred$KnownTargetCorrP = sapply(1:nrow(Target.Pred), function(x)
        errHandle(corrMatPval[Target.Pred[x,3], Target.Pred[x,1]]) )
    # Known Target Significance - FDR corrected
    Target.Pred$KnownTargetCorrFDR = sapply(1:nrow(Target.Pred), function(x)
        errHandle(  corrMatFDR[Target.Pred[x,3], Target.Pred[x,1]]) )
    Target.Pred
}


