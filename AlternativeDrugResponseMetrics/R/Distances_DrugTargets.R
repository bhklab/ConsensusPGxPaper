#' Calculation of Distances to Drug Targets
#'
#' Drug targets retrieved from GDSC and CHEMBL
#'
#' @param metric ANOVA results corresponding to a specific metric
#' @return A dataframe with the drug, FDR gene (Gene that is identified as significant in FDR and has FDR<10%), Original target and distances
#' @export
#'
Distances_DrugTargets<-function(metric,write.File=T,FileName)
{
  outFile=paste("Calculated_Values/",FileName,".txt",sep="")

  AllFeatures<-c()

  #Extract those drug gene associations with FDR < 10% and get the list of genes which are either mutated or have Copy Number Variations (CNV)
  Metric_extract_FDR25<-metric[metric$ANOVA_FEATURE_FDR < 10,]
  for(i in 1:nrow(Metric_extract_FDR25))
  {
    Feature_extract<-sapply(as.vector(Metric_extract_FDR25$"FEATURE"[i]), function(x) unlist(strsplit(x, "_")))
    Feature_extract_clean<-Feature_extract[!(Feature_extract %in% c("gain","cna","loss","mut"))]
    Features_Drugs<-cbind.data.frame(Metric_extract_FDR25$"DRUG_NAME"[i],Feature_extract_clean)
    AllFeatures<-rbind(AllFeatures,Features_Drugs)
  }
  colnames(AllFeatures) <- c("Drug","FDR_Gene")

  #Merge FDR Genes with actual Drug Targets
  FDR_Gene_Targets_merged<-merge(AllFeatures,Orig_Target,by=c("Drug"))

  #Load Omnipath networks
  NW_extract<-OmniPathNW[,colnames(OmniPathNW) %in% c("Node1","Node2")]
  g=graph.data.frame(NW_extract,directed=FALSE)
  Alldist<-c()

  #Find the distance between the actual drug targets and FDR genes
  for(j in 1:nrow(FDR_Gene_Targets_merged))
  {
    Source<-as.vector(FDR_Gene_Targets_merged$FDR_Gene[j])
    results<-tryCatch({
      d<-distances(g,v=Source,to=as.vector(FDR_Gene_Targets_merged$Target[j]))
      dist_info<-cbind.data.frame(FDR_Gene_Targets_merged[j,], Path_length=d[,1])
      Alldist<-rbind(Alldist,dist_info)
    }, error = function(err) {
      #print(paste("MY_ERROR:  ",err))

    } )

  }
 if(write.File) (write.table(Alldist, file=outFile, sep="\t",row.names=FALSE, quote=FALSE))
 return(Alldist)
}
