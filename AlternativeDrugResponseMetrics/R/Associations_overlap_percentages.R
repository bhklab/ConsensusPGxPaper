#' Association overlap between different metrics
#'
#' Calculation of overlaps between the conventional and GR/NDR metrics based on varying thresholds
#'
#'
#'
#' @param  metric1 Conventional metric
#' @param  metric2 GR/NDR metric
#' @return A dataframe with the Association thresholds and Percentage overlap
#' @export
#'
#'
Associations_overlap_percentages<-function(metric1,metric2,metric1_name,metric2_name)
{
  #Extract the top 500 associations; Increase the number of associations by 10 and calculate the cummulative overlap percentages
  Association_groups<-seq(10,500,by=10)

  Overlaps_all<-data.frame(nrows=length(Association_groups))
  Associations_Overlap<-data.frame(matrix(ncol=6))
  colnames(Associations_Overlap)<-c("Metric1","Metric2","Assoc_Treshold","overlap")
  Associations_Overlap<-Associations_Overlap[complete.cases(Associations_Overlap),]

  for(i in 1:length(Association_groups))
  {
    #Extract the gene drug associations corresponding to different thresholds
    Metric1_Metric2_overlap<-merge(metric1[1:Association_groups[i],],metric2[1:Association_groups[i],],by=c("FEATURE","DRUG_NAME"))
    Overlaps_all<-cbind.data.frame(Metric1=metric1_name,Metric2=metric2_name,Assoc_Threshold=Association_groups[i],overlap=(nrow(Metric1_Metric2_overlap)/Association_groups[i])*100)
    colnames(Associations_Overlap)<-colnames(Overlaps_all)
    Associations_Overlap<-rbind.data.frame(Associations_Overlap,Overlaps_all)

  }
  return(Associations_Overlap)
}
