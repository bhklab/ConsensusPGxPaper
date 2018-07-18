#' Simulate overlap coincidence between the different response metrics
#'
#' Estimate the probability of finding overlaps between 2 different metrics after random sampling of a set of gene drug associations;
#' Different association thresholds are used 
#'
#' @param ANOVA_res ANOVA results of one of the metrics
#' @return A data frame with the number of associations and the random overlap percentages
#' @export
#'


Random_coincidence<-function(ANOVA_res){

total_genes = nrow(ANOVA_res)
associations = seq(10,500,by=10)

Metric_overlap_coincidence_all<-c()

#Expected value of the number of genes that overlap between the 2 sets (n1 and n2) drawn randomly from a total of N genes can be calculated by (n1*n2)/N, as it follows a hypergeometric distribution 
for(i in 1:length(associations)){

#Calculate the random percentage overlap
overlap_coincidence<-(((associations[i]*associations[i])/total_genes)/associations[i])*100
Metric_overlap_coincidence_all<-c(Metric_overlap_coincidence_all,overlap_coincidence)
}
Random_overlap_percentages<-cbind.data.frame(Assoc_Threshold=associations,overlap=Metric_overlap_coincidence_all)
return(Random_overlap_percentages)
}
