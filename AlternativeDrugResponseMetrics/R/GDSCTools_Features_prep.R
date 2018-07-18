#' Generate response matrix in GDSC format
#'
#' @param metric IC50/AUC/Emax/GR50/GRAOC/GRmax/NDR50/NDRAOC/PI50/PIAOC
#' @return A dataframe with drugs as columns and cell lines as rows; Drug response metrics values are filled up in matrices
#' @export
#'
GDSCTools_Features_prep<-function(GDSC_IC50_GR50,COSMIC_ID,metric,dataset)
{
  #Map cosmic ids to cell lines
  GDSC_IC50_GR50_CosmicIDs_merged<-merge(COSMIC_ID,GDSC_IC50_GR50,by="CellLine")

  #Generate matrices of different metrics
  GDSC_IC50_GR50_extract<-GDSC_IC50_GR50_CosmicIDs_merged[,which(colnames(GDSC_IC50_GR50_CosmicIDs_merged) %in% c("COSMIC_ID","Drug",metric))]
  DR_matrix<-dcast(GDSC_IC50_GR50_extract,COSMIC_ID~Drug)

  #Change column names to match with GDSC tools format
  #Drugs_unique<-unique(as.vector(GDSC_IC50_GR50_extract$Drug))
  colnames_matrix<-paste("Drug_",colnames(DR_matrix)[2:ncol(DR_matrix)],"_IC50",sep="")
  colnames(DR_matrix)<-c("COSMIC_ID",colnames_matrix)

  write.table(DR_matrix,file=paste("ANOVA/GDSC_",metric,dataset,"_matrix.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
}
