#' Calculate NDR, PI and GR values based on the approach described in Gupta et al., 2018
#'
#'Negative control, positive control and viabilities at starting (0h) and end points (72h) of the assay are needed for NDR calculations. This information is not available in GDSC.
#'Therefore, the End points (Measurements after 72h) and Doubling times (GCSI) are used to estimate the starting point values values (Negative control,NC(t)=NC0*(2^t/DT); NC0=NCt/(2^t/DT))
#'
#'
#' @param GDSC_rawData Raw data file with information about fold dilutions, maximum concentration, viabilities, positive and negative control
#' @param Cell_division_times Cell Division times extracted from GCSI data
#' @return A dataframe with the cell line, drug, concentration, relative viability, NDR, PI and GR metrics
#' @export
#'
#'
NDR_PI_GRmetrics<-function(GDSC_rawData,CellDivisionTimes,writeFile=T,outFile="Calculated_Values/AllCellLines_Drugs_RelVia_NDR_GR_PI.txt")
{
#Compute Concentrations and Relative viability
GDSC_SelectiveData_Conc<-calculate_concentrations(GDSC_rawData)

Cell_division_times<-CellDivisionTimes[which(CellDivisionTimes$CellLine %in% GDSC_rawData$CellLine),]
Cell_Drug_allMetrics_combined<-c()

Cell_Drug_allMetrics_headers<-cbind.data.frame("CellLine","Drug","Conc","Via","NDR","PI","GR")
if(writeFile) write.table(Cell_Drug_allMetrics_headers,file=outFile,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)

    result<-tryCatch({
    for(i in 1:nrow(GDSC_rawData))
    {
     GDSC_SelectiveData<-GDSC_rawData[i,]
     DoublingTime_Cell<-Cell_division_times[Cell_division_times$CellLine %in% GDSC_rawData$CellLine[i],"Doubling_Time"]
     Conc<-GDSC_SelectiveData_Conc[which((GDSC_SelectiveData_Conc$CellLine %in% GDSC_SelectiveData$CellLine) & (GDSC_SelectiveData_Conc$Drug %in% GDSC_SelectiveData$Drug ) ),"Conc"]

     #Assay time
     t=72

      #Calculate average positive and negative controls
      Neg_control<-GDSC_SelectiveData[grep("^control",colnames(GDSC_SelectiveData),value=TRUE)]
      Neg_control_avg<-mean(as.numeric(Neg_control),na.rm=TRUE)
      Pos_control<-GDSC_SelectiveData[grep("^blank",colnames(GDSC_SelectiveData),value=TRUE)]
      Pos_control_avg<-mean(as.numeric(Pos_control),na.rm=TRUE)

      N72<-Neg_control_avg
      P72<-Pos_control_avg
      N0=N72/(2^(t/DoublingTime_Cell))


      #Effect of drugs after 72h
      via_drugs_72<-GDSC_SelectiveData[grep("^raw",colnames(GDSC_SelectiveData),value=TRUE)]

      #At the beginning, number of cells should be the same for positive control, negative control and treated wells. It's ideal to use the values calculated for negative control (N0), as these are untreated measurements.
      Normalized_NC<-N72/N0
      Normalized_PC<-P72/N0
      Normalized_rel_via<-via_drugs_72/N0

      NDR_fold_change_nr<-sapply(Normalized_rel_via, function(x) 2^(log2(x)/log2(Normalized_PC)))
      NDR_fold_change_dr<-2^(log2(Normalized_NC)/log2(Normalized_PC))
      NDR<-sapply(NDR_fold_change_nr, function(x) max(-1,((1-x)/(1-NDR_fold_change_dr))))

      #Compute relative viability
      Rel_via<-sapply(via_drugs_72, function(x) (x-P72)/(N72-P72))

      #Recompute GR based on the formula in Gupta etal., 2018
      GR<-sapply(Normalized_rel_via, function(x) (2^(log2(x)/log2(Normalized_NC)))-1)

      #Calculate Percent Inhibition
      PI<-sapply(via_drugs_72, function(x) ((N72-x)/(N72-P72)))

      Cell_Drug_allMetrics<-cbind.data.frame(CellLine=GDSC_SelectiveData$CellLine,Drug=GDSC_SelectiveData$Drug,Conc=Conc,Rel_via=Rel_via,NDR=NDR,PI=PI,GR=GR)
      if(writeFile) write.table(Cell_Drug_allMetrics,file=outFile,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)
      Cell_Drug_allMetrics_combined<-rbind.data.frame(Cell_Drug_allMetrics_combined,Cell_Drug_allMetrics)
    }
  }, error = function(err) {
    print(paste("MY_ERROR:  ",err))
  } )

  return(Cell_Drug_allMetrics_combined)

}
