#' Calculate all the associated drug response metrices
#'
#' Fit sigmoid curves and compute all associated GR metrics (GR50,GEC50,GRmax,GRinf,GRhill,GRAOC,GRR2) and other metrics NDR / PI for GDSC
#' Significance of curve fits assessed based on the conditions specified in Marc Hafner's paper
#' Concentrations are capped to the highest concentration value, when the GEC50s are 10 times more than the highest concentration; Same also applies, when the p-values resulting from the significance of fits is >0.5
#'
#'
#' @param All_GDSC_GR_data A data frame with CellLine, Drug, Concentration, GR/NDR/PI and Relative Viability of all GDSC data
#' @return A data frame with all the GR/NDR/PI associated metrics
#' @export
#'

AllMetrics_Calc<-function(All_GDSC_GR_data,metric,Metrics_headers,writeFile=T,fileName)
{
  metric<-noquote(metric)
  GDSC_Cells<-unique(as.vector(All_GDSC_GR_data$CellLine))

  if(writeFile) write.table(Metrics_headers,file=fileName,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)
  Cell_Drug_allMetrics<-c()
  for(i in 1:length(GDSC_Cells))
  {
    CellLine_data<-All_GDSC_GR_data[All_GDSC_GR_data$CellLine==GDSC_Cells[i],]
    Drugs_specific_cell<-unique(as.vector(CellLine_data$Drug))
    for(j in 1:length(Drugs_specific_cell))
    {
      CellLine_drugs_data<-CellLine_data[CellLine_data$Drug==Drugs_specific_cell[j],]
      CellLine_drugs_data_conc_filter<-CellLine_drugs_data
      CellLine_drugs_data_ordered<-CellLine_drugs_data_conc_filter[order(CellLine_drugs_data_conc_filter$Conc),]
      Highest_conc<-max(CellLine_drugs_data_ordered$Conc)

      #All these parameters are set based on the limits specified in Hafner's paper
      c = unique(CellLine_drugs_data_conc_filter$Conc)
      priors=c(2, 0.1, stats::median(c))
      lower = c(0.1, -1, min(c) * 0.01)
      upper = c(5, 1, max(c) * 100)

      controls = drc::drmc()
      controls$relTol = 1e-06
      controls$errorm = FALSE
      controls$noMessage = TRUE
      controls$rmNA = TRUE

      result<-tryCatch({
        Sigmoid_curve<-drm(
          CellLine_drugs_data_conc_filter[,metric]~Conc,data=CellLine_drugs_data_conc_filter,
          fct=drc::LL.3u(names = c("h_Metric",
                                   "Metric_inf", "Metric_EC50")),start=priors,lowerl = lower,
          upperl = upper, control = controls, na.action = na.omit)

        #Perform f test for estimating the significance of fits
        RSS2 = sum(stats::residuals(Sigmoid_curve)^2, na.rm = TRUE)
        RSS1=sum(( CellLine_drugs_data_conc_filter[,metric] - mean(CellLine_drugs_data_conc_filter[,metric], na.rm = TRUE))^2,na.rm = TRUE)


        Npara=3
        Npara_flat=1
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(CellLine_drugs_data_conc_filter[,metric])) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = pf(f_value, df1, df2, lower.tail = FALSE)
        pval = f_pval

        #All Alternative metrics; Code from Hafner's Github used
        Metric_R2<-1-(RSS2/RSS1)
        Metric_mean=mean(na.omit(CellLine_drugs_data_conc_filter[,metric]))

        if(pval<=0.05)
        {
          Metric_EC50<-Sigmoid_curve$coefficients[3]
          H_Metric<-Sigmoid_curve$coefficients[1]
          Metric_inf<-Sigmoid_curve$coefficients[2]
          Metric_EC50<-ifelse(Metric_EC50>(Highest_conc*10),Highest_conc,Sigmoid_curve$coefficients[3])
          if((Metric_inf>0.5) | (Metric_EC50==Inf))
          {
            Metric_50= Highest_conc
            Metric_EC50= Highest_conc
          }else{
           Metric_50<-Metric_EC50 * ((1 - Metric_inf)/(0.5 -
                                           Metric_inf) - 1)^(1/H_Metric)
          }
        }else{
          Metric_EC50= Highest_conc
          H_Metric=0.01
          Metric_50= Highest_conc
          Metric_inf=Metric_mean
        }

        Metric_avg<-CellLine_drugs_data_ordered[,metric]
        concs<-CellLine_drugs_data_ordered$Conc
        Metric_AOC<- sum((1 - (Metric_avg[1:(length(Metric_avg) - 1)] + Metric_avg[2:length(Metric_avg)])/2) * diff(log10(concs), lag = 1), na.rm = TRUE)/(log10(concs[length(concs)]) - log10(concs[1]))
        Metric_max<-min(CellLine_drugs_data_ordered[nrow(CellLine_drugs_data_ordered),metric],CellLine_drugs_data_ordered[(nrow(CellLine_drugs_data_ordered)-1),metric])


        Metric_50_final<-ifelse(Metric_50>(Highest_conc*10),Highest_conc,Metric_50)
        Cell_allMetrics<-cbind.data.frame(unique(as.vector(CellLine_drugs_data_conc_filter$CellLine)),unique(as.vector(CellLine_drugs_data_conc_filter$Drug)),Metric_50_final,Metric_EC50,Metric_max,Metric_inf,H_Metric,Metric_AOC,Metric_R2,log(Metric_50_final))
        Cell_Drug_allMetrics<-rbind(Cell_Drug_allMetrics,Cell_allMetrics)

        if(writeFile) write.table(Cell_allMetrics,file=fileName,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)
      }, error = function(err) {
        print(paste("MY_ERROR:  ",err))
        print(paste("Calculations for CellLine",GDSC_Cells[i], "and Drug", Drugs_specific_cell[j],"failed"))

      } )
    }
  }
  return(Cell_Drug_allMetrics)
}
