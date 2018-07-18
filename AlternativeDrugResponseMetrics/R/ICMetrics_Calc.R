#' Calculate the conventional metrics
#'
#' Fit sigmoid curves to compute the conventional metrics (IC50,EC50,Emax,ICinf,IChill,AAC,ICR2) for GDSC
#' Significance of curve fits assessed based on the conditions specified in Marc Hafner's paper
#' Concentrations are capped to the highest concentration value, when the GEC50s are 10 times more than the highest concentration; Same also applies, when the p-values resulting from the significance of fits is >0.5
#'
#'
#' @param All_GDSC_IC_data A data frame with CellLine, Drug, Concentration and Relative Viability of all GDSC data
#' @return A data frame with all the conventional metrics
#' @export
#'
ICMetrics_Calc<-function(All_GDSC_IC_data,writeFile=T,fileName)
{
  GDSC_Cells<-unique(as.vector(All_GDSC_IC_data$CellLine))
  ICmetrics_headers<-cbind.data.frame("CellLine","Drug","IC50","EC50","Emax","IC_inf","IC_HillCoeff","AUC","IC_R2","LN_IC50","AAC")
  if(writeFile) write.table(ICmetrics_headers,file=fileName,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
  Cell_Drug_IC_allMetrics<-c()

  for(i in 1:length(GDSC_Cells))
  {
    CellLine_data<-All_GDSC_IC_data[All_GDSC_IC_data$CellLine==GDSC_Cells[i],]
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
      lower=c(.1, 0, min(c)*1e-2)
      upper = c(5, 1, max(c) * 100)

      controls = drc::drmc()
      controls$relTol = 1e-06
      controls$errorm = FALSE
      controls$noMessage = TRUE
      controls$rmNA = TRUE

      result<-tryCatch({
        Sigmoid_curve<-drm(
          Via~Conc,data=CellLine_drugs_data_conc_filter,
          fct=drc::LL.3u(names = c("h_IC",
                                   "ICinf", "EC50")),start=priors,lowerl = lower,
          upperl = upper, control = controls, na.action = na.omit)

        #Perform f test for estimating the significance of fits
        RSS2 = sum(stats::residuals(Sigmoid_curve)^2, na.rm = TRUE)
        RSS1=sum(( CellLine_drugs_data_conc_filter$Via - mean(CellLine_drugs_data_conc_filter$Via, na.rm = TRUE))^2,na.rm = TRUE)

        Npara=3
        Npara_flat=1
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(CellLine_drugs_data_conc_filter$Via)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = pf(f_value, df1, df2, lower.tail = FALSE)
        pval = f_pval

        #All IC metrics; Code from Hafner's Github used
        IC_R2<-1-(RSS2/RSS1)
        IC_mean=mean(na.omit(CellLine_drugs_data_conc_filter$Via))

        if(pval<=0.05)
        {
          IC50<-Sigmoid_curve$coefficients[3]
          H_IC<-Sigmoid_curve$coefficients[1]
          IC_inf<-Sigmoid_curve$coefficients[2]
          EC50<-ifelse(EC50>(Highest_conc*10),Highest_conc,Sigmoid_curve$coefficients[3])
          if((IC_inf>0.5) | (EC50==Inf))
          {
            IC50= Highest_conc
            EC50= Highest_conc
          }else{
            IC50<-EC50 * ((1 - IC_inf)/(0.5 -
                                          IC_inf) - 1)^(1/H_IC)
          }
        }else{
          EC50= Highest_conc
          H_IC=0.01
          IC50= Highest_conc
          IC_inf=IC_mean
        }

        ICavg<-CellLine_drugs_data_ordered$Via
        concs<-CellLine_drugs_data_ordered$Conc
        AUC = sum(((ICavg[1:(length(ICavg)-1)]+ICavg[2:length(ICavg)])/2)*diff(log10(concs), lag = 1), na.rm = TRUE)/(log10(concs[length(concs)]) - log10(concs[1]))
        AAC= sum((1 - (ICavg[1:(length(ICavg) - 1)] + ICavg[2:length(ICavg)])/2) * diff(log10(concs), lag = 1), na.rm = TRUE)/(log10(concs[length(concs)]) - log10(concs[1]))

        ICmax<-min(CellLine_drugs_data_ordered$Via[nrow(CellLine_drugs_data_ordered)],CellLine_drugs_data_ordered$Via[nrow(CellLine_drugs_data_ordered)-1])


        IC50_final<-ifelse(IC50>(Highest_conc*10),Highest_conc,IC50)
        Cell_allICmetrics<-cbind.data.frame(unique(as.vector(CellLine_drugs_data_conc_filter$CellLine)),unique(as.vector(CellLine_drugs_data_conc_filter$Drug)),IC50_final,EC50,ICmax,IC_inf,H_IC,AUC,IC_R2,log(IC50_final),AAC)
        Cell_Drug_IC_allMetrics<-rbind(Cell_Drug_IC_allMetrics,Cell_allICmetrics)

        if(writeFile) write.table(Cell_allICmetrics,file=fileName,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)
      }, error = function(err) {
        print(paste("MY_ERROR:  ",err))
        print(paste("Calculations for CellLine",GDSC_Cells[i], "and Drug", Drugs_specific_cell[j],"failed",sep=""))

      } )
    }
  }
  return(Cell_Drug_IC_allMetrics)
}
