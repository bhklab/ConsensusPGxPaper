#' Calculate interquartile means for Gemcitabine data across different concentrations
#'
#' Interquartile means correspond to the reference values used for calculating GRs for other cell lines and drugs in GDSC
#'
#' @param Gemcitabine_GR_filtered  Filtered Gemcitabine response data between concentration 0.1 and 32Mm
#' @param type  GCSI_GRRef or GDSC_ICRef (Specify whether the interquartile means are to be calculated for the Gemcitabine data from GCSI or GDSC)
#' @param metric GR or Rel_via (Specify whether GR or Relative Viability values are to be used for interquartile men calculations)
#' @return A dataframe with cell lines and interquartile means of GRs or Relative viabilities
#' @export
InterquartileMean <- function(Gemcitabine_GR_filtered,type,metric=c("Rel_via","GR"),writeFile=T,outFile=paste("Calculated_Values/Gemcitabine_",type,"_allCellLines.txt",sep=""))
{
	metric = match.arg(metric)

	Gemcitabine_Ref = plyr::ddply(Gemcitabine_GR_filtered,.variables ="CellLine", function(Gemcitabine_data){


		GR_sorted<-sort(as.vector(Gemcitabine_data[,metric]))

		#If there is only one data point satisfying the concentration conditions, GRref would be be the same as the relative viability; Otherwise, it's calculated as shown below.
		if(nrow(Gemcitabine_data)>1)
		{
			#If the total number of datapoints satisfying the concentration conditions is even, interquartile mean is the arithmretic mean of second and third quartiles.
			if((length(GR_sorted) %% 4)==0 )
			{
				Quartiles_entry<-floor(length(GR_sorted)/4)
				GR_sorted_2nd_3rd_quartile<-GR_sorted[(Quartiles_entry+1):(length(GR_sorted)-(Quartiles_entry))]
				GR_sorted_IQR<-mean(GR_sorted_2nd_3rd_quartile)
			}else{
				#If the total number of datapoints satisfying the concentration conditions is odd, interquartile mean is the weighted average of the quartiles and interquartiles
				Quartiles_entry<-floor(length(GR_sorted)/4)
				First_fourth_quartile<-(length(GR_sorted)/4)*2
				Interquartiles<-length(GR_sorted)-(Quartiles_entry*2)
				Interquartiles_fract<-(First_fourth_quartile-(Interquartiles-2))/2
				Index<-which(GR_sorted==median(GR_sorted))
				WholeNums<-seq((Quartiles_entry+2),(length(GR_sorted)-(Quartiles_entry+1)),by=1)
				Sums_wholenums<-sum(GR_sorted[WholeNums])
				FractNums<-range(WholeNums[1]-1,WholeNums[length(WholeNums)]+1)
				Sums_FractNums<-sum(GR_sorted[FractNums])
				GR_sorted_IQR<-sum(Sums_wholenums,(Sums_FractNums* Interquartiles_fract))/ First_fourth_quartile
			}
		}else{
			GR_sorted_IQR=GR_sorted
		}
		GRRef<-cbind.data.frame(CellLine=Gemcitabine_data[1,"CellLine"],Ref=GR_sorted_IQR)
	},.progress = plyr::progress_text() )
	if(writeFile) write.table(Gemcitabine_Ref,file=outFile,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)
	return(Gemcitabine_Ref)
}
