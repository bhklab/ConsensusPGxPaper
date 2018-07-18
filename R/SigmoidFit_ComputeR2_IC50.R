#' Fit Sigmoid curve; Compute R2 & IC50
#'
#' @param Gemcitabine_conc_filt  A data frame with concentrations and Relative viabilities (Concentration ranges: < 32Mm)
#' @return A data frame with cell lines together with R2 and IC50 extracted from curve fitting
#' @export
SigmoidFit_ComputeR2_IC50<-function(Gemcitabine_conc_filt, writeFile=T, outFile="Calculated_Values/Gemcitabine_GDSC_FittedSigmoid_R2_IC50.txt")
{

	Cell_headers <- cbind.data.frame("CellLine","R2","IC50")
	if(writeFile)  write.table(Cell_headers,file=outFile,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)

	IC50data = plyr::ddply(Gemcitabine_conc_filt,.variables ="CellLine", function(CellLine_data){
		result<-tryCatch({
			Sigmoid_curve<-drm(
				Rel_via~Conc,data=CellLine_data,
				fct=drc::LL.3u())
			RSS2 = sum(stats::residuals(Sigmoid_curve)^2, na.rm = TRUE)
			RSS1=sum((CellLine_data$Rel_via - mean(CellLine_data$Rel_via, na.rm = TRUE))^2,na.rm = TRUE)
			R2<-1-(RSS2/RSS1)
			IC50<-Sigmoid_curve$coefficients[3]
			Cell_IC50<-cbind.data.frame(CellLine=as.vector(CellLine_data[1,"CellLine"]),R2=R2,IC50=IC50)
			if(writeFile)  write.table(Cell_IC50,file=outFile,sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)
			return(Cell_IC50)
		}, error = function(err) {
			print(paste("MY_ERROR:  ",err))
			print(paste("calculation for ",CellLine_data[1,"CellLine"], " failed. "))
			return(data.frame())
		} )

	},.progress = plyr::progress_text() )
	return(IC50data)
}
