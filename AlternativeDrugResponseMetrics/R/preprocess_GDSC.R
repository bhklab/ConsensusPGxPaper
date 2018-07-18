#' Preprocess GDSc data
#'
#' @param RawData read \code{\link{AllDrugs_CellLines_GDSC}}  for data format
#' @param writeFile (default TRUE) write results into file
#' @param outFile destination of results
#' @return data.frame with CellLine, Drug and corresponding Concentration and cell viability values
#' @export

preprocess_GDSC = function(RawData,writeFile=T,outFile="Calculated_Values/AllCellLines_Drugs_Conc_Via.txt"){

	Conc = calculate_concentrations(RawData)
	CellVia = calculate_relativeViability(RawData)

	prepData = cbind(Conc,Rel_via=CellVia)
	if(writeFile){
	  write.table(prepData,file=outFile,sep="\t",row.names=FALSE,quote=FALSE,col.names=T)
	}
	return(invisible(prepData))
}


#' Calculate_concentrations
#'
#' Calculate concentrations based on fold dilution and Max. conc.
#'
#' @param RawData read  \code{\link{AllDrugs_CellLines_GDSC}} for data format.
#' @param writeFile (default FALSE) write results into file
#' @param outFile destination of results
#' @return data.frame with CellLine, Drug and corresponding Concentration values.
#' @export


calculate_concentrations<-function(RawData,writeFile=F,outFile="AllCellLines_Drugs_Conc.txt")
{

	if(!all(c("MAX_CONC","FOLD_DILUTION","CellLine","Drug") %in% colnames(RawData))){
		stop(paste0("Required column names are: ", paste0(c("MAX_CONC","FOLD_DILUTION","CellLine","Drug"), collapse = ", ")))
	}

	RawData$id = 1:nrow(RawData)

	Conc_data = plyr::ddply(.data = RawData, .variables = "id", .fun = function(r){
		All_Conc = r$MAX_CONC *(r$FOLD_DILUTION ^ seq(0,-8))
		data.frame(CellLine=r$CellLine, Drug=r$Drug,Conc=All_Conc)
	},.progress = plyr::progress_text())

	Conc_data$id = NULL

	if(writeFile){
		write.table(Conc_data,file=outFile,sep="\t",row.names=FALSE,quote=FALSE,col.names=T)
	}
	return(invisible(Conc_data))
}



#' @title Calculate Relative viabilities for GDSC data
#'
#' @description  Calculate relative viabilities (Relative Viability = (Raw intensity - Positive Control) / (Negative Control - Positive Control))
#'
#' @param RawData read \code{\link{AllDrugs_CellLines_GDSC}} for data format.
#' @param Drug_Conc data.frame with Drug, CellLine, and Concetration columns.
#' @param writeFile (default FALSE) write results into file
#' @param outFile destination of results
#' @return vector with relative viabilities
#' @export
calculate_relativeViability<-function(RawData)
{
	Rel_via_all<-c()

	#Calculate relative viabilities
	Neg_control<-RawData[grep("^control",colnames(RawData),value=TRUE)]
	Neg_control_avg<-apply(as.matrix(Neg_control),1,mean,na.rm=TRUE)
	Pos_control<-RawData[grep("^blank",colnames(RawData),value=TRUE)]
	Pos_control_avg<-apply(as.matrix(Pos_control),1,mean,na.rm=TRUE)

	via<-RawData[grep("^raw",colnames(RawData),value=TRUE)]
	Rel_via<-plyr::adply(via, .margin=2, function(x){
	  t((x-Pos_control_avg)/(Neg_control_avg-Pos_control_avg))
	  },.progress = plyr::progress_text())
	Rel_via = t(Rel_via[,-1])
	Rel_via_all<-as.vector(t(Rel_via))

	return(Rel_via_all)
}


