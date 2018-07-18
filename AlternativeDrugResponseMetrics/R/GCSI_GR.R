#' Calculate GR values for Gemcitabine (Reference compound);
#'
#' Division times and viabilities from GCSI used
#' @param Gemcitabine_GCSI_Via_Doses  Dose response data of Gemcitabine
#' @param CellDivisionRates  Division times of the cell lines
#' @param writeFile (default TRUE) Save results to a file
#' @param outFile file destination for results
#' @return data.frame with CellLine, Conc  and GR columns
#' @export
#'
GCSI_GR<-function(Gemcitabine_GCSI_Via_Doses, CellDivisionRates, Assay_Time, writeFile=T, outFile="Calculated_Values/Gemcitabine_GCSI_GR_allCellLines.txt")
{

	#Extract raw viabilities and Doses
	GCSI_data_viability<-grep("^Median",colnames(Gemcitabine_GCSI_Via_Doses),value=TRUE)
	if(length(GCSI_data_viability)==0) stop("Median Raw Viability data not found among the columns...")
	GCSI_via_all_cellLines<-Gemcitabine_GCSI_Via_Doses[,GCSI_data_viability]

	GCSI_data_Conc<-grep("^Dose",colnames(Gemcitabine_GCSI_Via_Doses),value=TRUE)
	if(length(GCSI_data_Conc)==0) stop("Dose data not found among the columns...")
	GCSI_Conc_all_cellLines<-Gemcitabine_GCSI_Via_Doses[,GCSI_data_Conc]

	AllCellLines_GR_merged<-c()

	GCSI_unique_cellLines<-unique(as.vector(Gemcitabine_GCSI_Via_Doses$CellLine))
	for(i in 1:nrow(GCSI_via_all_cellLines))
	{
		CellLines_GR_merged<-c()
		Cell<-as.vector(CellDivisionRates$CellLine[i])
		Td<-CellDivisionRates$Doubling_Time[i]
		T_Td<-Assay_Time/Td
		Conc<-t(data.frame(GCSI_Conc_all_cellLines[i,]))
		Relative_via<-GCSI_via_all_cellLines[i,]

		GR_allconc<-sapply(Relative_via, function(x) ((2^(1+(log2(x)/T_Td)))-1))
		CellLines_GR_merged<-cbind.data.frame(Cell,Conc,as.data.frame(GR_allconc))
		colnames(CellLines_GR_merged)<-c('CellLine','Conc','GR')
		AllCellLines_GR_merged<-rbind(AllCellLines_GR_merged,CellLines_GR_merged)
	}
	if(writeFile) write.table(AllCellLines_GR_merged,file=outFile,sep="\t",row.names=FALSE,quote=FALSE)

	return(AllCellLines_GR_merged)
}
