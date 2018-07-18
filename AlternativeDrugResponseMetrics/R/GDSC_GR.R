#' Calculate GR for GDSC
#'
#' @param GDSC_finalCellLines Set of cell lines which have passed the quality checks and for which GR metrics is to be calculated
#' @param IGDSC_ref  Reference values of cell lines calculated from relative viabilities of Gemcitabine (GDSC data)
#' @param GR_GCSI_ref Reference values of cell lines calculated from GRs of Gemcitabine (GCSI data)
#' @param GDSC_via_complete Relative viabilities of GDSC cell lines
#' @return A data frame with cellLine, Drug, Concentration and GR values for final set of GDSC cell lines
#' @export
GDSC_GR<-function(GDSC_finalCellLines,IGDSC_ref,GR_GCSI_ref,GDSC_via_complete,writeFile=T,outFile="Calculated_Values/AllGR_GDSCdata.txt")
{
	#Calculate X0(GDSC)
	IGDSC_ref_overlap<-IGDSC_ref[IGDSC_ref$CellLine %in% GDSC_finalCellLines,]
	colnames(IGDSC_ref_overlap)<-c("CellLine","ICRef")
	GR_ref_overlap<-GR_GCSI_ref[GR_GCSI_ref$CellLine %in% GDSC_finalCellLines,]
	colnames(GR_ref_overlap)<-c("CellLine","GRref")
	GR_IGDSC_merge<-merge(IGDSC_ref_overlap,GR_ref_overlap,by="CellLine")
	X0_GDSC<-apply(GR_IGDSC_merge[,c('ICRef','GRref')],1, function(x) x[1]^(1/(1-log2(x[2]+1))))
	X0_GDSC_cellLine<-cbind.data.frame(CellLine=GR_IGDSC_merge$CellLine,X0_GDSC)


	#Calculate GR for GDSC
	GDSC_CellLines_list<-unique(as.vector(X0_GDSC_cellLine$CellLine))
	All_GR<-c()

	for(i in 1:length(GDSC_CellLines_list))
	{
		Cell<-GDSC_CellLines_list[i]
		X0_GDSC_specific_CellLine<-X0_GDSC_cellLine[X0_GDSC_cellLine$CellLine==Cell,"X0_GDSC"]
		GDSC_extract_specific_CellLine<-GDSC_via_complete[GDSC_via_complete$CellLine==Cell,]

		#Set the relative viabilities with negative values to 0
		Rel_via_Neg_Zero<-ifelse(GDSC_extract_specific_CellLine$Rel_via<0,0,GDSC_extract_specific_CellLine$Rel_via)
		GR_GDSC_final<-sapply( Rel_via_Neg_Zero, function(x) ((2^(log2(x/X0_GDSC_specific_CellLine)/log2(1/X0_GDSC_specific_CellLine)))-1))
		AllData_merged<-cbind.data.frame(CellLine=Cell,Drug= GDSC_extract_specific_CellLine$Drug,Conc=GDSC_extract_specific_CellLine$Conc,Via=GDSC_extract_specific_CellLine$Rel_via,GR=GR_GDSC_final)
		All_GR<-rbind(All_GR,AllData_merged)
	}

	if(writeFile) write.table(All_GR,file=outFile,sep="\t",quote=FALSE,row.names=FALSE)

	return(All_GR)
}


