#' Filter the cell lines based on the criteria described in Hafner's publication
#'
#' @param Gemcitabine_GR_Ref A dataframe with cell lines and Reference values computed from Gemcitabine responses (GCSI data)
#' @param Gemcitabine_GR_var A dataframe with cell lines and variance values computed from Gemcitabine responses (GCSI data)
#' @param Gemcitabine_IC_Ref  A dataframe with cell lines and Reference values computed from Gemcitabine responses (GDSC data)
#' @param Gemcitabine_IC_Var A dataframe with cell lines and variance values computed from Gemcitabine responses (GDSC data)
#' @return A vector with final set of cell lines to be used for GR calculations
#' @export
CellLines_filter<-function(Gemcitabine_GR_Ref,Gemcitabine_GR_var,Gemcitabine_IC_Ref,Gemcitabine_IC_Var,GDSC_R2,GCSI_published_Gemcitabine,writeFile=F,outFile="Calculated_Values/GDSC_final_CellLines_toCalculate_GR.txt")
{
	#Pick CellLines with GRref(GCSI) below 0.7
	Gemcitabine_GCSI_GR_0.7<-Gemcitabine_GR_Ref[which(Gemcitabine_GR_Ref$Ref<0.7),]
	Gemcitabine_GCSI_GR_0.7_uniqueCellLines<-unique(as.vector(Gemcitabine_GCSI_GR_0.7$CellLine))

	#Pick cell lines with GRvar below 0.01
	Gemcitabine_GCSI_GRvar_0.01<-Gemcitabine_GR_var[which(Gemcitabine_GR_var$Response_var<0.01),]
	Gemcitabine_GCSI_GRref_0.7_GRvar_0.01<-Gemcitabine_GCSI_GR_0.7[Gemcitabine_GCSI_GR_0.7$CellLine %in% Gemcitabine_GCSI_GRvar_0.01$CellLine,]

	#Pick CellLines with Iref(GDSC) below 0.85 and Ivar below 0.01
	GDSC_Iref_0.85<-Gemcitabine_IC_Ref[which(Gemcitabine_IC_Ref$Ref<0.85),]
	Gemcitabine_GDSC_Ivar_0.01<-Gemcitabine_IC_Var[which(Gemcitabine_IC_Var$Response_var<0.01),]
	Gemcitabine_GDSC_IRef_0.85_Ivar_0.01<-GDSC_Iref_0.85[GDSC_Iref_0.85$CellLine %in% Gemcitabine_GDSC_Ivar_0.01$CellLine,]

	#Pick cell lines with fitted R2>0.5 (GDSC data)
	GDSC_R2_gt0.5<-GDSC_R2[which(GDSC_R2$R2>0.5),]

	#Pick those cell lines with difference between GEC50 in GCSI (Downloaded from GR50 website) and IC50 in GDSC data < 1.5
	GCSI_GDSC_merge<-merge(GCSI_published_Gemcitabine,GDSC_R2,by="CellLine")
	GEC50_IC50_diff<-abs(GCSI_GDSC_merge$GEC50-GCSI_GDSC_merge$IC50)
	GCSI_GDSC_merge_diff<-cbind(GCSI_GDSC_merge,GEC50_IC50_diff)
	GCSI_GDSC_merge_diff_lt1.5<-GCSI_GDSC_merge_diff[GCSI_GDSC_merge_diff$GEC50_IC50_diff<1.5,]

	#GCSI GDSC overlap
	GCSI_GDSC_overlap<-Gemcitabine_GR_Ref[GCSI_Gemcitabine_published$CellLine %in% Gemcitabine_IC_Ref$CellLine,]

	#GCSI filters
	GEC50diff_CellLines<-GCSI_GDSC_merge_diff_lt1.5[GCSI_GDSC_merge_diff_lt1.5$CellLine %in% GCSI_GDSC_overlap$CellLine,]
	GEC50_Gref_0.7_Gvar0.01<-GEC50diff_CellLines[GEC50diff_CellLines$CellLine %in% Gemcitabine_GCSI_GRref_0.7_GRvar_0.01$CellLine,]

	#GDSC filters
	GDSC_Iref_GCSI_filter<-Gemcitabine_GDSC_IRef_0.85_Ivar_0.01[Gemcitabine_GDSC_IRef_0.85_Ivar_0.01$CellLine %in% GEC50_Gref_0.7_Gvar0.01$CellLine,]
	GDSC_Iref_GCSI_filter_GDSC_0.5<-GDSC_Iref_GCSI_filter[GDSC_Iref_GCSI_filter$CellLine %in% GDSC_R2_gt0.5$CellLine,]


	#Final GDSC cell line list for GR calculations
	GDSC_finalCellLines<-unique(as.vector(GDSC_Iref_GCSI_filter_GDSC_0.5$CellLine))
	if(writeFile) write.table(GDSC_finalCellLines,file=outFile,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

	return(GDSC_finalCellLines)
}
