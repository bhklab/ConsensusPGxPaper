#' Calculate variance of GR and Rel_via for each cell line
#'
#' Only responses corresponding to concentrations<32Mm considered
#'
#' @param Gemcitabine_conc_filt Gemcitabine responses (GR/Rel_via)
#' @param metric GR/Rel_via
#' @export
#'
CellLines_var<-function(Gemcitabine_conc_filt,metric=c("GR","Rel_via"))
{
  metric = match.arg(metric)
  Gemcitabine_var = plyr::ddply(Gemcitabine_conc_filt,.variables ="CellLine", function(Gemcitabine_variance){
	Response_var<-var(Gemcitabine_variance[,metric])
	Var_CellLine<-cbind.data.frame(CellLine=unique(as.vector(Gemcitabine_variance[,"CellLine"])),Response_var)
	},.progress = plyr::progress_text())
	return(Gemcitabine_var)
}
