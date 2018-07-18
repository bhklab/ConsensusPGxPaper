#' All associated conventional and GR metrics for 173 cell lines and 265 drugs
#'
#' All associated conventional metrics:IC50,EC50,Emax,ICHill,IC_inf,IC_R2,AUC
#' All associated alternative metrics:GR50,GEC50,GRmax,GRHill,GR_inf,GR_R2,GRAOC
#' 
#' A dataframe with all associated metrics calculated by AllMetrics_Calc and ICMetrics_Calc function
#'
#'
#' @format A data frame with 38393 rows and 19 variables:
#' \describe{
#'   \item{CellLine}{Name of the cell line}
#'   \item{Drug}{Name of the drug, can be also Drug ID.}
#'   \item{IC50/GR50}{Half maximal inhibitory concentrations}
#'   \item{LN_IC50/LN_GR50}{log Half maximal inhibitory concentrations}
#'   \item{AUC/AAC/GRAOC}{Area Under/Above/Over the Curves}
#'   \item{EC50/GEC50}{Half maximal effective concentrations}
#'   \item{Emax/GRmax}{Maximal response}
#'   \item{IC_inf/GR_inf}{Response at infinite concentrations}
#'   \item{IC_HillCoeff,GR_HillCoeff}{Slope of the dose response curves}
#'   \item{IC_R2,GR_R2}{Goodness of fit of sigmoidal curves}
#'   }
#'
"Conventional_GR_allCalculatedMetrics_merged"