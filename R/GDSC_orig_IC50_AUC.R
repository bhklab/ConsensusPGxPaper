#' Original IC50s and AUCs of GDSC data
#'
#' Drug response information downloaded from GDSC \url{https://www.cancerrxgene.org/downloads}
#' IC50s and AUCs in this correspond to the values generated based on non-linear mixed effects model. For details, refer \url{https://www.ncbi.nlm.nih.gov/pubmed/27180993}
#'
#' A dataframe with cell lines, drugs, LNIC50s and AUCs
#'
#' @format A data frame with 224202 rows and 4 variables:
#' \describe{
#'   \item{CellLine}{Name of the cell line}
#'   \item{Drug}{Name of the drug}
#'   \item{LN_IC50}{Half maximal Inhibitory Concentrations}
#'   \item{AUC}{Area Under the Curve}
#'   }
#' @source \url{https://www.cancerrxgene.org/downloads}
"GDSC_orig_IC50_AUC"
