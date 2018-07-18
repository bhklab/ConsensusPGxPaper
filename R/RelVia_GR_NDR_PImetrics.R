#' Conventional and alternative response metrics for 339 cell lines and 265 drugs
#'
#' A dataframe with concentrations and different metrics calculated by NDR_PI_GRmetrics function
#'
#' For more details concerning the different metrics, refer \url {https://www.ncbi.nlm.nih.gov/pubmed/27180993} and \url {https://www.biorxiv.org/content/early/2018/02/08/262568}
#'
#' @format A data frame with 681543 rows and 7 variables:
#' \describe{
#'   \item{CellLine}{Name of the cell line}
#'   \item{Drug}{Name of the drug, can be also Drug ID.}
#'   \item{Conc}{ Concentrations calculated by the calculate_concentrations function}
#'   \item{Rel_via}{Relative viabilities}
#'   \item{NDR}{Normalized Drug Response}
#'   \item{PI}{Percentage Inhibition}
#'   \item{GR}{Growth Rates}
#'   }
#'
"RelVia_GR_NDR_PImetrics"
