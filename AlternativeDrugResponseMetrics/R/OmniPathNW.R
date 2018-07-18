#' Omnipath network
#'
#' Signed and weighted protein-protein interaction networks downloaded from \url{http://pypath.omnipathdb.org}
#'
#' A dataframe with information about the proteins that are connected, directionality and strength of the connection
#'
#' @format A data frame with 8858 rows and 4 variables:
#' \describe{
#'   \item{Node1}{Protein1}
#'   \item{Sign}{Directionality of the protein-protein interaction}
#'   \item{Node2}{Protein2}
#'   \item{Weight}{Strength of the interaction}
#'   }
#' @source \url{http://pypath.omnipathdb.org}
"OmniPathNW"