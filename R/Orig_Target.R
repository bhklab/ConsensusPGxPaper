#' Primary targets of GDSC drugs
#'
#' Target information downloaded from GDSC \url{https://www.cancerrxgene.org/downloads} and ChEMBL 23 \url{https://www.ebi.ac.uk/chembl/}
#' Only those bioactivities with target confidence of 7 and above were retained; Further the bioactivities were filtered based on the following criteria:
#' Kd/Ki/IC50 < 100nM; Inhibition % > 80% and Residual activities < 20%
#'
#' A dataframe with drugs and the corresponding primary targets
#'
#' @format A data frame with 1711 rows and 2 variables:
#' \describe{
#'   \item{Drug}{Name of the drug}
#'   \item{Target}{Primary target of the drug}
#'   }
#' @source \url{https://www.cancerrxgene.org/downloads}
#' @source \url{https://www.ebi.ac.uk/chembl/}
"Orig_Target"
