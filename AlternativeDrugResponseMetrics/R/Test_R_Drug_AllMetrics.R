#' Test set predictions of drug response models
#'
#' Test set predicts correspond to the 5 different test sets generated by random splitting; Each test set constitutes 20% of the data
#'
#' A dataframe with the predictions from models based on GRAOC, AUC, GR50 and IC50
#'
#' @format A data frame with 5 rows and 4 variables:
#' \describe{
#'   \item{GRAOC/AUC/GR50/IC50}{Pearson correlation between the observed and predicted values}
#'   }
#' @source \url{https://www.cancerrxgene.org/downloads}
"Test_R_Drug_AllMetrics"