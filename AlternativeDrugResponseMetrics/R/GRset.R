#' ANOVA results corresponding to all metrics in GR set
#'
#' A dataframe with the drug gene associations, Effect size and FDR thresholds
#' For more details regarding the ANOVA results, refer \url {http://gdsctools.readthedocs.io/en/master/}
#'
#' @name GRset
#' @rdname GRset
#' @aliases results_GRAOC_GRset results_GDSC_AUC_GRset results_RecomputedAUC_GRset results_GRmax_GRset results_Emax_GRset results_LNGR50_GRset results_GDSC_LNIC50_GRset results_RecomputedLNIC50_GRset
#' @format A data frame with 49438 /49533 rows(only in case of GDSC IC50s and GDSC AUCs) and 21 variables:
#' \describe{
#'   \item{ASOC_ID}{Association number}
#'   \item{Feature}{Mutation / CNV status of the gene}
#'   \item{DRUG_ID}{IDs of the drugs}
#'   \item{DRUG_NAME}{Name of the drug}
#'   \item{DRUG_TARGET}{Drug target information included, if available in GDSC}
#'   \item{N_FEATURE_neg}{Number of cell lines that do not have a specific mutation or CNV}
#'   \item{N_FEATURE_pos}{Number of cell lines that have a specific mutation or CNV}
#'   \item{FEATURE_pos_logIC50_MEAN}{Mean logIC50s of the cell lines that have a specific gene feature}
#'   \item{FEATURE_neg_logIC50_MEAN}{Mean logIC50s of the cell lines that do not have a specific gene feature}
#'   \item{FEATURE_delta_MEAN_IC50}{Differences in means of logIC50s between the cell lines that have and do not have a specific gene feature}
#'   \item{FEATURE_IC50_effect_size}{Effect size computed based on Cohen's metric}
#'   \item{FEATURE_neg_Glass_delta}{(Mean logIC50s / SD logIC50s) for negative sets}
#'   \item{FEATURE_pos_Glass_delta}{(Mean logIC50s / SD logIC50s) for positive sets}
#'   \item{FEATURE_neg_IC50_sd}{Standard deviation of the logIC50s for negative sets}
#'   \item{FEATURE_pos_IC50_sd}{Standard deviation of the logIC50s for positive sets}
#'   \item{ANOVA_FEATURE_pval}{p value representing the significance of drug-gene association}
#'   \item{ANOVA_TISSUE_pval}{p value representing the significance of drug-gene association and tissue type of the cell lines}
#'   \item{ANOVA_MSI_pval}{p value representing the significance of drug-gene association and Micro Satellite Instability status of the cell lines}
#'   \item{ANOVA_MEDIA_pval}{p value representing the significance of drug-gene association and medium of the cell lines}
#'   \item{ANOVA_FEATURE_FDR}{False Discovery Rate resulting from the correction of p-values in "ANOVA_FEATURE_pval" column}
#'   }
#' @source \url{http://gdsctools.readthedocs.io/en/master/}


"results_GRAOC_GRset"
#' @name GRset
#' @rdname GRset
#'
"results_GDSC_AUC_GRset"
#' @name GRset
#' @rdname GRset
#'
"results_RecomputedAUC_GRset"
#' @name GRset
#' @rdname GRset
#'
"results_GRmax_GRset"
#' @name GRset
#' @rdname GRset
#'
"results_Emax_GRset"
#' @name GRset
#' @rdname GRset
#'
"results_LNGR50_GRset"
#' @name GRset
#' @rdname GRset
#'
"results_GDSC_LNIC50_GRset"
#' @name GRset
#' @rdname GRset
#'
"results_RecomputedLNIC50_GRset"
#' @name GRset
#' @rdname GRset
#'
