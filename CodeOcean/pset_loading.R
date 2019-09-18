
## ADD PEARSON SPEARMAN FOR BIOMARKER DISCOVERY
##

source("foo.R")
library(PharmacoGx)
library(Biobase)

path.result <- "../results_Cello"
psets.path <- "~/Documents/Cello_Pubchem_PSets"
#psets.path <- "../Psets"

if("AZ" %in% datasets){
  load(file.path(psets.path, "AZ.RData"), verbose=T)
  AZ.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=AZ, sensitivity.measure=sensitivity_metric["AZ"], summary.stat="median"))
  AZ.drug.sensitivity <- AZ.drug.sensitivity[,apply(!is.na(AZ.drug.sensitivity), 2, sum, na.rm=T) != 0]
  AZ_drug_response_var <- apply(AZ.drug.sensitivity, 2, var, na.rm=T)
  drug.mad.AZ <- apply(AZ.drug.sensitivity, 2, mad, na.rm=TRUE)
  AZ_narrow_effect <- names(which(drug.mad.AZ <= .13))
  AZ_broad_effect <- names(which(drug.mad.AZ > .13))
}
if("GDSC1" %in% datasets){
  load(file.path(psets.path, "GDSCv1.RData"), verbose=T)
  GDSCv1 <- GDSC
  GDSCv1.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSCv1, sensitivity.measure=sensitivity_metric["GDSC1"], summary.stat="median"))/sensitivity_factor["GDSC1"]
}
if("GDSC2" %in% datasets){
  load(file.path(psets.path, "GDSCv2.RData"), verbose=T)
  GDSCv2 <- GDSC
  GDSCv2.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSCv2, sensitivity.measure=sensitivity_metric["GDSC2"], summary.stat="median"))/sensitivity_factor["GDSC2"]
  gdsc_drug_response_var <- apply(GDSCv2.drug.sensitivity, 2, var, na.rm=T)
  drug.mad.gdsc <- apply(GDSCv2.drug.sensitivity, 2, mad, na.rm=TRUE)
  gdsc_narrow_effect <- names(which(drug.mad.gdsc <= .13))
  gdsc_broad_effect <- names(which(drug.mad.gdsc > .13))
}

if("gCSI" %in% datasets){
  load(file.path(psets.path, "gCSI_V1.RData"), verbose=T)
  #gCSI <- gCSI_2018_merged
  gcsi.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=gCSI, sensitivity.measure=sensitivity_metric["gCSI"], summary.stat="median"))/sensitivity_factor["gCSI"]
  gcsi_drug_response_var <- apply(gcsi.drug.sensitivity, 2, var, na.rm=T)
  drug.mad.gcsi <- apply(gcsi.drug.sensitivity, 2, mad, na.rm=TRUE)
  gcsi_narrow_effect <- names(which(drug.mad.gcsi <= .13))
  gcsi_broad_effect <- names(which(drug.mad.gcsi > .13))
}
if("GRAY" %in% datasets){
  #load(file.path(psets.path, "GRAY_kallisto.RData"), verbose=T)
  load(file.path(psets.path, "GRAY2017_updated.RData"), verbose=T)
  #load(file.path(psets.path, "GRAY2013_updated.RData"), verbose=T)
  GRAY <- GRAY2017_updated
  gray.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GRAY, sensitivity.measure=sensitivity_metric["GRAY"], summary.stat="median"))/sensitivity_factor["GRAY"]
}
if("CTRPv2" %in% datasets){
  load(file.path(psets.path, "CTRPv2_updated.RData"), verbose=T)
  ctrpv2.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=CTRPv2, sensitivity.measure=sensitivity_metric["CTRPv2"], summary.stat="median"))/sensitivity_factor["CTRPv2"]
  ctrpv2_drug_response_var <- apply(ctrpv2.drug.sensitivity, 2, var, na.rm=T)
  drug.mad.ctrpv2 <- apply(ctrpv2.drug.sensitivity, 2, mad, na.rm=TRUE)
  ctrpv2_narrow_effect <- names(which(drug.mad.ctrpv2 <= .13))
  ctrpv2_broad_effect <- names(which(drug.mad.ctrpv2 > .13))
}
if("CCLE" %in% datasets){
  load(file.path(psets.path, "CCLE.RData"), verbose=T)
  ccle.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE, sensitivity.measure=sensitivity_metric["CCLE"], summary.stat="median"))/sensitivity_factor["CCLE"]
  ccle_drug_response_var <- apply(ccle.drug.sensitivity, 2, var, na.rm=T)
  drug.mad.ccle <- apply(ccle.drug.sensitivity, 2, mad, na.rm=TRUE)
  ccle_narrow_effect <- names(which(drug.mad.ccle <= .13))
  ccle_broad_effect <- names(which(drug.mad.ccle > .13))
}
#load(file.path(psets.path, "GDSC.RData"), verbose=T)
gdsc.mutation <- exprs(summarizeMolecularProfiles(pSet=GDSC, mDataType="mutation", summary.stat="or"))
gdsc.cnv <- exprs(summarizeMolecularProfiles(pSet=GDSC, mDataType="cnv", summary.stat="median"))
gdsc.fusion <- exprs(summarizeMolecularProfiles(pSet=GDSC, mDataType="fusion", summary.stat="or"))

xx <- featureInfo(GDSC, "fusion")
xx$Symbol <-rownames(featureInfo(GDSC, "fusion"))
featureInfo(GDSC, "fusion") <- xx
#load(file.path(psets.path, "GRAY_kallisto.RData"), verbose=T)

#molecular profiles

chemo <- c("5-FU", "Bortezomib", "Doxorubicin", "etoposide", "Methotrexate", "Docetaxel", "paclitaxel", "Gemcitabine")


hist(GRAY@sensitivity$profiles$auc_recomputed, main="local GRAY", xlab="AAC")
hist(GRAY2017_updated@sensitivity$profiles$auc_recomputed, main="GRAY 2017", xlab="AAC")
load(file.path(psets.path, "GRAY2013_updated.RData"), verbose=T)
hist(GRAY2013_updated@sensitivity$profiles$auc_recomputed, main="GRAY 2013", xlab="AAC")


