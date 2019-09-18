source("foo.R")
library(PharmacoGx)
library(Biobase)

path.result <- "../results"
psets.path <- "~/Google Drive/Psets"
#psets.path <- "../Psets"

if("CCLE" %in% datasets){
  load(file.path(psets.path, "CCLE_kallisto.RData"), verbose=T)
  ccle.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE, sensitivity.measure=sensitivity_metric, summary.stat="median"))
  ccle_no_effect <- names(which(apply(ccle.drug.sensitivity, 2, function(x){length(which(x > .25))}) <= 5))
  drug.mad.ccle <- apply(ccle.drug.sensitivity, 2, mad, na.rm=TRUE)
  ccle_narrow_effect <- names(which(drug.mad.ccle <= .13))
  ccle_broad_effect <- names(which(drug.mad.ccle > .13))
  ccle.drug.sensitivity <- ccle.drug.sensitivity[,-which(colnames(ccle.drug.sensitivity) %in% ccle_no_effect)]
  not_hemato <- which(!rownames(ccle.drug.sensitivity) %in% rownames(CCLE@cell)[which(CCLE@cell$tissueid == "haematopoietic_and_lymphoid_tissue")])
  ccle_drug_response_var <- apply(ccle.drug.sensitivity, 2, var, na.rm=T)
  #ccle.drug.sensitivity <- ccle.drug.sensitivity[not_hemato, ]
  
  ccle.rna <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="rna", summary.stat="median"))
  ccle.rnaseq <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", summary.stat="median"))
  ccle.mutation <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="mutation", summary.stat="or"))
  ccle.cnv <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="cnv", summary.stat="median"))
}
if("GDSC1000" %in% datasets){
  load(file.path(psets.path, "GDSC1000.RData"), verbose=T)
  gdsc1000.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC1000, sensitivity.measure=sensitivity_metric, summary.stat="median"))
  gdsc1000_no_effect <- names(which(apply(gdsc1000.drug.sensitivity, 2, function(x){length(which(x > .25))}) <= 5))
  drug.mad.gdsc1000 <- apply(gdsc1000.drug.sensitivity, 2, mad, na.rm=TRUE)
  gdsc1000_narrow_effect <- names(which(drug.mad.gdsc1000 <= .13))
  gdsc1000_broad_effect <- names(which(drug.mad.gdsc1000 > .13))
  gdsc1000.drug.sensitivity <- gdsc1000.drug.sensitivity[,-which(colnames(gdsc1000.drug.sensitivity) %in% gdsc1000_no_effect)]
  not_hemato <- which(!rownames(gdsc1000.drug.sensitivity) %in% rownames(GDSC1000@cell)[which(GDSC1000@cell$tissueid == "haematopoietic_and_lymphoid_tissue")])
  gdsc1000_drug_response_var <- apply(gdsc1000.drug.sensitivity, 2, var, na.rm=T)
  #gdsc1000.drug.sensitivity <- gdsc1000.drug.sensitivity[not_hemato, ]
  
  gdsc.rna <- exprs(summarizeMolecularProfiles(pSet=GDSC1000, mDataType="rna", summary.stat="median"))
  #load("../data/GDSC1000_mutation_gene_level.RData", verbose=T)
  #gdsc1000.mutation <- vcf_exprs_gene_level
  #gdsc1000.mutation.fData <- vcf_fData_gene_level
  load("../data/GDSC1000_mutation_published_lorio.RData", verbose=T)
  gdsc1000.mutation <- mmc4_mutations
  gdsc1000.mutation.fData <- mmc4_fData
}
if("gCSI_V1" %in% datasets){
  load(file.path(psets.path, "gCSI_V1.RData"), verbose=T)
  gcsi.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=gCSI, sensitivity.measure=sensitivity_metric, summary.stat="median"))
  if(sensitivity_metric == "GReff"){
    gcsi.drug.sensitivity <- -gcsi.drug.sensitivity
  }
  gcsi_no_effect <- names(which(apply(gcsi.drug.sensitivity, 2, function(x){length(which(x > .25))}) <= 5))
  drug.mad.gcsi <- apply(gcsi.drug.sensitivity, 2, mad, na.rm=TRUE)
  gcsi_narrow_effect <- names(which(drug.mad.gcsi <= .13))
  gcsi_broad_effect <- names(which(drug.mad.gcsi > .13))
  gcsi_drug_response_var <- apply(gcsi.drug.sensitivity, 2, var, na.rm=T)
  if(length(gcsi_no_effect > 0)){
    gcsi.drug.sensitivity <- gcsi.drug.sensitivity[,-which(colnames(gcsi.drug.sensitivity) %in% gcsi_no_effect)]
  }
  not_hemato <- which(!rownames(gcsi.drug.sensitivity) %in% rownames(gCSI@cell)[which(gCSI@cell$tissueid == "haematopoietic_and_lymphoid_tissue")])
  #gcsi.drug.sensitivity <- gcsi.drug.sensitivity[not_hemato, ]
  
  gcsi.rnaseq <- exprs(summarizeMolecularProfiles(pSet=gCSI, mDataType="rnaseq", summary.stat="median"))
  gcsi.cnv <- exprs(summarizeMolecularProfiles(pSet=gCSI, mDataType="cnv", summary.stat="median"))
  #gcsi.mutation <- exprs(summarizeMolecularProfiles(pSet=gCSI, mDataType="mutation", summary.stat="or"))
  gcsi.mutation <- exprs(gCSI@molecularProfiles$mutation)
  #  if(!"CCLE" %in% ls()){
  #    load(file.path(psets.path, "CCLE_kallisto.RData"), verbose=T)
  #    ccle.mutation <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="mutation", summary.stat="or"))
  #  }
}

if("CTRPv2" %in% datasets){
  load(file.path(psets.path, "CTRPv2.RData"), verbose=T)
  ctrpv2.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=CTRPv2, sensitivity.measure=sensitivity_metric, summary.stat="median"))
  ctrpv2_no_effect <- names(which(apply(ctrpv2.drug.sensitivity, 2, function(x){length(which(x > .25))}) <= 5))
  drug.mad.ctrpv2 <- apply(ctrpv2.drug.sensitivity, 2, mad, na.rm=TRUE)
  ctrpv2_narrow_effect <- names(which(drug.mad.ctrpv2 <= .13))
  ctrpv2_broad_effect <- names(which(drug.mad.ctrpv2 > .13))
  ctrpv2.drug.sensitivity <- ctrpv2.drug.sensitivity[,-which(colnames(ctrpv2.drug.sensitivity) %in% ctrpv2_no_effect)]
  not_hemato <- which(!rownames(ctrpv2.drug.sensitivity) %in% rownames(CTRPv2@cell)[which(CTRPv2@cell$tissueid == "haematopoietic_and_lymphoid_tissue")])
  ctrpv2_drug_response_var <- apply(ctrpv2.drug.sensitivity, 2, var, na.rm=T)
  #ctrpv2.drug.sensitivity <- ctrpv2.drug.sensitivity[not_hemato, ]
  
  if(!"CCLE" %in% ls()){
    load(file.path(psets.path, "CCLE_kallisto.RData"), verbose=T)
    ccle.rna <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="rna", summary.stat="median"))
    ccle.rnaseq <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", summary.stat="median"))
    ccle.mutation <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="mutation", summary.stat="or"))
    ccle.cnv <- exprs(summarizeMolecularProfiles(pSet=CCLE, mDataType="cnv", summary.stat="median"))
  }
}
if("GDSC1000_updatedGRset" %in% datasets){
  load(file.path(psets.path, "GDSC1000_updatedGRset.RData"), verbose=T)
  gdsc1000_GR.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC1000_updatedGRset, sensitivity.measure=sensitivity_metric, summary.stat="median"))
  gdsc_GR.rna <- exprs(summarizeMolecularProfiles(pSet=GDSC1000_updatedGRset, mDataType="rna", summary.stat="median"))
  load("../data/GDSC1000_mutation_published_lorio.RData", verbose=T)
  gdsc1000.mutation <- mmc4_mutations
  gdsc1000.mutation.fData <- mmc4_fData
}

if("AZ" %in% datasets){
  load(file.path(psets.path, "AZ.RData"), verbose=T)
  AZ.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=AZ, sensitivity.measure=sensitivity_metric, summary.stat="median"))
}

load(file.path(psets.path, "GDSC.RData"), verbose=T)
gdsc.mutation <- exprs(summarizeMolecularProfiles(pSet=GDSC, mDataType="mutation", summary.stat="or"))
gdsc.cnv <- exprs(summarizeMolecularProfiles(pSet=GDSC, mDataType="cnv", summary.stat="median"))
gdsc.fusion <- exprs(summarizeMolecularProfiles(pSet=GDSC, mDataType="fusion", summary.stat="or"))

xx <- featureInfo(GDSC, "fusion")
xx$Symbol <-rownames(featureInfo(GDSC, "fusion"))
featureInfo(GDSC, "fusion") <- xx
#load(file.path(psets.path, "GRAY_kallisto.RData"), verbose=T)

#molecular profiles

chemo <- c("5-FU", "Bortezomib", "Doxorubicin", "etoposide", "Methotrexate", "Docetaxel", "paclitaxel", "Gemcitabine")
