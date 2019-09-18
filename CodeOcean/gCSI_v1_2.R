options(stringsAsFactors = FALSE)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

gCSI_V1_2_raw <- read.csv("../data/gCSI_GRvalues_v1.2.tsv", sep="\t")
#gCSI_v1_1_metric <- read.csv("../data/gCSI_GRmetrics_v1.tsv", sep="\t")
gCSI_v1_2_metric <- read.csv("../data/gCSI_GRmetrics_v1.2.tsv", sep="\t")

fn <- "~/Google Drive/Psets/gCSI_V1.RData"
if(file.exists(fn)){
  load(fn)
}else{
  load("~/Google Drive/Psets/gCSI_kallisto.RData")
}
cell_annotation_all <- read.csv("~/Documents/GitHub/PharmacoGx-private/inst/extdata/cell_annotation_all.csv", na.string=c("", " "))
drug_annotation_all <- read.csv("~/Documents/GitHub/PharmacoGx-private/inst/extdata/drug_annotation_all.csv", na.string=c("", " "))

drugid <- drug_annotation_all$unique.drugid
drugs <- unique(gCSI_V1_2_raw$DrugName)
min(length(drugid), length(drugs)) - length(intersect(drugs, drugid))

drugs_curated <- tolower(gsub(badchars, "", drugs))
drugid_curated <- tolower(gsub(badchars, "", drugid))
min(length(drugs_curated), length(drugid_curated)) - length(intersect(drugs_curated, drugid_curated))
sort(drugs[which(drugs_curated %in% setdiff(drugs_curated, drugid_curated))])

drugs_new <- drugs[which(drugs_curated %in% setdiff(drugs_curated, drugid_curated))]
gcsi_drugid <- drugid[match(drugs_curated, drugid_curated)]
gcsi_drugid[which(is.na(gcsi_drugid))] <- drugs_new

drug_slot <- gCSI@drug
length(setdiff(rownames(drug_slot), gcsi_drugid))
xx <- matrix(NA, ncol=ncol(drug_slot), nrow=length(setdiff(gcsi_drugid, rownames(drug_slot))),
             dimnames=list(setdiff(gcsi_drugid, rownames(drug_slot)), colnames(drug_slot)))
xx[,"originalDrugId"] <- drugs[which(gcsi_drugid %in% rownames(xx))]
xx[,"drugid"] <- rownames(xx)
drug_slot <- rbind(drug_slot, xx)
dd <- drugs[match(rownames(drug_slot), gcsi_drugid)]
dd[which(is.na(dd))] <- rownames(drug_slot)[which(is.na(dd))]
drug_slot[,"originalDrugId"] <- dd

drug_curation <- gCSI@curation$drug
xx_curation <- cbind(rownames(xx), xx[,"originalDrugId"])
rownames(xx_curation) <- rownames(xx)
colnames(xx_curation) <- colnames(drug_curation)
drug_curation <- rbind(drug_curation, xx_curation)
rownames(drug_curation) <- rownames(drug_slot)
drug_curation[,"gCSI.drugid"] <- drug_slot[,"originalDrugId"]

cellid <- cell_annotation_all$unique.cellid#rownames(gCSI@cell)
cells <- unique(gCSI_V1_2_raw$CellLineName)
min(length(cellid), length(cells)) - length(intersect(cells, cellid))

cells_curated <- tolower(gsub(badchars, "", cells))
cellid_curated <- tolower(gsub(badchars, "", cellid))
min(length(cells_curated), length(cellid_curated)) - length(intersect(cells_curated, cellid_curated))
##two cell lines has been used with different names in Marc's new data
cells[which(cells_curated %in% cells_curated[which(duplicated(cells_curated))])]
#[1] "HCT-116" "JeKo-1"  "HCT 116" "JEKO-1"
#cells[which(cells == "HCT 116")] <- "HCT-116"
#cells[which(cells == "JEKO-1")] <- "JeKo-1"
#cells <- unique(cells)
cells_curated <- tolower(gsub(badchars, "", cells))
cellid_curated <- tolower(gsub(badchars, "", cellid))
min(length(cells_curated), length(cellid_curated)) - length(intersect(cells_curated, cellid_curated))
gcsi_cellid <- cellid[match(cells_curated, cellid_curated)]

cell_slot <- gCSI@cell
xx <- matrix(NA, ncol=ncol(cell_slot), nrow=length(setdiff(gcsi_cellid, rownames(cell_slot))),
             dimnames=list(setdiff(gcsi_cellid, rownames(cell_slot)), colnames(cell_slot)))
xx[,"CellLineName"] <- cells[which(gcsi_cellid %in% rownames(xx))]
xx[,"tissueid"] <- cell_annotation_all$unique.tissueid[which(cell_annotation_all$unique.cellid %in% rownames(xx))]
xx[,"TissueMetaclass"] <- gCSI_V1_2_raw$PrimaryTissue[match(xx[,"CellLineName"], gCSI_V1_2_raw$CellLineName)]
cell_slot <- rbind(cell_slot, xx)
cc <- cells[match(rownames(cell_slot), gcsi_cellid)]
cc[which(is.na(cc))] <- rownames(cell_slot)[which(is.na(cc))]
cell_slot[,"CellLineName"] <- cc

cell_curation <- gCSI@curation$cell
xx_curation <- cbind(rownames(xx), xx[,"CellLineName"], NA)
rownames(xx_curation) <- rownames(xx)
colnames(xx_curation) <- colnames(cell_curation)
cell_curation <- rbind(cell_curation, xx_curation)
cell_curation[,"gCSI.cellid"] <- cell_slot[,"CellLineName"]

tissue_curation <- gCSI@curation$tissue
xx_curation <- cbind(xx[,"tissueid"], xx[,"TissueMetaclass"], NA)
rownames(xx_curation) <- rownames(xx)
colnames(xx_curation) <- colnames(tissue_curation)
tissue_curation <- rbind(tissue_curation, xx_curation)

load("../data/gCSI_auc_v1_2.RData", verbose=T)
gCSI_AUC <- cbind(gCSI_AUC, "cellid"=NA, "drugid"=NA)
gCSI_AUC$drugid <- drug_curation$unique.drugid[match(gCSI_AUC$DrugName, drug_curation$gCSI.drugid)]
#gCSI_AUC$CellLineName[which(gCSI_AUC$CellLineName == "HCT 116")] <- "HCT-116"
#gCSI_AUC$CellLineName[which(gCSI_AUC$CellLineName == "JEKO-1")] <- "JeKo-1"
gCSI_AUC$cellid <- cell_curation$unique.cellid[match(gCSI_AUC$CellLineName, cell_curation$gCSI.cellid)]
rownames(gCSI_AUC) <- sprintf("%s_%s_%s", gCSI_AUC$cellid, gCSI_AUC$drugid, gCSI_AUC$experimentid)
old_experimentds <- which(gCSI@sensitivity$info$drugid %in% setdiff(rownames(drug_slot), gcsi_drugid))
sensitivity_info <- rbind(gCSI_AUC[, c("cellid", "drugid", "experimentid")], cbind(gCSI@sensitivity$info[old_experimentds, c("cellid", "drugid")], "experimentid"=NA))
sensitivity_profiles <-  rbind(gCSI_AUC[,"auc_recomputed", drop=FALSE], gCSI@sensitivity$profiles[old_experimentds, "auc_recomputed", drop=FALSE])
gCSI_V1_2 <- PharmacoSet(molecularProfiles=gCSI@molecularProfiles,
                       name="gCSI",
                       cell=cell_slot,
                       drug=drug_slot,
                       sensitivityInfo=sensitivity_info,
                       sensitivityRaw=NULL,
                       sensitivityProfiles=sensitivity_profiles,
                       sensitivityN=NULL,
                       curationCell=cell_curation,
                       curationDrug=drug_curation,
                       curationTissue=tissue_curation,
                       datasetType="sensitivity")
warnings()
gCSI <- gCSI_V1_2
gCSI@molecularProfiles <- gCSI@molecularProfiles[c("rnaseq", "cnv")]
gCSI@sensitivity$profiles$auc_recomputed <- as.numeric(gCSI@sensitivity$profiles$auc_recomputed)
save(gCSI, file="~/Google Drive/Psets/gCSI_V1.RData")
load("~/Google Drive/Psets/gCSI_V1.RData")
gCSI@sensitivity$n <- PharmacoGx:::.summarizeSensitivityNumbers(gCSI)
save(gCSI, file="~/Google Drive/Psets/gCSI_V1.RData")

##For now it is not possible to integrate other metrics into the pset due to replicates issue
drug_curation <- gCSI@curation$drug
cell_curation <- gCSI@curation$cell
gCSI_V1_2_metrics <- read.csv("../data/gCSI_GRmetrics_v1.2.tsv", sep="\t")
gCSI_V1_2_metrics <- cbind(gCSI_V1_2_metrics, "cellid"=NA, "drugid"=NA)

gCSI_V1_2_metrics$drugid <- gCSI@curation$drug$unique.drugid[match(gCSI_V1_2_metrics$DrugName, gCSI@curation$drug$gCSI.drugid)]
table(is.na(gCSI_V1_2_metrics$drugid))

#gCSI_V1_2_metrics$CellLineName[which(gCSI_V1_2_metrics$CellLineName == "HCT 116")] <- "HCT-116"
#gCSI_V1_2_metrics$CellLineName[which(gCSI_V1_2_metrics$CellLineName == "JEKO-1")] <- "JeKo-1"
gCSI_V1_2_metrics$cellid <- gCSI@curation$cell$unique.cellid[match(gCSI_V1_2_metrics$CellLineName, gCSI@curation$cell$gCSI.cellid)]
table(is.na(gCSI_V1_2_metrics$cellid))
rownames(gCSI_V1_2_metrics) <- sprintf("%s_%s_%s", gCSI_V1_2_metrics$cellid, gCSI_V1_2_metrics$drugid, gCSI_V1_2_metrics$ExperimentNumber)
##temporary solution:remove duplicates
#dd <- which(!duplicated(xx))
#xx <- xx[dd]
#gCSI_V1_2_metrics <- gCSI_V1_2_metrics[dd,]
ss <- sprintf("%s_%s_%s", gCSI@sensitivity$info$cellid, gCSI@sensitivity$info$drugid, gCSI@sensitivity$info$experimentid)
gCSI@sensitivity$profiles[rownames(gCSI_V1_2_metrics), "GR_AOC"] <- gCSI_V1_2_metrics$GR_AOC
gCSI@sensitivity$profiles[rownames(gCSI_V1_2_metrics), "GReff"] <- gCSI_V1_2_metrics$GR_05uM_fit
save(gCSI, file="~/Google Drive/Psets/gCSI_V1.RData")

###update doubling time
if(!"DoublingTime" %in% colnames(gCSI@cell)){
  gCSI_V1_2_raw <- read.csv("../data/gCSI_GRvalues_v1.2.tsv", sep="\t")
  gCSI_V1_2_raw <- gCSI_V1_2_raw[!duplicated(gCSI_V1_2_raw$CellLineName),]
  dim(gCSI_V1_2_raw)
  mm <- match(tolower(gsub(badchars, "", gCSI_V1_2_raw$CellLineName)), tolower(gsub(badchars, "", gCSI@curation$cell$gCSI.cellid)))
  gCSI@cell$DoublingTime[mm] <- gCSI_V1_2_raw$doublingtime
  save(gCSI, file="~/Google Drive/Psets/gCSI_V1.RData")
}
