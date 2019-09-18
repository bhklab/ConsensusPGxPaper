if(!require(xtable)){install.packages("xtable");library(xtable)}
path.results <- "../results"
load("~/Google Drive/PSets/GDSC.RData")
ids <- unlist(strsplit(GDSC@drug[which(rownames(GDSC@drug) == "AZD6482"), "drugid"], split="/"))

tt <- GDSC@sensitivity$profiles[grep(paste0("drugid_", ids[1]), rownames(GDSC@sensitivity$info)), "auc_recomputed"]
ss <- GDSC@sensitivity$profiles[grep(paste0("drugid_", ids[2]), rownames(GDSC@sensitivity$info)),"auc_recomputed"]

names(tt) <- GDSC@sensitivity$info[grep(paste0("drugid_", ids[1]), rownames(GDSC@sensitivity$info)), "cellid"]
names(ss) <- GDSC@sensitivity$info[grep(paste0("drugid_", ids[2]), rownames(GDSC@sensitivity$info)), "cellid"]

GDSC1000_rep <- length(intersect(names(tt), names(ss)))

dd <- apply(CTRPv2@sensitivity$n, 2, function(x){length(which(x > 1))})
ctrpv2_rep <- sum(dd)

dd <- apply(gCSI@sensitivity$n, 2, function(x){length(which(x > 1))})
gCSI_rep <- sum(dd)

load("~/Google Drive/Psets/GRAY_kallisto.RData")
dd <- apply(GRAY@sensitivity$n, 2, function(x){length(which(x > 1))})
GRAY_rep <- sum(dd)

xx <- as.data.frame(matrix(NA,
                           ncol=7, nrow=5,
                           dimnames=(list(c("CCLE", "CTRPv2", "gCSI", "GDSC1000", "GRAY"),
                                          c("Dataset", "Compounds", "Cell lines", "Experiments", "Replicates", "Assay", "Molecular Profiles")))))
xx["CCLE", ] <- c(pSetName(CCLE), length(drugNames(CCLE)), length(cellNames(CCLE)), nrow(CCLE@sensitivity$info), "_", "CellTiter-Glo", paste(setdiff(names(CCLE@molecularProfiles), c("isoforms", "rnaseq.counts", "isoforms.counts")), collapse=", "))
xx["CTRPv2", ] <- c(pSetName(CTRPv2), length(drugNames(CTRPv2)), length(cellNames(CTRPv2)), nrow(CTRPv2@sensitivity$info), ctrpv2_rep, "CellTiter-Glo", paste(names(CTRPv2@molecularProfiles), collapse=", "))
xx["gCSI", ] <- c(pSetName(gCSI), length(drugNames(gCSI)), length(cellNames(gCSI)), nrow(gCSI@sensitivity$info), gCSI_rep, "CellTiter-Glo", paste(names(gCSI@molecularProfiles), collapse=", "))
xx["GDSC1000", ] <- c(pSetName(GDSC1000), length(drugNames(GDSC1000)), length(cellNames(GDSC1000)), nrow(GDSC1000@sensitivity$info), GDSC1000_rep, "Syto60", paste(names(GDSC1000@molecularProfiles), collapse=", "))
xx["GRAY", ] <- c(pSetName(GRAY), length(drugNames(GRAY)), length(cellNames(GRAY)), nrow(GRAY@sensitivity$info), GRAY_rep, "CellTiter-Glo", paste(setdiff(names(GRAY@molecularProfiles), c("isoforms", "rnaseq.counts", "isoforms.counts")), collapse=", "))
xtable::print.xtable(xtable::xtable(xx), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.results, "datasets.tex"), append=FALSE)

gCSI_v1_2_metric <- read.csv("../data/gCSI_GRmetrics_v1.2.tsv", sep="\t")
length(unique(gCSI_v1_2_metric$DrugName))
#[1] 39
length(unique(gCSI_v1_2_metric$CellLineName))
#[1] 569
nrow(gCSI_v1_2_metric)
#[1] 14284


#table 2
datasets <- c("CCLE", "gCSI_V1", "CTRPv2", "GDSC1000")
sensitivity_metric="auc_recomputed"
source("pset_loading.R")

ll <- length(datasets)
drugs_all <- drugs <- NULL
test <- list()
j <- NULL
for(i in 1:ll){
  j <- combn(ll, ll - 1)[,i]
  drugs <- union(drugs, intersectList(sapply(j, function(x){colnames(drug.sensitivity(datasets[x]))})))
  drugs_all <- union(drugs, intersectList(sapply(j, function(x){drugNames(pset(datasets[x]))})))
  test[[i]] <- intersectList(sapply(j, function(x){colnames(drug.sensitivity(datasets[x]))}))
  names(test)[i] <- paste(datasets[j], collapse="_")
}
setdiff(drugs_all, drugs)
xx <- as.data.frame(matrix(NA,
                           ncol=3, nrow=length(drugs)))
colnames(xx) <- c("Compond", "Target", "Chemo")
rownames(xx) <- xx[,1] <- drugs
xx[,"Chemo"] <- ""
xx[which(rownames(xx) %in% chemo), "Chemo"] <- "X"
ii <- intersect(drugs, rownames(CTRPv2@drug))
xx[ii, 2] <- CTRPv2@drug[ii,"target_or_activity_of_compound"]
ii <- setdiff(drugs, rownames(CTRPv2@drug))
xx[ii, 2] <- CCLE@drug[ii, "Mechanism.of.action"]
xx[,2] <- gsub("inhibitor of ", "", xx[,2])
xx[16, 2] <- "DNA replication"
xx[22, 2] <- "thymidylate synthase"
xtable::print.xtable(xtable::xtable(xx), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.results, "drugs.tex"), append=FALSE)



##Supplementary Table 1
xx <- read.csv("../results/comparison_rnaseq.csv", row.names=1)
nrow(xx)
length(which(!is.na(xx$rCI)))
xx <- xx[which(xx$rCI_pvalue <= 0.05), ]
dim(xx)
#colnames(xx)[4:ncol(xx)] <-c("rCI_pval", "NO_vs_ALL", "NO_vs_SOME", "")
xtable::print.xtable(xtable::xtable(xx, digits=c(0, 0, 0, 2, -1, -1, -1, -1, -1)), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(path.results, "supp_table_1.tex"), append=FALSE)
