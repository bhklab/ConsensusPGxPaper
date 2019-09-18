options(stringsAsFactors = FALSE)
serial <- True
gCSI_GR_AOC <- read.csv("../data/gCSI_GRvalues_v1.2.tsv", sep="\t")

uids <- unique(sprintf("%s__%s__%s", gCSI_GR_AOC$DrugName,
          gCSI_GR_AOC$CellLineName,
          gCSI_GR_AOC$ExperimentNumber))
gCSI_AUC <- matrix(NA, ncol=5, nrow=length(uids))
colnames(gCSI_AUC) <- c("DrugName", "CellLineName", "experimentid", "rep_no", "auc_recomputed")
rownames(gCSI_AUC) <- uids
if(!serial){
  library(parallel)
  nthread <- parallel::detectCores()
  nthread
  mres <- parallel::mclapply(1:length(uids), function(i, gCSI_GR_AOC){
    uid <- unlist(strsplit(uids[i], "__"))
    xx <- gCSI_GR_AOC[which(gCSI_GR_AOC$DrugName==uid[1] &
                              gCSI_GR_AOC$CellLineName==uid[2] &
                              gCSI_GR_AOC$ExperimentNumber==uid[3]),]
    #PharmacoGx::drugDoseResponseCurve(concentrations=xx$log10Concentration,
    #                                  viabilities=xx$relative_cell_count*100,
    #                                  conc_as_log=T)

    auc <- PharmacoGx::computeAUC(concentration=xx$log10Concentration,
                                  viability=xx$relative_cell_count,
                                  conc_as_log=T,
                                  viability_as_pct=F)
    return(c(uid, auc, 1))
  }, gCSI_GR_AOC=gCSI_GR_AOC, mc.cores=nthread)
}else{
  for(i in which(is.na(gCSI_AUC[,1]))){#1:length(uids)){
    uid <- unlist(strsplit(uids[i], "__"))
    xx <- gCSI_GR_AOC[which(gCSI_GR_AOC$DrugName==uid[1] &
                              gCSI_GR_AOC$CellLineName==uid[2] &
                              gCSI_GR_AOC$ExperimentNumber==uid[3]),]
    #PharmacoGx::drugDoseResponseCurve(concentrations=xx$log10Concentration,
    #                                  viabilities=xx$relative_cell_count*100,
    #                                  conc_as_log=T)

    auc <- PharmacoGx::computeAUC(concentration=xx$log10Concentration,
                                  viability=xx$relative_cell_count,
                                  conc_as_log=T,
                                  viability_as_pct=F)
    gCSI_AUC[uids[i],] <- c(uid, , auc, 1)


  }
}
gCSI_AUC <- as.data.frame(gCSI_AUC, stringsAsFactors=FALSE)
save(gCSI_AUC, file="../data/gCSI_auc_v1_2.RData")

