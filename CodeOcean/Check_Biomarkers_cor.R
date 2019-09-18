path.results <- "../results_gcsi_gr"
sensitivity_metric <- "GR_AOC"
delta=0.22
#source("pset_loading.R")

load(file.path("../results/G2P_clarified_hits_cnv.RData"), verbose=T)
#load(file.path(pset.download.path, "GDSC1000.RData"), verbose=T)

plot_mutation <- function(drug_sensitivity_vec, molecular_profile_vec, cell.lines, metric, dataset, drug, feature, type, ylim){
  boxplot(drug_sensitivity_vec~molecular_profile_vec, outline=F,
          ylim=ylim, ylab=sprintf("%s %s", drug, metric), xlab=sprintf("%s %s", feature, type), main=sprintf("%s", dataset))
  stripchart(drug_sensitivity_vec~molecular_profile_vec , vertical = TRUE,
             method = "jitter", add = TRUE, pch = 20, col = 'blue')
  mci <- paired.concordance.index(predictions=molecular_profile_vec,
                                  observations=drug_sensitivity_vec,
                                  delta.pred=0,
                                  delta.obs=delta,
                                  alternative="greater",
                                  logic.operator="and")["cindex"]
  ci <- paired.concordance.index(predictions=molecular_profile_vec,
                                 observations=drug_sensitivity_vec,
                                 delta.pred=0,
                                 delta.obs=0,
                                 alternative="greater",
                                 logic.operator="and")["cindex"]
  legend("topright", legend=sprintf("mCI=%.2f\nCI=%.2f", mci, ci), bty="n")
}
plot_expr <- function(drug_sensitivity_vec, molecular_profile_vec, cell.lines, metric, dataset, drug, feature, type, ylim){

  mci <- paired.concordance.index(predictions=molecular_profile_vec,
                                  observations=drug_sensitivity_vec,
                                  delta.pred=0,
                                  delta.obs=delta,
                                  alternative="greater",
                                  logic.operator="and")["cindex"]
  ci <- paired.concordance.index(predictions=molecular_profile_vec,
                                 observations=drug_sensitivity_vec,
                                 delta.pred=0,
                                 delta.obs=0,
                                 alternative="greater",
                                 logic.operator="and")["cindex"]
  xlim=range(molecular_profile_vec, na.rm=T)
  xlim[2] <- xlim[2] + 1
  xlim[1] <- xlim[1] - 1
  if(sensitivity_metric=="GR_AOC"){
    dup <- as.numeric(sapply(cell.lines, function(x){strsplit(x, ":")[[1]][2]}))

    mycol2=RColorBrewer::brewer.pal(name="Set1", n=5)
    cc <- colorRampPalette(c(mycol2[2], mycol2[1]))
    i=10
    ii <- which(dup > quantile(dup, .2) & dup < quantile(dup, .8))
    cc2 <- cc(i)[as.numeric(cut((dup-min(dup))/(max(dup)-min(dup)), breaks=i))]
    cc3 <- vector(length=length(dup))
    cc3[which(dup <= quantile(dup, .2))] <- mycol2[2]
    cc3[which(dup > quantile(dup, .2) & dup < quantile(dup, .8))] <- cc(i)[as.numeric(cut((dup[ii]-min(dup[ii]))/(max(dup[ii])-min(dup[ii])), breaks=i))]
    cc3[which(dup >= quantile(dup, .8))] <- mycol2[1]
    #barplot(rep(1,length(dup)), col=cc2)
    #corr <- cor(dup, molecular_profile_vec, use="pairwise.complete.obs", method="spearman")
    pheno_corr <- paired.concordance.index(predictions=molecular_profile_vec,
                                   observations=dup,
                                   delta.pred=0,
                                   delta.obs=0,
                                   alternative="greater",
                                   logic.operator="and")[c("cindex", "p.value")]
    sensitivity_corr <- paired.concordance.index(predictions=dup,
                                           observations=drug_sensitivity_vec,
                                           delta.pred=0,
                                           delta.obs=delta,
                                           alternative="greater",
                                           logic.operator="and")[c("cindex", "p.value")]
    myScatterPlot(y=drug_sensitivity_vec,
                  x=molecular_profile_vec,
                  legend.label=sprintf("mCI=%.2f", mci), method="plain", col=alpha(cc2, 0.8),
                  ylim=ylim, xlim=xlim, ylab=sprintf("%s %s", drug, metric), xlab=sprintf("%s %s", feature, type), main=paste0(sprintf("cells doubling time ~ %s %s, CI: %.2f, %2.e\n", feature, type, pheno_corr$cindex, pheno_corr$p.value),
                                                                                                                               sprintf("response to %s ~ cells doubling time, mCI: %.2f, %2.e\n", drug, sensitivity_corr$cindex, sensitivity_corr$p.value)), cex.main=.8)
    legend("topright", legend=sprintf("mCI=%.2f", mci), cex=1.1, bty="n")

    xx <- order(drug_sensitivity_vec, decreasing=T)
    xxx <- order(molecular_profile_vec, decreasing=T)
    l <- length(xx)
    ii <- intersect(xx[1:3], xxx[1:3])
    if(length(ii)==0){ii <- xx[1:3]}
    iii <- intersect(xx[(l-3):l], xxx[(l-3):l])
    if(length(iii)==0){iii <- xx[(l-3):l]}
    xx <- unionList(ii, iii)

    #text(molecular_profile_vec[xx], drug_sensitivity_vec[xx], labels=cell.lines[xx], cex= 0.7, pos=rep(c(1, 3), length(xx)%/%2 + 1))
    }else{
      myScatterPlot(y=drug_sensitivity_vec,
                    x=molecular_profile_vec,
                    legend.label=sprintf("mCI=%.2f", mci), method="plain",#, method="transparent", transparency=.2,
                    ylim=ylim, xlim=xlim, ylab=sprintf("%s %s", drug, metric), xlab=sprintf("%s %s", feature, type))#, main=sprintf("%s", dataset))

      legend("topright", legend=sprintf("mCI=%.2f", mci), cex=1.1, bty="n")
    }

}
plot_feature_drug_rel <- function(drug_sensitivity_vec, molecular_profile_vec, cell.lines, metric, dataset, drug, feature, type, ylim){
  if(sensitivity_metric=="GR_AOC"){
    ylim=c(-.5, 1.5)
  }else{
    ylim=c(0,1)
  }
  if(type=="mutation"){
    plot_mutation(drug_sensitivity_vec=drug_sensitivity_vec, molecular_profile_vec=molecular_profile_vec, cell.lines=cell.lines, metric=metric, dataset=dataset, drug=drug, feature=feature, type=type, ylim=ylim)
  }else{
    plot_expr(drug_sensitivity_vec=drug_sensitivity_vec, molecular_profile_vec=molecular_profile_vec, cell.lines=cell.lines, metric=metric, dataset=dataset, drug=drug, feature=feature, type=type, ylim=ylim)
  }
}

for(hit in 1:nrow(G2P_clarified_hits)){
  drug <- G2P_clarified_hits$compound[hit]
  feature <- G2P_clarified_hits$gene[hit]
  type <- G2P_clarified_hits$ttype[hit]

  xx <- G2P_clarified_hits[hit,grep("_mCI$", colnames(G2P_clarified_hits))]
  if(length(xx)==1){names(xx) <- grep("_mCI$", colnames(G2P_clarified_hits), value = T)}
  dd <- sapply(names(xx)[which(!is.na(xx))], function(x){strsplit(x, "_mCI")[[1]]})

  pdf(file.path(path.results, sprintf("Plot_%s_%s_%s.pdf", drug, feature, type)), width=5 * 2, height=5 * (length(dd) %/% 2 + length(dd) %% 2))
  par(mfrow=c(length(dd) %/% 2 + length(dd) %% 2, 2))
  for(dataset.name in dd){
    id <- rownames(featureInfoType(dataset.name, type))[match(feature, featureInfoType(dataset.name, type)[, "Symbol"])]
    molecular_profile <- molecular.profile(dataset.name, type)
    #par(mfrow=c(1, 3))
    drug_sensitivity <- drug.sensitivity(dataset.name)
    cell.lines <- intersect(colnames(molecular_profile), rownames(drug_sensitivity))
    cell.lines <- cell.lines[which(!is.na(drug_sensitivity[cell.lines, drug]))]
    cell.lines <- cell.lines[which(!is.na(molecular_profile[id ,cell.lines]))]
    sensitivity_metric="GR_AOC"
    delta=.22
    if(sensitivity_metric=="GR_AOC"){
      dup <- gCSI@cell[cell.lines, "DoublingTime"]
      labels=paste(cell.lines, dup, sep=":")
    }else{
      labels=cell.lines
    }

    plot_feature_drug_rel(molecular_profile_vec=as.numeric(molecular_profile[id ,cell.lines]),
                          drug_sensitivity_vec=as.numeric(drug_sensitivity[cell.lines, drug]),
                          cell.lines=labels,
                          metric=sensitivity_metric, dataset=strsplit(dataset.name, "_")[[1]][1], drug=drug, feature=feature, type=type, ylim=ylim)
    gcsi.drug.sensitivity.auc <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=gCSI, sensitivity.measure="auc_recomputed", summary.stat="median"))
    sensitivity_metric="AAC"
    delta=.13
    plot_feature_drug_rel(molecular_profile_vec=as.numeric(molecular_profile[id ,cell.lines]),
                          drug_sensitivity_vec=as.numeric(gcsi.drug.sensitivity.auc[cell.lines, drug]),
                          cell.lines=cell.lines,
                          metric="AAC", dataset=strsplit(dataset.name, "_")[[1]][1], drug=drug, feature=feature, type=type, ylim=ylim)
  }
  dev.off()
  # drug_sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=pset(dataset.name), sensitivity.measure="AAC"))
  # cell.lines <- intersect(colnames(molecular_profile), rownames(drug_sensitivity))
  # plot_feature_drug_rel(drug_sensitivity_vec=as.numeric(molecular_profile[id ,cell.lines]),
  #                       molecular_profile_vec=as.numeric(drug_sensitivity[cell.lines, drug]),
  #                       metric="AAC", dataset=dataset.name, drug=drug, feature=feature, type=type)
  #
  #
  # drug_sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=pset("GDSC1000"), sensitivity.measure="auc_recomputed"))
  # cell.lines <- intersect(colnames(molecular_profile), rownames(drug_sensitivity))
  # plot_feature_drug_rel(drug_sensitivity_vec=as.numeric(molecular_profile[id ,cell.lines]),
  #                       molecular_profile_vec=as.numeric(drug_sensitivity[cell.lines, drug]),
  #                       metric="auc_recomputed", dataset="GDSC1000", drug=drug, feature=feature, type=type)
}
