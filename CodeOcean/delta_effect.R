delta.effect <- function(first.dataset.name, second.dataset.name, first.sensitivity, second.sensitivity, drugs, logic.operator="or", mus=seq(0, 0.3, 0.01), mycol=RColorBrewer::brewer.pal(7, "Set1")[c(3, 4)]){
  library(PharmacoGx)
  all.time <- NULL
  all.time2 <- NULL
  cindex.list <- list()
  pairs.no.list <- list()
  pvalue.list <- list()
  cell.lines <- intersect(rownames(first.sensitivity), rownames(second.sensitivity))
  drugs <- PharmacoGx::intersectList(colnames(first.sensitivity), colnames(second.sensitivity), drugs)
  pdf(sprintf("../results/%s_%s_delta_effect.pdf", first.dataset.name, second.dataset.name), height=length(drugs)*4, width=10)
  par(mfrow=c(length(drugs), 2))
  for(drug in drugs){
    pp <- proc.time()[["elapsed"]]
    pred <- first.sensitivity[cell.lines, drug]
    obs <- second.sensitivity[cell.lines, drug]

    cindex.delta <- NULL
    pairs.no.delta <- NULL
    pvalue.delta <- NULL
    for(delta in mus){
      all.time <- c(all.time, system.time(ttt <- mCI::paired.concordance.index(predictions=pred,
                                                                          observations=obs,
                                                                          alternative="greater",
                                                                          delta.pred=delta,
                                                                          delta.obs=delta,
                                                                          logic.operator=logic.operator))[["elapsed"]])
      cindex.delta <- c(cindex.delta, ttt$cindex)
      pairs.no.delta <- c(pairs.no.delta, ttt$relevant.pairs.no)
      pvalue.delta <- c(pvalue.delta, ttt$p.value)
    }
    cindex.list[[drug]] <- cindex.delta
    pairs.no.list[[drug]] <- pairs.no.delta
    pvalue.list[[drug]] <- pvalue.delta

    names(pvalue.list[[drug]]) <-names(cindex.list[[drug]]) <- names(pairs.no.list[[drug]]) <- names(cindex.delta) <- names(pairs.no.delta) <- mus
    plot(mus, cindex.delta, xlab="delta", ylab="cindex", pch=20, col=mycol[1], main=sprintf("%s, mCI", drug))
    plot(mus, pairs.no.delta, xlab="delta", ylab="#pairs", pch=20, col=mycol[2], main=sprintf("%s, pairs no", drug))
    all.time2 <- c(all.time2, proc.time()[["elapsed"]] - pp)
  }
  dev.off()
  return(list("cindex"=cindex.list, "pairs.no"=pairs.no.list, "pvalue"= pvalue.list, "time"=all.time, "time2"=all.time2))
}
source("pset_loading.R")
i <- 3; j <- 4
dd <- c("Nilotinib", "PLX4720", "paclitaxel", "Docetaxel")
datasets.comparison <- list()
#load("../data/3_drugs", verbose=T)
xx <- delta.effect(first.dataset.name=datasets[i],
                   second.dataset.name=datasets[j],
                   first.sensitivity=drug.sensitivity(datasets[i]),
                   second.sensitivity=drug.sensitivity(datasets[j]),
                   drugs=dd)
ii <- sprintf("%s_%s", datasets[i], datasets[j])
datasets.comparison[[ii]] <- list()
datasets.comparison[[ii]][["cindex"]] <- xx[["cindex"]]
datasets.comparison[[ii]][["pairs.no"]] <- xx[["pairs.no"]]
datasets.comparison[[ii]][["pvalue"]] <- xx[["pvalue"]]
all.time <- rbind(all.time, xx[["time"]])

myf <- "../results/datasets_comparison_deltas.RData"
if(!exists(myf)){
  save(datasets.comparison, file=myf)
}else{
  load(myf, verbose = T)
}
i <- 3; j <- 4
ii <- sprintf("%s_%s", datasets[i], datasets[j])
first.sensitivity=drug.sensitivity(datasets[i])
second.sensitivity=drug.sensitivity(datasets[j])
cell.lines <- intersect(rownames(first.sensitivity), rownames(second.sensitivity))
#dd <- PharmacoGx::intersectList(colnames(first.sensitivity), colnames(second.sensitivity), drugs)

mycol=RColorBrewer::brewer.pal(7, "Set1")[c(2, 3, 4)]
pdf(sprintf("../results/%s_%s_delta_effect_specific_drugs2.pdf", datasets[i], datasets[j]), height=length(dd)*4, width=4.5)
par(mfrow=c(length(dd), 1))
for (drug in dd){#names(datasets.comparison[[ii]][["cindex"]])){
  mus <- as.numeric(names(datasets.comparison[[ii]][["cindex"]][[drug]]))
  par(mar = c(5,5,2,5))
  pv <- datasets.comparison[[ii]][["pvalue"]][[drug]]
  pch <- vector(length = length(pv))
  pch[which(pv < 0.05)] <- 17
  #pch[which(pv >= 0.01 & pv < 0.05)] <- "o"
  pch[which(pv>= 0.05)] <- 20
  plot(mus, datasets.comparison[[ii]][["cindex"]][[drug]], xlab="delta", ylab="rCI", pch=pch, col=mycol[2], main=drug, las=1, cex.lab=1, font.lab=2, ylim=c(.45, 1))
  #abline(v=0.2, lty=2, col="red")
  #lines(mus, datasets.comparison[[ii]][["cindex"]][[drug]], col=mycol[2], lty=2)
  par(new=T)
  pp <- datasets.comparison[[ii]][["pairs.no"]][[drug]]
  plot(mus, round(pp/pp[1], digits=2) * 100,
                  pch=20, col=mycol[3], axes=F, xlab=NA, ylab=NA, cex=1.2, ylim=c(0,100))
  test <- axis(side = 4, las=T)
  mtext(side = 4, line = 3, 'Percentage of pairs kept', cex=0.9, font=2)
  mtext(side = 4, line = 4, sprintf('Original #pairs=%s', format(pp[1], big.mark=",", scientific=FALSE)), cex=0.9, font=2)

  #legend("topright", bty="n", legend=sprintf("AUC=%s",round(ff[drug, "AUC"], digits=2)), cex=0.8, text.font=2)
  #plot(x=first.sensitivity[cell.lines, drug], y=second.sensitivity[cell.lines, drug], pch=20, col=mycol[1], xlab=datasets[i], ylab=datasets[j])
#  abline(h=0.2, lty=2, col="red")
#  abline(v=0.2, lty=2, col="red")
}
dev.off()
pdf(sprintf("results/%s_%s_delta_effect_specific_drugs3.pdf", datasets[i], datasets[j]), height=length(dd)*4, width=4)
par(mfrow=c(length(dd), 1))
for(drug in dd){
  plot(x=first.sensitivity[cell.lines, drug],
       y=second.sensitivity[cell.lines, drug],
       pch=20, col=mycol[1], xlab=sprintf("%s (AAC)", datasets[i]), ylab=sprintf("%s (AAC)", datasets[j]))
}
dev.off()

library(caTools)
xx <- datasets.comparison[[ii]]
mci <- list()
for(i in 1:length(xx[["cindex"]])){
  gg <- c(1, 1, 1, 1)
  aa <- c(0, 0, 0, 0)
  for(j in 1:length(xx[["cindex"]][[i]])){
    gg[1] <- gg[1] * (xx[["cindex"]][[i]][[j]] ^ (j))
    aa[1] <- aa[1] + (j)

    gg[2] <- gg[2] * (xx[["cindex"]][[i]][[j]])
    aa[2] <- aa[2] + 1

    gg[3] <- gg[3] + xx[["cindex"]][[i]][[j]]
    aa[3] <- aa[3] + 1

    gg[4] <- gg[4] + (xx[["cindex"]][[i]][[j]] * j)
    aa[4] <- aa[4] + j
  }
  mci[["geom.mean"]] <- c(mci[["geom.mean"]], gg[2] ^ (1/aa[2]))
  mci[["weighted.geom.mean"]] <- c(mci[["weighted.geom.mean"]], gg[1] ^ (1/aa[1]))
  mci[["arith.mean"]] <- c(mci[["arith.mean"]], gg[3]/aa[3])
  mci[["weighted.arith.mean"]] <- c(mci[["weighted.arith.mean"]], gg[4]/aa[4])
  mci[["cutoff.0.1"]] <- c(mci[["cutoff.0.1"]], xx[["cindex"]][[i]][["0.1"]])
  mci[["cutoff.0.2"]] <- c(mci[["cutoff.0.2"]], xx[["cindex"]][[i]][["0.2"]])
  mci[["AUC"]] <- c(mci[["AUC"]], trapz(x=as.numeric(names(xx[["cindex"]][[i]])), y=xx[["cindex"]][[i]]))
}
ff=sapply(mci, cbind)
rownames(ff) <- names(xx[["cindex"]])

pdf("results/mci_trend_summarization_metrics.pdf", height=5, width=8)
par(mar = c(7,4,2,1))
barcols <- RColorBrewer::brewer.pal(10, "Paired")[c(3:6, 9:10)]
barplot(t(ff[,1:6]),beside=T, col=barcols, las=2,  border=NA, main="Metrics to summarize mCI trend")
legend("topleft", col=barcols, legend=c("Geometric mean", "Wighted Geometric mean", "Arithmetic mean", "Wighted Arithmetic mean", "cut off=0.1", "cut off=0.2"), pch=15, bty="n", cex=.7,text.font=2)
dev.off()

ll <- c("Geometric mean", "Wighted Geometric mean", "Arithmetic mean", "Wighted Arithmetic mean", "cut off=0.1", "cut off=0.2", "AUC")
CI.metrics <- matrix(NA, ncol=ncol(ff), nrow=ncol(ff), dimnames=list(colnames(ff), colnames(ff)))
for(i in 1:ncol(ff)){
  for(j in 1:ncol(ff)){
    CI.metrics[i, j] <- paired.concordance.index(ff[,i], ff[,j], delta.pred=0, delta.obs=0)$cindex
    print(sprintf("%s with %s: %s", colnames(ff)[i], colnames(ff)[j], round(CI.metrics[i, j], digits=2)))
  }
}

pdf("results/mci_trend_summarization_metrics_concordance.pdf", height=8, width=8)
par(mar = c(10,10,3,1))
image(CI.metrics, col=RColorBrewer::brewer.pal(9,"BuPu"), axes = FALSE)
axis(side=2, labels=ll, at=seq(0, 1, by=1/(nrow(CI.metrics)-1)), las=T, line=NA, cex.axis=.8, font=2)
axis(side=1, labels=ll, at=seq(0, 1, by=1/(nrow(CI.metrics)-1)), las=2, line=NA, cex.axis=.8, font=2)
mtext(side=3, line=1, "Concordance of summarization metrics", font=2)
dev.off()

k <- 0
yy <- matrix(NA, nrow=1, ncol=nrow(ff) * 2)
for (i in 1:nrow(ff)){
  k <- k+1
  yy[1, c(k, k+1)] <- c(ff[i,1], ff[i,2])
  k <-k+1
}
mybarplot(Filename = file.path(path.result, "mCI_WeightedmCI"),
          data=yy, barsNo=2, groupNo=length(mci), group.labels=names(xx[["cindex"]]),
          ylab.label="cindex", legend.lables=c("mCI", "WeightedmCI"), cex=.7, yaxis="Regular",
          main.label="mCI_WeightedmCI",
          hh=.8,
          plot.highlight="")




