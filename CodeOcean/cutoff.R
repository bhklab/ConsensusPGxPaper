#' ---
#' title: "Estimating the cutoff"
#' author: "Zhaleh"
#' date: '2018-11-15'
#' output: pdf_document
#' ---
#'
## ----setup, include=FALSE------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)
library(PharmacoGx)
library(wCI)
source("foo.R")
xlim <- c(-.8, .8)
ci_values <- NULL

#' In this document we check the distribution of difference between AAC values between biological replicates across several pharmacogenomics dataset to estimate the reliability of the pairs
#' # GDSC
## ----gdsc, echo=FALSE----------------------------------------------------
load("~/Google Drive/PSets/GDSC.RData")
ids <- unlist(strsplit(GDSC@drug[which(rownames(GDSC@drug) == "AZD6482"), "drugid"], split="/"))

tt <- GDSC@sensitivity$profiles[grep(paste0("drugid_", ids[1]), rownames(GDSC@sensitivity$info)), "auc_recomputed"]
ss <- GDSC@sensitivity$profiles[grep(paste0("drugid_", ids[2]), rownames(GDSC@sensitivity$info)),"auc_recomputed"]

names(tt) <- GDSC@sensitivity$info[grep(paste0("drugid_", ids[1]), rownames(GDSC@sensitivity$info)), "cellid"]
names(ss) <- GDSC@sensitivity$info[grep(paste0("drugid_", ids[2]), rownames(GDSC@sensitivity$info)), "cellid"]

nn <- intersect(names(tt), names(ss))
gdsc.delta.auc <- abs(ss[nn]-tt[nn])
qq.gdsc <- quantile(gdsc.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]
ci <- wCI::paired.concordance.index(ss[nn], tt[nn], delta.pred=0, delta.obs=0, logic.operator="or")
mci <- wCI::paired.concordance.index(ss[nn], tt[nn], delta.pred=0.2, delta.obs=0.2, logic.operator="or")
kci <- wCI::paired.concordance.index.weighted.version(ss[nn], tt[nn], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                            weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
kcil <- wCI::paired.concordance.index.weighted.version(ss[nn], tt[nn], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_laplace,#weighting_gaussian,
                                                            weightingFun_pred=kernel_laplace,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
ci_values <- rbind(ci_values, c("GDSC", ci$cindex, mci$cindex, kci$cindex, kcil$cindex))
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
x <- densCols(ss[nn], tt[nn], colramp=colorRampPalette(c("black", "white")))
dens <- col2rgb(x)[1,] + 1L
col <- cols[dens]

pdf("../results/GDSC_reps.pdf", height=5, width=5)
myScatterPlot(x=ss[nn], y=tt[nn], xlab="Replicate 1 (AAC)", ylab="Replicate 2 (AAC)", main="All replicates in GDSC", pch=20, method="plain", legend.label="", xlim=c(0,1), ylim=c(0,1), col=col)
legend("topright", legend=sprintf("CI=%.2f\nkCI=%.2f\nwCI=%.2f", ci$cindex, kci$cindex, mci$cindex), bty="n")
dev.off()
all1 <- ss[nn]
all2 <- tt[nn]
d_name <- rep("GDSC", length(nn))

delta_rep_aac <- ss[nn] - tt[nn]
gdsc_delta_rep_aac <- delta_rep_aac

pdf("../results/GDSC_cutoff.pdf", height=5, width=5)
par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
hist(delta_rep_aac, breaks=100, col="gray", main="All replicates in GDSC", xlab="Delta AAC", xlim=xlim, axes=FALSE, tck=-.01)
axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
axis(2)
CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
text(x=CI_95+.1, y=25, labels=round(CI_95, digits=2), col="red", cex=.9)
abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
#abline(v=-CI_95, col="red", lty=2)
legend("topright", inset=c(-.30,0), legend=sprintf("%s replicates\n(AZD6482 in two sites)", length(nn)), bty="n", cex=.8)
dev.off()


#' # CTRPv2
#' 536 out of 544 compounds in CTRPv2 have replicates for at least one cell so I wanted to limited the list to have replicates for at least 15 cells
## ----ctrpv2, echo=FALSE--------------------------------------------------
load("~/Google Drive/PSets/CTRPv2.RData")
dd <- apply(CTRPv2@sensitivity$n, 2, function(x){length(which(x > 1))})

dd.rep <- names(which(dd > 0))
qq.ctrp <- NULL
ctrp.delta.auc <- NULL
all.rep.pairs <- NULL
for(drug in dd.rep){
  cells <- names(which(CTRPv2@sensitivity$n[,drug] > 1))
  ii <- sensitivityInfo(CTRPv2)[which(sensitivityInfo(CTRPv2)$drugid==drug & sensitivityInfo(CTRPv2)$cellid %in% cells),]
  delta.auc <- NULL
  for(cell in cells){
    xx <- sensitivityProfiles(CTRPv2)[rownames(ii)[which(ii$cellid == cell)], "auc_recomputed", drop=T]
    tmp.delta <- NULL
    for(i in 1:(length(xx)-1)){
      for(j in (i+1):length(xx)){
        all.rep.pairs <- rbind(all.rep.pairs, c(xx[i], xx[j]))
        tmp.delta <- c(abs(xx[i]-xx[j]), tmp.delta)
      }
    }
    delta.auc <- c(delta.auc, mean(tmp.delta))
  }
  qq.ctrp <- rbind(qq.ctrp, quantile(delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")])
  ctrp.delta.auc <- c(ctrp.delta.auc, delta.auc)
}

mci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0.2, delta.obs=0.2, logic.operator="or")
ci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or")
kci <- wCI::paired.concordance.index.weighted.version(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                            weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
kcil <- wCI::paired.concordance.index.weighted.version(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_laplace,#weighting_gaussian,
                                                            weightingFun_pred=kernel_laplace,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
ci_values <- rbind(ci_values, c("CTRPv2", ci$cindex, mci$cindex, kci$cindex, kcil$cindex))
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
x <- densCols(all.rep.pairs[,1], all.rep.pairs[,2], colramp=colorRampPalette(c("black", "white")))
dens <- col2rgb(x)[1,] + 1L
col <- cols[dens]



pdf("../results/CTRPv2_reps.pdf", height=5, width=5)
myScatterPlot(x=all.rep.pairs[,1], y=all.rep.pairs[,2], xlab="Replicate 1 (AAC)", ylab="Replicate 2 (AAC)", main="All replicates in CTRPv2", pch=20, method="plain", legend.label="", xlim=c(0,1), ylim=c(0,1), col=col)
legend("topright", legend=sprintf("CI=%.2f\nkCI=%.2f\nwCI=%.2f", ci$cindex, kci$cindex, mci$cindex), bty="n")
dev.off()

all1 <- c(all1, all.rep.pairs[,1])
all2 <- c(all2, all.rep.pairs[,2])
d_name <- c(d_name, rep("CTRPv2", nrow(all.rep.pairs)))


#' ## Distribution of delta AAC values across all replicates in CTRPv2
## ----ctrpv2_rep----------------------------------------------------------
pdf("../results/CTRPv2_cutoff.pdf", height=5, width=5)
par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
delta_rep_aac <- all.rep.pairs[,1] - all.rep.pairs[,2]
ctrpv2_delta_rep_aac <- all.rep.pairs[,1] - all.rep.pairs[,2]
hist(delta_rep_aac, breaks=100, col="gray", main="All replicates in CTRPv2", xlab="Delta AAC", xlim=xlim, axes=FALSE, tck=-.01)
axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
axis(2)
CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
text(x=CI_95+.1, y=800, labels=round(CI_95, digits=2), col="red", cex=.9)
abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
#abline(v=-CI_95, col="red", lty=2)
legend("topright", inset=c(-.30,0), legend=sprintf("%s replicates\n(across %s drugs)", nrow(all.rep.pairs), length(dd.rep)), bty="n")
dev.off()


#' ## Threshold distribution across drugs in CTRPv2
## ----ctrpv2_drugs--------------------------------------------------------
qq.ctrp.avg <- apply(qq.ctrp, 2, mean)
qq.ctrp.max <- apply(qq.ctrp, 2, max)
ctrp.drug.max.delta.auc <- apply(qq.ctrp, 2, which.max)
ctrp.quantile.over.combined.drugs.delta.auc <- quantile(ctrp.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]

hist(qq.ctrp, col="gray", breaks=100, xlab="delta AAC", main="CTRPv2")
legend("topright", legend=paste(c("90%","95%"), round(ctrp.quantile.over.combined.drugs.delta.auc, digits=2), sep=":"), bty="n")


#' # GRAY
## ----gray, echo=FALSE----------------------------------------------------
load("~/Google Drive/PSets/GRAY_kallisto.RData")
dd <- apply(GRAY@sensitivity$n, 2, function(x){length(which(x > 1))})
dd.rep <- names(which(dd > 0))
qq.gray <- NULL
gray.delta.auc <- NULL
all.rep.pairs <- NULL
for(drug in dd.rep){
  cells <- names(which(GRAY@sensitivity$n[,drug] > 1))
  ii <- sensitivityInfo(GRAY)[which(sensitivityInfo(GRAY)$drugid==drug & sensitivityInfo(GRAY)$cellid %in% cells),]
  delta.auc <- NULL
  for(cell in cells){
    xx <- sensitivityProfiles(GRAY)[rownames(ii)[which(ii$cellid == cell)], "auc_recomputed", drop=T]
    tmp.delta <- NULL
    for(i in 1:(length(xx)-1)){
      for(j in (i+1):length(xx)){
        all.rep.pairs <- rbind(all.rep.pairs, c(xx[i], xx[j]))
        tmp.delta <- c(abs( xx[i]-xx[j]), tmp.delta)
      }
    }
    delta.auc <- c(delta.auc, mean(tmp.delta))
  }
  qq.gray <- rbind(qq.gray, quantile(delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")])
  gray.delta.auc <- c(gray.delta.auc, delta.auc)
}
library(wCI)
mci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0.2, delta.obs=0.2, logic.operator="or")
ci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or")
kci <- wCI::paired.concordance.index.weighted.version(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                            weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
kcil <- wCI::paired.concordance.index.weighted.version(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_laplace,#weighting_gaussian,
                                                            weightingFun_pred=kernel_laplace,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
ci_values <- rbind(ci_values, c("GRAY", ci$cindex, mci$cindex, kci$cindex, kcil$cindex))
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
x <- densCols(all.rep.pairs[,1], all.rep.pairs[,2], colramp=colorRampPalette(c("black", "white")))
dens <- col2rgb(x)[1,] + 1L
col <- cols[dens]

pdf("../results/GRAY_reps.pdf", height=5, width=5)
myScatterPlot(x=all.rep.pairs[,1], y=all.rep.pairs[,2], xlab="Replicate 1 (AAC)", ylab="Replicate 2 (AAC)", main="All replicates in GRAY", pch=20, method="plain", legend.label="", xlim=c(0,1), ylim=c(0,1), col=col)
legend("topright", legend=sprintf("CI=%.2f\nkCI=%.2f\nwCI=%.2f", ci$cindex, kci$cindex, mci$cindex), bty="n")
dev.off()

all1 <- c(all1, all.rep.pairs[,1])
all2 <- c(all2, all.rep.pairs[,2])
d_name <- c(d_name, rep("GRAY", nrow(all.rep.pairs)))


#' ## Distribution of delta AAC values across all replicates in GRAY
## ----gray_rep------------------------------------------------------------
pdf("../results/gray_cutoff.pdf", height=5, width=5)
par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
delta_rep_aac <- all.rep.pairs[,1] - all.rep.pairs[,2]
gray_delta_rep_aac <- delta_rep_aac
hist(delta_rep_aac, breaks=100, col="gray", main="All replicates in GRAY", xlab="Delta AAC", xlim=xlim, axes=FALSE, tck=-.01)
axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
axis(2)
CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
text(x=CI_95+.1, y=800, labels=round(CI_95, digits=2), col="red", cex=.9)
abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
#abline(v=-CI_95, col="red", lty=2)
legend("topright", inset=c(-.30,0), legend=sprintf("%s replicates\n(across %s drugs)", nrow(all.rep.pairs), length(dd.rep)), bty="n")
dev.off()

#' ## Threshold distribution across drugs in GRAY
## ----gray_drugs----------------------------------------------------------

qq.gray.avg <- apply(qq.gray, 2, mean)
qq.gray.max <- apply(qq.gray, 2, max)
gray.drug.max.delta.auc <- apply(qq.gray, 2, which.max)
gray.quantile.over.combined.drugs.delta.auc <- quantile(gray.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]

hist(qq.gray, col="gray", breaks=100, xlab="delta AAC", main="GRAY")
legend("topright", legend=paste(c("90%","95%"), round(gray.quantile.over.combined.drugs.delta.auc, digits=2), sep=":"), bty="n")

#' # gCSI
## ----gCSI, echo=FALSE----------------------------------------------------
load("~/Google Drive/PSets/gCSI_V1.RData")
dd <- apply(gCSI@sensitivity$n, 2, function(x){length(which(x > 1))})
dd.rep <- names(which(dd > 0))
qq.gcsi <- NULL
gcsi.delta.auc <- NULL
all.rep.pairs <- NULL
t1 <- t2 <- 0
for(drug in dd.rep){
  cells <- names(which(gCSI@sensitivity$n[,drug] > 1))
  t1 <- t1 + length(cells)
  ii <- sensitivityInfo(gCSI)[which(sensitivityInfo(gCSI)$drugid==drug & sensitivityInfo(gCSI)$cellid %in% cells),]
  t2 <- t2 + nrow(ii)
  delta.auc <- NULL
  for(cell in cells){
    xx <- sensitivityProfiles(gCSI)[rownames(ii)[which(ii$cellid == cell)], "auc_recomputed", drop=T]
    tmp.delta <- NULL
    for(i in 1:(length(xx)-1)){
      for(j in (i+1):length(xx)){
        all.rep.pairs <- rbind(all.rep.pairs, c(xx[i], xx[j]))
        tmp.delta <- c(abs( xx[i]-xx[j]), tmp.delta)
      }
    }
    delta.auc <- c(delta.auc, mean(tmp.delta))
  }
  qq.gcsi <- rbind(qq.gcsi, quantile(delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")])
  gcsi.delta.auc <- c(gcsi.delta.auc, delta.auc)
}
library(wCI)
mci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0.2, delta.obs=0.2, logic.operator="or")
ci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or")
kci <- wCI::paired.concordance.index.weighted.version(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or",                                                         weightingFun_obs=kernel_laplace,#weighting_gaussian,
weightingFun_pred=kernel_laplace,#weighting_gaussian,
alternative="greater", CPP=TRUE)
kcil <- wCI::paired.concordance.index.weighted.version(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_laplace,#weighting_gaussian,
                                                            weightingFun_pred=kernel_laplace,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
ci_values <- rbind(ci_values, c("gCSI", ci$cindex, mci$cindex, kci$cindex, kcil$cindex))
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
x <- densCols(all.rep.pairs[,1], all.rep.pairs[,2], colramp=colorRampPalette(c("black", "white")))
dens <- col2rgb(x)[1,] + 1L
col <- cols[dens]

pdf("../results/gCSI_reps.pdf", height=5, width=5)
myScatterPlot(x=all.rep.pairs[,1], y=all.rep.pairs[,2], xlab="Replicate 1 (AAC)", ylab="Replicate 2 (AAC)", main="All replicates in gCSI", pch=20, method="plain", legend.label="", xlim=c(0,1), ylim=c(0,1), col=col)
legend("topright", legend=sprintf("CI=%.2f\nkCI=%.2f\nwCI=%.2f", ci$cindex, kci$cindex, mci$cindex), bty="n")
dev.off()

all1 <- c(all1, all.rep.pairs[,1])
all2 <- c(all2, all.rep.pairs[,2])
d_name <- c(d_name, rep("gCSI", nrow(all.rep.pairs)))

#' ## Distribution of delta AAC values across all replicates in gCSI
## ----gcsi_rep------------------------------------------------------------
pdf("../results/gCSI_cutoff.pdf", height=5, width=5)
par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
delta_rep_aac <- all.rep.pairs[,1] - all.rep.pairs[,2]
gcsi_delta_rep_aac <- delta_rep_aac
hist(delta_rep_aac, breaks=100, col="gray", main="All replicates in gCSI", xlab="Delta AAC", xlim=xlim, axes=FALSE, tck=-.01)
axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
axis(2)
CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
text(x=CI_95+.1, y=250, labels=round(CI_95, digits=2), col="red", cex=.9)
abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
legend("topright", inset=c(-.30,0), legend=sprintf("%s replicates\n(across %s drugs)", nrow(all.rep.pairs), length(dd.rep)), bty="n")
dev.off()




#' ##  Threshold distribution across drugs in gCSI
## ----gcsi_drugs----------------------------------------------------------
qq.gcsi.avg <- apply(qq.gcsi, 2, mean)
qq.gcsi.max <- apply(qq.gcsi, 2, max)
gcsi.drug.max.delta.auc <- apply(qq.gcsi, 2, which.max)
gcsi.quantile.over.combined.drugs.delta.auc <- quantile(gcsi.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]

hist(qq.gcsi, col="gray", breaks=100, xlab="delta AAC", main="gCSI")
legend("topright", legend=paste(c("90%","95%"), round(gcsi.quantile.over.combined.drugs.delta.auc, digits=2), sep=":"), bty="n")

#' # AZ
## ----AZ, echo=FALSE----------------------------------------------------
load("~/Google Drive/PSets/AZ.RData")
dd <- apply(AZ@sensitivity$n, 2, function(x){length(which(x > 1))})
dd.rep <- names(which(dd > 0))
qq.AZ <- NULL
AZ.delta.auc <- NULL
all.rep.pairs <- NULL
for(drug in dd.rep){
  cells <- names(which(AZ@sensitivity$n[,drug] > 1))
  ii <- sensitivityInfo(AZ)[which(sensitivityInfo(AZ)$drugid==drug & sensitivityInfo(AZ)$cellid %in% cells),]
  delta.auc <- NULL
  for(cell in cells){
    xx <- sensitivityProfiles(AZ)[rownames(ii)[which(ii$cellid == cell)], "auc_recomputed", drop=T]
    tmp.delta <- NULL
    for(i in 1:(length(xx)-1)){
      for(j in (i+1):length(xx)){
        all.rep.pairs <- rbind(all.rep.pairs, c(xx[i], xx[j]))
        tmp.delta <- c(abs( xx[i]-xx[j]), tmp.delta)
      }
    }
    delta.auc <- c(delta.auc, mean(tmp.delta))
  }
  qq.AZ <- rbind(qq.AZ, quantile(delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")])
  AZ.delta.auc <- c(AZ.delta.auc, delta.auc)
}
set.seed(42)
rand_id <- sample(1:nrow(all.rep.pairs), size=10000)
all.rep.pairs <- all.rep.pairs[rand_id,]
library(wCI)
mci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0.2, delta.obs=0.2, logic.operator="or")
ci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or")
kci <- wCI::paired.concordance.index.weighted.version(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                      weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                      weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                      alternative="greater", CPP=TRUE)
kcil <- wCI::paired.concordance.index.weighted.version(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                       weightingFun_obs=kernel_laplace,#weighting_gaussian,
                                                       weightingFun_pred=kernel_laplace,#weighting_gaussian,
                                                       alternative="greater", CPP=TRUE)
ci_values <- rbind(ci_values, c("AZ", ci$cindex, mci$cindex, kci$cindex, kcil$cindex))
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
x <- densCols(all.rep.pairs[,1], all.rep.pairs[,2], colramp=colorRampPalette(c("black", "white")))
dens <- col2rgb(x)[1,] + 1L
col <- cols[dens]

pdf("../results/AZ_reps.pdf", height=5, width=5)
myScatterPlot(x=all.rep.pairs[,1], y=all.rep.pairs[,2], xlab="Replicate 1 (AAC)", ylab="Replicate 2 (AAC)", main="All replicates in AZ", pch=20, method="plain", legend.label="", xlim=c(0,1), ylim=c(0,1), col=col)
legend("topright", legend=sprintf("CI=%.2f\nkCI=%.2f\nwCI=%.2f", ci$cindex, kci$cindex, mci$cindex), bty="n")
dev.off()

all1 <- c(all1, all.rep.pairs[,1])
all2 <- c(all2, all.rep.pairs[,2])
d_name <- c(d_name, rep("AZ", nrow(all.rep.pairs)))


#' ## Distribution of delta AAC values across all replicates in AZ
## ----AZ_rep------------------------------------------------------------
pdf("../results/AZ_cutoff.pdf", height=5, width=5)
par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
delta_rep_aac <- all.rep.pairs[,1] - all.rep.pairs[,2]
AZ_delta_rep_aac <- delta_rep_aac
hist(delta_rep_aac, breaks=100, col="AZ", main="All replicates in AZ", xlab="Delta AAC", xlim=xlim, axes=FALSE, tck=-.01)
axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
axis(2)
CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
text(x=CI_95+.1, y=800, labels=round(CI_95, digits=2), col="red", cex=.9)
abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
#abline(v=-CI_95, col="red", lty=2)
legend("topright", inset=c(-.30,0), legend=sprintf("%s replicates\n(across %s drugs)", nrow(all.rep.pairs), length(dd.rep)), bty="n")
dev.off()

#' ## Threshold distribution across drugs in AZ
## ----AZ_drugs----------------------------------------------------------

qq.AZ.avg <- apply(qq.AZ, 2, mean)
qq.AZ.max <- apply(qq.AZ, 2, max)
AZ.drug.max.delta.auc <- apply(qq.AZ, 2, which.max)
AZ.quantile.over.combined.drugs.delta.auc <- quantile(AZ.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]

hist(qq.AZ, col="AZ", breaks=100, xlab="delta AAC", main="AZ")
legend("topright", legend=paste(c("90%","95%"), round(AZ.quantile.over.combined.drugs.delta.auc, digits=2), sep=":"), bty="n")


#' # Distribution of delta AAC values across all replicates in GRAY, gCSI, CTRPv2
## ----all datasets--------------------------------------------------------
all_replicates <- cbind("rep1"=all1, "rep2"=all2, "dataset"=d_name)
colnames(ci_values) <- c("dataset", "CI", "rCI", "kCI", "kCIl")
save(ci_values, all_replicates, file="../data/all_replicates.RData")
load("../data/all_replicates.RData")
devtools::install_github("bhklab/mci", force=T)

library(wCI)
mci <- wCI::paired.concordance.index(all_replicates[,1], all_replicates[,2], delta.pred=0.2, delta.obs=0.2, logic.operator="or")
ci <- wCI::paired.concordance.index(all_replicates[,1], all_replicates[,2], delta.pred=0, delta.obs=0, logic.operator="or")
kci <- wCI::paired.concordance.index.weighted.version(all_replicates[,1], all_replicates[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                            weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
kcil <- wCI::paired.concordance.index.weighted.version(all_replicates[,1], all_replicates[,2], delta.pred=0, delta.obs=0, logic.operator="or",
                                                            weightingFun_obs=kernel_laplace,#weighting_gaussian,
                                                            weightingFun_pred=kernel_laplace,#weighting_gaussian,
                                                            alternative="greater", CPP=TRUE)
ci_values <- rbind(ci_values, c("all", ci$cindex, mci$cindex, kci$cindex, kcil$cindex))

cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
x <- densCols(all_replicates[,1], all_replicates[,2], colramp=colorRampPalette(c("black", "white")))
dens <- col2rgb(x)[1,] + 1L
col <- cols[dens]

pdf("../results/all_reps.pdf", height=5, width=5)
myScatterPlot(x=all_replicates[,1], y=all_replicates[,2], xlab="Replicate 1 (AAC)", ylab="Replicate 2 (AAC)", main="All replicates", pch=20, method="plain", legend.label="", xlim=c(0,1), ylim=c(0,1), col=col)
#legend("topright", legend=sprintf("CI=%.2f\nkCI=%.2f\nwCI=%.2f", ci$cindex, kci$cindex, mci$cindex), bty="n")
legend("topright", legend=sprintf("CI=%.2f\nwCI=%.2f", ci$cindex, mci$cindex), bty="n")
dev.off()

pdf("../results/all_cutoff.pdf", height=5, width=5)
par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
delta_rep_aac <- c(gray_delta_rep_aac, ctrpv2_delta_rep_aac, gcsi_delta_rep_aac, gdsc_delta_rep_aac)
hist(delta_rep_aac, breaks=100, col="gray", main="", xlab="Delta AAC", xlim=xlim, axes=FALSE, tck=-.01)
axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
axis(2)
CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
text(x=CI_95+.1, y=2500, labels=round(CI_95, digits=2), col="red", cex=.9)
abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
legend("topright", inset=c(-.30,0), legend=sprintf("delta=%s (95%%)", round(CI_95, digits=2)), bty="n")
dev.off()



## ----fit_gaussian--------------------------------------------------------
delta_rep_aac=as.numeric(all_replicates[,1])- as.numeric(all_replicates[,2])
x <- delta_rep_aac
gaussian_fit <- function(par)
{
 -sum(dnorm(x, par[1], par[2], log=TRUE))
}

laplace_fit <- function(par)
{
    m <- par[1]
    b <- par[2]
    xhat <- (1/(2*b)) * exp(-abs(x - m) / b)
    -sum(log(xhat))
    #sum((r - rhat) ^ 2)
}
model_gaussian_fit <- optim(c(0, 1), gaussian_fit)    #
model_laplace <- optim(c(0, 1), laplace_fit)

library(mclust)
model_gaussian_mclust_em_fit <- Mclust(x, G=1)

hist(x, breaks=100, prob=T, xlab="delta AAC", main="Fitted distributions")
lines(density(x), lwd=1)
rr$parameters
#points(x, dnorm(x, model_gaussian_mclust_em_fit$parameters$mean, sqrt(model_gaussian_mclust_em_fit$parameters$variance$sigmasq)), col="green", pch=".")
points(x, dnorm(x, model_gaussian_fit$par[1], model_gaussian_fit$par[2]), col="blue", pch=".")
points(x, (1/(2*model_laplace$par[2])) * exp(-abs(r - model_laplace$par[1]) / model_laplace$par[2]), col="red", pch=".")#, add=T) # optim fit
legend("topright", col=c("red", "blue"), legend=c("laplace", "gaussian"), bty="n", pch=20)
library(mclust)
Mclust(x)
#Apply function nls

#' # GR AOC values distribution (gCSI)
## ----gCSI_GR, echo=FALSE-------------------------------------------------
load("~/Google Drive/PSets/gCSI_V1.RData")
dd <- apply(gCSI@sensitivity$n, 2, function(x){length(which(x > 1))})
dd.rep <- names(which(dd > 0))
qq.gcsi <- NULL
gcsi.delta.auc <- NULL
all.rep.pairs <- NULL
for(drug in dd.rep){
  cells <- names(which(gCSI@sensitivity$n[,drug] > 1))
  ii <- sensitivityInfo(gCSI)[which(sensitivityInfo(gCSI)$drugid==drug & sensitivityInfo(gCSI)$cellid %in% cells),]
  delta.auc <- NULL
  for(cell in cells){
    xx <- sensitivityProfiles(gCSI)[rownames(ii)[which(ii$cellid == cell)], "GR_AOC", drop=T]
    xx <- xx[which(!is.na(xx))]
    if(length(xx) > 1){
      tmp.delta <- NULL
    for(i in 1:(length(xx)-1)){
      for(j in (i+1):length(xx)){
        all.rep.pairs <- rbind(all.rep.pairs, c(xx[i], xx[j]))
        tmp.delta <- c(abs( xx[i]-xx[j]), tmp.delta)
      }
    }
    delta.auc <- c(delta.auc, mean(tmp.delta))
    }
  }
  qq.gcsi <- rbind(qq.gcsi, quantile(delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")])
  gcsi.delta.auc <- c(gcsi.delta.auc, delta.auc)
}
library(wCI)
mci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0.2, delta.obs=0.2, logic.operator="or")
ci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or")

myScatterPlot(x=all.rep.pairs[,1], y=all.rep.pairs[,2], xlab=sprintf("Replicate 1 (%s)", expression('GR'['AOC'])), ylab=sprintf("Replicate 2 (%s)", expression('GR'['AOC'])), main="All replicates in gCSI", pch=16, method="transparent")
legend("topleft", legend=sprintf("wCI=%.3f\nCI=%.3f", mci$cindex, ci$cindex), bty="n")


#' ## Distribution of delta GR_AOC values across all replicates in gCSI
## ----gcsi_GR_rep---------------------------------------------------------


pdf("../results/gCSI_GR_AOC.pdf", height=5, width=5)
par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
delta_rep_aac <- all.rep.pairs[,1] - all.rep.pairs[,2]
hist(delta_rep_aac, breaks=100, col="gray", main="All replicates in gCSI", xlab=expression('Delta GR'['AOC']), xlim=xlim, axes=FALSE, tck=-.01)
axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
axis(2)
CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
text(x=CI_95+.1, y=250, labels=round(CI_95, digits=2), col="red", cex=.9)
abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
legend("topright", inset=c(-.30,0), legend=sprintf("%s replicates\n(across %s drugs)", nrow(all.rep.pairs), length(dd.rep)), bty="n")
dev.off()


#'
#' # GReff(0.5 uM) values distribution (gCSI)
## ----gCSI_GR_eff, echo=FALSE---------------------------------------------
load("~/Google Drive/PSets/gCSI_V1.RData")
dd <- apply(gCSI@sensitivity$n, 2, function(x){length(which(x > 1))})
dd.rep <- names(which(dd > 0))
qq.gcsi <- NULL
gcsi.delta.auc <- NULL
all.rep.pairs <- NULL
for(drug in dd.rep){
  cells <- names(which(gCSI@sensitivity$n[,drug] > 1))
  ii <- sensitivityInfo(gCSI)[which(sensitivityInfo(gCSI)$drugid==drug & sensitivityInfo(gCSI)$cellid %in% cells),]
  delta.auc <- NULL
  for(cell in cells){
    xx <- sensitivityProfiles(gCSI)[rownames(ii)[which(ii$cellid == cell)], "GReff", drop=T]
    xx <- xx[which(!is.na(xx))]
    if(length(xx) > 1){
      tmp.delta <- NULL
    for(i in 1:(length(xx)-1)){
      for(j in (i+1):length(xx)){
        all.rep.pairs <- rbind(all.rep.pairs, c(xx[i], xx[j]))
        tmp.delta <- c(abs( xx[i]-xx[j]), tmp.delta)
      }
    }
    delta.auc <- c(delta.auc, mean(tmp.delta))
    }
  }
  qq.gcsi <- rbind(qq.gcsi, quantile(delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")])
  gcsi.delta.auc <- c(gcsi.delta.auc, delta.auc)
}
library(wCI)
mci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0.2, delta.obs=0.2, logic.operator="or")
ci <- wCI::paired.concordance.index(all.rep.pairs[,1], all.rep.pairs[,2], delta.pred=0, delta.obs=0, logic.operator="or")

myScatterPlot(x=all.rep.pairs[,1], y=all.rep.pairs[,2], xlab="Replicate 1 (GReff)", ylab="Replicate 2 (GReff)", main="All replicates in gCSI", pch=16, method="transparent")
legend("topleft", legend=sprintf("wCI=%.3f\nCI=%.3f", mci$cindex, ci$cindex), bty="n")


#' ## Distribution of delta GReff values across all replicates in gCSI
## ----gcsi_GR_eff_rep-----------------------------------------------------


pdf("../results/gCSI_GReff.pdf", height=5, width=5)
par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
delta_rep_aac <- all.rep.pairs[,1] - all.rep.pairs[,2]
gcsi_delta_rep_aac <- delta_rep_aac
hist(delta_rep_aac, breaks=100, col="gray", main="All replicates in gCSI", xlab="Delta GReff", xlim=c(-2,2), axes=FALSE, tck=-.01)
axis(1, at=c(-2, -1.5, -1, -.5, -.2, 0, .2 , .5, 1, 1.5, 2), labels=c(-2, -1.5, -1, -.5, -.2, 0, .2 , .5, 1, 1.5, 2), las=2)
axis(2)
CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
text(x=CI_95+.2, y=250, labels=round(CI_95, digits=2), col="red", cex=.9)
abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
legend("topright", inset=c(-.30,0), legend=sprintf("%s replicates\n(across %s drugs)", nrow(all.rep.pairs), length(dd.rep)), bty="n")
dev.off()


#' ##  Threshold distribution across drugs in gCSI
## ----gcsi_GR_eff_drugs---------------------------------------------------
qq.gcsi.avg <- apply(qq.gcsi, 2, mean)
qq.gcsi.max <- apply(qq.gcsi, 2, max)
gcsi.drug.max.delta.auc <- apply(qq.gcsi, 2, which.max)
gcsi.quantile.over.combined.drugs.delta.auc <- quantile(gcsi.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]

hist(qq.gcsi, col="gray", breaks=100, xlab="delta GReff", main="gCSI")
legend("topright", legend=paste(c("90%","95%"), round(gcsi.quantile.over.combined.drugs.delta.auc, digits=2), sep=":"), bty="n")


