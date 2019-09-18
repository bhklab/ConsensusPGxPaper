source("foo.R")
load("../data/all_replicates.RData")
ci_values <- cbind(ci_values, "y"=c(25, 800, 800, 250, 2500))
#devtools::install_github("bhklab/mci", force=T)
library(mCI)

mci <- mCI::paired.concordance.index(as.numeric(all_replicates[,1]), as.numeric(all_replicates[,2]), delta.pred=0.2, delta.obs=0.2, logic.operator="or")
ci <- mCI::paired.concordance.index(as.numeric(all_replicates[,1]), as.numeric(all_replicates[,2]), delta.pred=0, delta.obs=0, logic.operator="or")
kci <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[,1]), as.numeric(all_replicates[,2]), delta.pred=0, delta.obs=0, logic.operator="or",
                                                      weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                      weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                      alternative="greater", CPP=TRUE)
kcil <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[,1]), as.numeric(all_replicates[,2]), delta.pred=0, delta.obs=0, logic.operator="or",
                                                       weightingFun_obs=kernel_laplace,#weighting_gaussian,
                                                       weightingFun_pred=kernel_laplace,#weighting_gaussian,
                                                       alternative="greater", CPP=TRUE)
ci_values <- rbind(ci_values, c("all", ci$cindex, mci$cindex, kci$cindex, kcil$cindex))
#ci_values[which(ci_values[,"dataset"]=="all"),] <- c("all", 0.7821853, 0.9181906, 0.8811281, 0.8277402)
#rownames(ci_values) <- ci_values[,"dataset"]
pdf("../results/all_reps.pdf", height=5, width=5)
myScatterPlot(x=all_replicates[,1], y=all_replicates[,2], xlab="Replicate 1 (AAC)", ylab="Replicate 2 (AAC)", main="All replicates", pch=16, method="transparent", legend.label="", xlim=c(0,1), ylim=c(0,1))
#legend("topright", legend=sprintf("CI=%.2f\nkCI=%.2f\nmCI=%.2f", ci$cindex, kci$cindex, mci$cindex), bty="n")
legend("topright", legend=sprintf("CI=%.2f\nkCI=%.2f\nmCI=%.2f", .79, .97, 0.92), bty="n")
dev.off()

compute <- F
xlim <- c(-.8, .8)

cols <-  colorRampPalette(c("#2b8cbe", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
for(dataset in rownames(ci_values)){
  if(dataset != "all"){
    d <- which(all_replicates[,"dataset"]==dataset)
    main_label <- dataset
    border=NA
  }else{
    d <- 1:nrow(all_replicates)
    main_label <- ""
    border=TRUE
  }
  if(compute){
    mci <- mCI::paired.concordance.index(as.numeric(all_replicates[d,1]), as.numeric(all_replicates[d,2]), delta.pred=0.2, delta.obs=0.2, logic.operator="or")
    ci <- mCI::paired.concordance.index(as.numeric(all_replicates[d,1]), as.numeric(all_replicates[d,2]), delta.pred=0, delta.obs=0, logic.operator="or")
    kci <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[d,1]), as.numeric(all_replicates[d,2]), delta.pred=0, delta.obs=0, logic.operator="or",
                                                          weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                          weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                          alternative="greater", CPP=TRUE)
    kcil <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[d,1]), as.numeric(all_replicates[d,2]), delta.pred=0, delta.obs=0, logic.operator="or",
                                                           weightingFun_obs=kernel_laplace,#weighting_gaussian,
                                                           weightingFun_pred=kernel_laplace,#weighting_gaussian,
                                                           alternative="greater", CPP=TRUE)
    ci_values[which(ci_values[,"dataset"]==dataset),] <- c(dataset, ci$cindex, mci$cindex, kci$cindex, kcil$cindex)

  }
  x <- densCols(all_replicates[d, 1], all_replicates[d, 2], colramp=colorRampPalette(c("black", "white")))
  dens <- col2rgb(x)[1,] + 1L
  col <- cols[dens]
  cols_lab <- RColorBrewer::brewer.pal(7,"Set1")
  col_hist <- cols_lab[2]
  main_label_col <- cols_lab[2]
  if(dataset=="all"){
    col_hist <- "gray"
    main_label_col <- "black"
  }
  if(dataset=="GDSC"){
    col_hist <- cols_lab[3]
    main_label_col <- cols_lab[3]
  }

  pdf(sprintf("../results/%s_reps.pdf", dataset), height=5, width=5)
  #if(dataset=="all"){
  #  par(mar=c(4.1, 4.1, 5.1, 5.1), xpd=TRUE)
  #}
  myScatterPlot(x=as.numeric(all_replicates[d,1]),
                y=as.numeric(all_replicates[d,2]),
                xlab="Replicate 1 (AAC)",
                ylab="Replicate 2 (AAC)",
                main=main_label,
                col.main=main_label_col,
                pch=20, method="plain", legend.label="", xlim=c(0,1), ylim=c(0,1), col=col)
  if(dataset=="all"){#, inset=c(-.30,0)
    legend("topleft", legend=sprintf("CI=%.2f, kCI=%.2f, rCI=%.2f",
                                      as.numeric(ci_values[dataset, "CI"]),
                                      as.numeric(ci_values[dataset, "kCI"]),
                                      as.numeric(ci_values[dataset, "rCI"])), bty="n")
  }else{
    legend("topleft", legend=sprintf("CI=%.2f, kCI=%.2f, rCI=%.2f",
                                      as.numeric(ci_values[dataset, "CI"]),
                                      as.numeric(ci_values[dataset, "kCI"]),
                                      as.numeric(ci_values[dataset, "rCI"])), bty="n")
  }
  dev.off()

  pdf(sprintf("../results/%s_cutoff.pdf", dataset), height=5, width=5)
  par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
  delta_rep_aac <- as.numeric(all_replicates[d,1]) - as.numeric(all_replicates[d,2])
  hist(delta_rep_aac, breaks=100, col=col_hist, main=main_label, col.main=main_label_col,
       xlab=expression(paste(Delta, "AAC")), xlim=xlim, axes=FALSE, tck=-.01, border=border)
  axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
  axis(2)
  CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
  abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
  if(dataset != "all"){
    text(x=CI_95+.1, y=as.numeric(ci_values[dataset,"y"]), labels=round(CI_95, digits=2), col="red", cex=.9)
  }
  abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
  if(dataset != "all"){
    legend("topright", inset=c(-.10,0), legend=sprintf("%s replicates", length(delta_rep_aac)), bty="n")
  }
  dev.off()

}

d <- which(all_replicates[,"dataset"]==dataset)
all_replicates <- all_replicates[-d,]
dim(all.rep.pairs)
all_replicates <- rbind(all_replicates, cbind(all.rep.pairs, rep("CTRPv2", nrow(all.rep.pairs))))

system.time(
  kci <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[,1]), as.numeric(all_replicates[,2]), delta.pred=0, delta.obs=0, logic.operator="or",
                                                        weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                        weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                        alternative="greater", CPP=TRUE, permute=F)
)

devtools::install_github("bhklab/mci", force=T)
library(mCI)
library(parallel)
load("../data/all_replicates.RData")
kCI <- mclapply(1:1000, function(x){ii <- sample(1:nrow(all_replicates), 1000)
  kk <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[ii,1]), as.numeric(all_replicates[ii,2]), delta.pred=0, delta.obs=0, logic.operator="or",
                                                      weightingFun_obs=kernel_gaussian,#weighting_gaussian,
                                                      weightingFun_pred=kernel_gaussian,#weighting_gaussian,
                                                      alternative="greater", CPP=TRUE, permute=F)$cindex
  cc <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[ii,1]), as.numeric(all_replicates[ii,2]), delta.pred=0, delta.obs=0, logic.operator="or",
                                                       alternative="greater", CPP=TRUE, permute=F)$cindex
  return(c("kCI"=kk, "CI"=cc))
}, mc.cores=parallel::detectCores())
save(kCI, file="../results/kCI_down_sampled2.RData")

load("../results/kCI_down_sampled2.RData")
kCI <- unlist(kCI)
kCI_down <- kCI[which(names(kCI) == "kCI")]
CI_down <- kCI[which(names(kCI) == "CI")]
load("../results/kCI_per.RData")
mycol=RColorBrewer::brewer.pal(7, "Dark2")

pdf("../results/kCI_permutation.pdf", width=5, height=5)
par(mar=c(4.1, 4.1, 4.1, 7), xpd=TRUE)
hist(kCI_per, xlim=c(.7,.95), col=rgb(0,0,1,1/4), xlab="CI", main="", border=F)
hist(CI_down, add=T, col=rgb(1,0,0,1/4), border=F)
hist(kCI_down, add=T, col=mycol[1], border=F)
abline(v=ci_values["all", "kCI"], col="red", lty=2, xpd=F)
text(x=as.numeric(ci_values["all", "kCI"])+.02, y=190, labels=round(as.numeric(ci_values["all", "kCI"]), digits=2), col="red", cex=.9)

abline(v=ci_values["all", "CI"], col="blue", lty=2, xpd=F)
text(x=as.numeric(ci_values["all", "CI"])+.02, y=190, labels=round(as.numeric(ci_values["all", "CI"]), digits=2), col="blue", cex=.9)

legend("topright", inset=c(-.5, 0), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4), mycol[1]), legend=c("Permuted weights", "Down-sampled CI", "Down-sampled kCI"), bty="n", pch=16)
dev.off()

load("../data/all_replicates.RData")
rCI <- mclapply(1:1000, function(x){ii <- sample(1:nrow(all_replicates), 1000)
rr <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[ii,1]), as.numeric(all_replicates[ii,2]), delta.pred=.2, delta.obs=.2, logic.operator="or",
                                                     alternative="greater", CPP=TRUE, permute=F)$cindex
pp <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[ii,1]), as.numeric(all_replicates[ii,2]), delta.pred=.2, delta.obs=.2, logic.operator="or",
                                                     alternative="greater", CPP=TRUE, permute=T)$cindex
cc <- mCI::paired.concordance.index.weighted.version(as.numeric(all_replicates[ii,1]), as.numeric(all_replicates[ii,2]), delta.pred=0, delta.obs=0, logic.operator="or",
                                                     alternative="greater", CPP=TRUE, permute=F)$cindex
return(c("rCI"=rr, "rCI_per"=pp, "CI"=cc))
}, mc.cores=parallel::detectCores())
save(rCI, file="../results/rCI_down_sampled.RData")
rCI <- unlist(rCI)
rCI_down <- rCI[which(names(rCI) == "rCI")]
rCI_per <- rCI[which(names(rCI) == "rCI_per")]
CI_down <- rCI[which(names(rCI) == "CI")]

pdf("../results/rCI_permutation.pdf", width=5, height=5)
par(mar=c(4.1, 4.1, 4.1, 7), xpd=TRUE)
hist(rCI_per, xlim=c(.7,.95), col=rgb(0,0,1,1/4), xlab="CI", main="", border=F)
hist(CI_down, add=T, col=rgb(1,0,0,1/4), border=F)
hist(rCI_down, add=T, col=mycol[1], border=F)
abline(v=ci_values["all", "rCI"], col="red", lty=2, xpd=F)
text(x=as.numeric(ci_values["all", "rCI"])+.02, y=190, labels=round(as.numeric(ci_values["all", "rCI"]), digits=2), col="red", cex=.9)

abline(v=ci_values["all", "CI"], col="blue", lty=2, xpd=F)
text(x=as.numeric(ci_values["all", "CI"])+.02, y=190, labels=round(as.numeric(ci_values["all", "CI"]), digits=2), col="blue", cex=.9)

legend("topright", inset=c(-.5, 0), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4), mycol[1]), legend=c("Permuted weights", "Down-sampled CI", "Down-sampled rCI"), bty="n", pch=16)
dev.off()
