#source("~/Documents/mCI/R/paired.concordance.index.R")
#devtools::install_github("bhklab/wci", force=T)
library(wCI)

source("foo.R")
library(PharmacoGx)
path.result <- "../results_Cello"
datasets <- c("AZ", "GDSC2", "gCSI", "GRAY", "CTRPv2", "CCLE")
sensitivity_metric <- c("auc_recomputed", "AAC", "auc_recomputed", "auc_recomputed", "auc_recomputed", "auc_recomputed")
sensitivity_factor <- c(1, 100, 1, 1, 100, 1)
names(sensitivity_metric) <- names(sensitivity_factor) <- datasets
source("pset_loading.R")

ll <- length(datasets)
drugs_all <- drugs <- NULL
test <- list()
j <- NULL
cc <- combn(ll, 3)
for(i in 1:ncol(cc)){
  j <- cc[,i]
  drugs <- union(drugs, intersectList(sapply(j, function(x){colnames(drug.sensitivity(datasets[x]))})))
  drugs_all <- union(drugs, intersectList(sapply(j, function(x){drugNames(pset(datasets[x]))})))
  test[[i]] <- intersectList(sapply(j, function(x){colnames(drug.sensitivity(datasets[x]))}))
  names(test)[i] <- paste(datasets[j], collapse="_")
}
setdiff(drugs_all, drugs)
mci.discordance.rate <- matrix(NA, ncol=length(datasets), nrow=length(datasets), dimnames=list(datasets, datasets))
cindex.discordance.rate <- wmci.discordance.rate <- mci.discordance.rate
wmci <- mci <- cindex <-
  wmci.pairs  <- mci.pairs <- cindex.pairs <-
  common.celllines <-
  wmci.discordant <- mci.discordant <- cindex.discordant <- list()
for(i in 1:ll){
  wmci.discordant[[datasets[i]]] <- mci.discordant[[datasets[i]]] <- cindex.discordant[[datasets[i]]] <-
    wmci[[datasets[i]]]  <- mci[[datasets[i]]] <- cindex[[datasets[i]]] <-
    wmci.pairs[[datasets[i]]]  <- mci.pairs[[datasets[i]]] <- cindex.pairs[[datasets[i]]] <-
    common.celllines[[datasets[i]]] <- list()
}
for(i1 in 1:(ll)){
  i2 <- i1
  ds1 <- datasets[i1]
  while(i2 < ll){
    i2 <- i2 + 1
    ds2 <- datasets[i2]
    #dd <- intersectList(colnames(drug.sensitivity(ds1)), colnames(drug.sensitivity(ds2)))
    dd <- intersectList(colnames(drug.sensitivity(ds1)), colnames(drug.sensitivity(ds2)), drugs)
    cell.lines <- intersect(rownames(drug.sensitivity(ds1)), rownames(drug.sensitivity(ds2)))
    #mci.temp <- matrix(NA, ncol=length(drugs), nrow=length(drugs), dimnames=list(drugs, drugs))
    mci.temp <- matrix(NA, ncol=length(dd), nrow=length(dd), dimnames=list(dd, dd))
    pairs.mci.temp <- pairs.cindex.temp <- pairs.wmci.temp <-
      cindex.temp <- wmci.temp <- mci.temp
    drug.sens1 <- drug.sensitivity(ds1)
    drug.sens2 <- drug.sensitivity(ds2)
    for(i in dd){
      #j <- i
        for(j in dd){
            print(sprintf("%s, %s, %s, %s, %s out of %s", ds1, ds2, i, j , which(dd==i), length(dd)))
            x <- drug.sens1[cell.lines,i]
            y <- drug.sens2[cell.lines,j]
            cc <- which(!is.na(x) & !is.na(y))
            if(length(cc) > 0){
              x <- x[cc]; y <- y[cc]
              tt <- wCI::paired.concordance.index(predictions=x,
                                                  observations=y,
                                                  delta.pred=0.2,
                                                  delta.obs=0.2,
                                                  alternative="greater",
                                                  logic.operator="or", CPP=TRUE)
              if(typeof(tt$cindex)=="list"){tt$cindex <- NA}
              mci.temp[i, j] <- tt$cindex
              pairs.mci.temp[i, j] <- tt$relevant.pairs.no
              tt <- wCI::paired.concordance.index(predictions=x,
                                                  observations=y,
                                                  delta.pred=0,
                                                  delta.obs=0,
                                                  alternative="greater", CPP=TRUE)
              if(typeof(tt$cindex)=="list"){tt$cindex <- NA}
              cindex.temp[i, j] <- tt$cindex
              pairs.cindex.temp[i, j] <- tt$relevant.pairs.no
              tt <- wCI::paired.concordance.index.weighted.version(predictions=x,
                                                                   observations=y,
                                                                   delta.pred=0,
                                                                   delta.obs=0,
                                                                   weightingFun_obs=wCI:::kernel_gaussian,#weighting_gaussian,
                                                                   weightingFun_pred=wCI:::kernel_gaussian,#weighting_gaussian,
                                                                   alternative="greater",
                                                                   logic.operator="or", CPP=TRUE)
              wmci.temp[i, j] <- tt$cindex
              pairs.wmci.temp[i, j] <- tt$relevant.pairs.no
            }

        }
    }
    if(nrow(mci.temp) >= 2){

    }
    plot_ci_heatmap(file_name=sprintf(file.path(path.result, "MCI_%s_%s.pdf"), ds1, ds2), data=mci.temp)



    mci.discordant[[ds1]][[ds2]] <- dd[which(!(1:length(dd) == apply(mci.temp, 1, which.max)))]
    mci.discordant[[ds2]][[ds1]] <- dd[which(!(1:length(dd) == apply(mci.temp, 2, which.max)))]
    mci.discordance.rate[ds1, ds2] <- length(mci.discordant[[ds1]][[ds2]])/length(dd)
    mci.discordance.rate[ds2, ds1] <- length(mci.discordant[[ds2]][[ds1]])/length(dd)

    plot_ci_heatmap(file_name=sprintf(file.path(path.result, "cindex_%s_%s.pdf"), ds1, ds2), data=cindex.temp)

    cindex.discordant[[ds1]][[ds2]] <- dd[which(!(1:length(dd) == apply(cindex.temp, 1, which.max)))]
    cindex.discordant[[ds2]][[ds1]] <- dd[which(!(1:length(dd) == apply(cindex.temp, 2, which.max)))]
    cindex.discordance.rate[ds1, ds2] <- length(cindex.discordant[[ds1]][[ds2]])/length(dd)
    cindex.discordance.rate[ds2, ds1] <- length(cindex.discordant[[ds2]][[ds1]])/length(dd)

    wmci.discordant[[ds1]][[ds2]] <- dd[which(!(1:length(dd) == apply(wmci.temp, 1, which.max)))]
    wmci.discordant[[ds2]][[ds1]] <- dd[which(!(1:length(dd) == apply(wmci.temp, 2, which.max)))]
    wmci.discordance.rate[ds1, ds2] <- length(wmci.discordant[[ds1]][[ds2]])/length(dd)
    wmci.discordance.rate[ds2, ds1] <- length(wmci.discordant[[ds2]][[ds1]])/length(dd)

    mci[[ds1]][[ds2]] <- mci[[ds2]][[ds1]] <- mci.temp
    wmci[[ds1]][[ds2]] <- wmci[[ds2]][[ds1]] <- wmci.temp
    cindex[[ds1]][[ds2]] <- cindex[[ds2]][[ds1]] <- cindex.temp
    mci.pairs[[ds1]][[ds2]] <- mci.pairs[[ds2]][[ds1]] <- pairs.mci.temp
    wmci.pairs[[ds1]][[ds2]] <- wmci.pairs[[ds2]][[ds1]] <- pairs.wmci.temp
    cindex.pairs[[ds1]][[ds2]] <- cindex.pairs[[ds2]][[ds1]] <- pairs.cindex.temp
    common.celllines[[ds1]][[ds2]] <- common.celllines[[ds2]][[ds1]] <- cell.lines

  }
}
save(common.celllines,
     mci,
     mci.pairs,
     mci.discordant,
     mci.discordance.rate,
     cindex,
     cindex.pairs,
     cindex.discordant,
     cindex.discordance.rate,
     wmci,
     wmci.pairs,
     wmci.discordant,
     wmci.discordance.rate,
     file=file.path(path.result, "mci_results.RData"))
table(mci.discordance.rate <= cindex.discordance.rate)


load(file.path(path.result, "mci_results.RData"))
#pairwise barplots of datasets comparing cindex and wCI
all.datasets.ci <- NULL
for(ii in 1:(length(datasets) - 1))
{
  ds1 <- datasets[ii]
  for(j in (ii + 1):length(datasets)){
    ds2 <- datasets[j]
    ff <- data.frame("CI"=diag(cindex[[ds1]][[ds2]]), "rCI"=diag(mci[[ds1]][[ds2]]), "kCI"=diag(wmci[[ds1]][[ds2]]))
    ff <- ff[which(!is.na(ff[,1])),]
    xx <- matrix(NA, nrow=1, ncol=nrow(ff) * 3)
    k <- 0
    for (i in 1:nrow(ff)){
      k <- k+1
      xx[1, c(k, k+1, k+2)] <- c(ff[i, 1], ff[i, 2], ff[i, 3])
      k <-k+2
    }
    mybarplot(Filename = file.path(path.result, sprintf("%s_%s_barplot", ds1, ds2)),
              data=xx, barsNo=ncol(ff), groupNo=nrow(ff), group.labels=rownames(ff),
              ylab.label="cindex", legend.lables=colnames(ff), cex=.7, yaxis="Regular",
              main.label=sprintf("%s, %s", ds1, ds2),
              hh=.8,
              plot.highlight=paste0(sprintf("%s out of %s drugs\n have rCI >= 0.80", length(which(ff[,2] >= .8)), nrow(ff)),
                                    sprintf("\n%s out of %s drugs\n have kCI >= 0.80", length(which(ff[,3] >= .8)), nrow(ff))))
    dd <- intersect(drugs, names(diag(cindex.pairs[[ds1]][[ds2]])))
    rr <- (diag(cindex.pairs[[ds1]][[ds2]])[dd]- diag(mci.pairs[[ds1]][[ds2]])[dd])/diag(cindex.pairs[[ds1]][[ds2]])[dd]
    all.datasets.ci <- rbind(all.datasets.ci, cbind(diag(cindex[[ds1]][[ds2]])[dd], "CI", names(diag(cindex[[ds1]][[ds2]])[dd]), sprintf("%s_%s", ds1, ds2), NA),
                                              cbind(diag(mci[[ds1]][[ds2]])[dd], "rCI", names(diag(mci[[ds1]][[ds2]])[dd]), sprintf("%s_%s", ds1, ds2), rr),
                                              cbind(diag(wmci[[ds1]][[ds2]])[dd], "kCI", names(diag(wmci[[ds1]][[ds2]])[dd]), sprintf("%s_%s", ds1, ds2), NA))


  }
}
broad_effect_drugs <- NULL
for(drug in drugs){
  i <- 1
  if(drug %in% ccle_broad_effect){
    i <- i + 1
  }
  if(drug %in% gcsi_broad_effect){
    i <- i + 1
  }
  if(drug %in% gdsc1000_broad_effect){
    i <- i + 1
  }
  if(drug %in% ctrpv2_broad_effect){
    i <- i + 1
  }
  if(i >=3){
    broad_effect_drugs <- c(broad_effect_drugs, drug)
  }
}

nn <- names(which(all.datasets.ci[,2]=="CI"))
ci <- as.numeric(all.datasets.ci[which(all.datasets.ci[,2]=="CI"), 1])
mci <- as.numeric(all.datasets.ci[which(all.datasets.ci[,2]=="rCI"), 1])
kci <- as.numeric(all.datasets.ci[which(all.datasets.ci[,2]=="kCI"), 1])
mci_p <- wilcox.test(mci, ci, alternative="greater")[["p.value"]]#1.410961e-27
kci_p <- wilcox.test(kci, ci, alternative="greater")[["p.value"]]#7.089215e-05

round(median(mci, na.rm=T)-median(ci, na.rm=T), 2)#0.18
round(median(kci, na.rm=T)-median(ci, na.rm=T), 2)#0.05


colnames(all.datasets.ci) <- c("cindex", "metric", "Drugs", "dataset", "pairs_lost")
rownames(all.datasets.ci) <- 1:nrow(all.datasets.ci)
all.datasets.ci <- as.data.frame(all.datasets.ci, stringsAsFactors=FALSE)

all.datasets.ci$cindex <- as.numeric(all.datasets.ci$cindex)
all.datasets.ci$metric <- factor(all.datasets.ci$metric)
all.datasets.ci$Drugs <- factor(all.datasets.ci$Drugs)
all.datasets.ci$dataset <- factor(all.datasets.ci$dataset)
xx <- which(all.datasets.ci$metric == "rCI"); yy <- which(all.datasets.ci$metric == "CI"); zz <- which(all.datasets.ci$metric == "kCI")
all.datasets.ci[xx, "drugs.value"] <- as.numeric(all.datasets.ci$Drugs)[xx] #+ runif(length(xx), 0.1, 0.3)
all.datasets.ci[yy, "drugs.value"] <- as.numeric(all.datasets.ci$Drugs)[yy] #+ runif(length(xx), -0.3, -.1)
all.datasets.ci[xx, "delta.cindex"] <- all.datasets.ci[yy, "delta.cindex"] <- all.datasets.ci[xx, "cindex"] - all.datasets.ci[yy, "cindex"]
all.datasets.ci[, "delta.alpha"] <- sapply(all.datasets.ci$delta.cindex, function(x){
  if(!is.na(x)){
    if(x < 0){return(0)}
    else if((x >= 0) & (x < .1)) {return(.1)}
    else if((x >= .1) & (x < .2)) {return(.4)}
    else if((x >= .2) & (x < .3)) {return(.6)}
    else if((x >= .3) & (x < .4)) {return(.8)}
    else {return(1)}
  }else{
    return(NA)
  }
})
all.datasets.ci[, "delta.line.type"] <- sapply(all.datasets.ci$delta.cindex, function(x){
  if(!is.na(x)){
    if(x < 0){return(4)}
    else if((x >= 0) & (x < .1)) {return(3)}
    else if((x >= .1) & (x < .2)) {return(3)}
    else if((x >= .2) & (x < .3)) {return(2)}
    else if((x >= .3) & (x < .4)) {return(2)}
    else {return(1)}
  }else{
    return(NA)
  }
})
all.datasets.ci$delta.line.type <- as.factor(all.datasets.ci$delta.line.type)

#load("data/18_drugs")
mm3 <- mm2 <- mm1 <- mm <- NULL
for(drug in drugs){
  d <- NULL
  if(drug %in% names(ccle_drug_response_var)){
    d <- c(d, ccle_drug_response_var[drug])
  }
  if(drug %in% names(gdsc_drug_response_var)){
    d <- c(d, gdsc_drug_response_var[drug])
  }
  if(drug %in% names(ctrpv2_drug_response_var)){
    d <- c(d, ctrpv2_drug_response_var[drug])
  }
  if(drug %in% names(gcsi_drug_response_var)){
    d <- c(d, gcsi_drug_response_var[drug])
  }
  if(drug %in% names(gcsi_drug_response_var)){
    d <- c(d, AZ_drug_response_var[drug])
  }
  mm <- c(mm, max(as.numeric(d), na.rm=T))
  mm1 <- c(mm1, median(all.datasets.ci[which(all.datasets.ci$Drugs==drug & all.datasets.ci$metric=="rCI"), "cindex"]-
                         all.datasets.ci[which(all.datasets.ci$Drugs==drug & all.datasets.ci$metric=="CI"), "cindex"]))
  mm2 <- c(mm2, min(all.datasets.ci[which(all.datasets.ci$Drugs==drug & all.datasets.ci$metric=="rCI"), "cindex"]-
                      all.datasets.ci[which(all.datasets.ci$Drugs==drug & all.datasets.ci$metric=="CI"), "cindex"]))
  mm3 <- c(mm3, max(all.datasets.ci[which(all.datasets.ci$Drugs==drug & all.datasets.ci$metric=="rCI"), "cindex"]-
                      all.datasets.ci[which(all.datasets.ci$Drugs==drug & all.datasets.ci$metric=="CI"), "cindex"]))
}
drugs_ordered <- drugs[order(mm)]
xx <- cbind("drugs"=drugs, "rCI"=mm1, "min"=mm2, "max"=mm3, "AAC variance"=mm, "col"="black")
xx <- xx[order(mm),]

xx[which(xx[,"drugs"] %in% chemo),"col"] <- RColorBrewer::brewer.pal(7,"Set1")[3]
pdf(sprintf("%s/rCI_var.pdf", path.result), height=6, width=6)
plot(NA, pch='', ylab='Delta CI', xlab='AAC variance', xlim=c(0, .06), ylim=c(-.5, .5))
for(i in 1:nrow(xx)){#which(apply(known.biomarkers, 1, function(x){length(which(is.na(x[grep("_mCI$", names(x))])))}) != length(datasets))){
  points(as.numeric(xx[i, "AAC variance"]), as.numeric(xx[i, "rCI"]), pch=19, col=xx[i, "col"])
  abline(h=0, col="gray", lty=2)
  #mci_lty=ifelse(meta_mci_pvalue < 0.05, 1, 2)
  #cindex_lty=ifelse(meta_cindex_pvalue < 0.05, 1, 2)
  mci_lty <- 1
  cindex_lty <- 1

  arrows(as.numeric(xx[i, "AAC variance"]),as.numeric(xx[i, "min"]),
         as.numeric(xx[i, "AAC variance"]),
                    as.numeric(xx[i, "max"])
                               , length=0.05, angle=90, code=3, lty=mci_lty, col=xx[i, "col"])
}
dev.off()
cor.test(abs(as.numeric(xx[,"max"])-as.numeric(xx[,"min"])), as.numeric(xx[,"AAC variance"]), method = "pearson", alternative="less")
#Pearson's product-moment correlation
#data:  abs(as.numeric(xx[, "max"]) - as.numeric(xx[, "min"])) and as.numeric(xx[, "AAC variance"])
#t = 1.0508, df = 18, p-value = 0.8464
#alternative hypothesis: true correlation is less than 0
#95 percent confidence interval:
#-1.0000000  0.5677111
#sample estimates:
#cor
#0.2404044
mm <- NULL
for(drug in drugs){
  #d <- all.datasets.ci[which(all.datasets.ci$Drugs==drug & all.datasets.ci$metric=="mCI"), "cindex"]
  d <- all.datasets.ci[which(all.datasets.ci$Drugs==drug & all.datasets.ci$metric=="rCI"), "pairs_lost"]
  mm <- c(mm, median(as.numeric(d), na.rm=T))
}
drugs_ordered2 <- drugs[order(mm, decreasing=T)]

broad_effect_drugs <- NULL
for(drug in drugs_ordered){
  i <- 1
  if(drug %in% ccle_broad_effect){
    i <- i + 1
  }
  if(drug %in% gcsi_broad_effect){
    i <- i + 1
  }
  if(drug %in% gdsc_broad_effect){
    i <- i + 1
  }
  if(drug %in% ctrpv2_broad_effect){
    i <- i + 1
  }
  if(drug %in% AZ_broad_effect){
    i <- i + 1
  }
  if(i >=3){
    broad_effect_drugs <- c(broad_effect_drugs, drug)
  }
}
#broad_effect_drugs <- intersectList(drugs_ordered, unionList(ccle_broad_effect, gcsi_broad_effect, gdsc1000_broad_effect, ctrpv2_broad_effect))
narrow_effect_drugs <- setdiff(drugs_ordered, broad_effect_drugs)
for(drugs_type in c("narrow", "broad", "all1", "all2", "all")){
  if(drugs_type == "narrow"){
    drug_effect <- NULL
    for(drug in narrow_effect_drugs){
      drug_effect <- rbind(drug_effect, all.datasets.ci[(all.datasets.ci$Drugs==drug), ])
    }
    #drug_effect <- all.datasets.ci[(all.datasets.ci$Drugs %in% narrow_effect_drugs), ]
    drug_effect$Drugs <- factor(drug_effect$Drugs, levels=narrow_effect_drugs, ordered = TRUE)
    w=10; h=5;t="Narrow effect therapies"
  }else if (drugs_type == "broad"){
    drug_effect <- NULL
    for(drug in broad_effect_drugs){
      drug_effect <- rbind(drug_effect, all.datasets.ci[(all.datasets.ci$Drugs==drug), ])
    }
    #drug_effect <- all.datasets.ci[(all.datasets.ci$Drugs %in% broad_effect_drugs), ]
    drug_effect$Drugs <- factor(drug_effect$Drugs, levels=broad_effect_drugs, ordered = TRUE)
    w=8; h=5;t="Broad effect therapies"
  }else if (drugs_type == "all1"){
    drug_effect <- NULL
    for(drug in drugs_ordered[1:ceil(length(drugs_ordered)/2)]){
      drug_effect <- rbind(drug_effect, all.datasets.ci[(all.datasets.ci$Drugs==drug), ])
    }
    #drug_effect <- all.datasets.ci[(all.datasets.ci$Drugs %in% broad_effect_drugs), ]
    drug_effect$Drugs <- factor(drug_effect$Drugs, levels=drugs_ordered[1:ceil(length(drugs_ordered)/2)], ordered = TRUE)
    w=10; h=5;t=""
  }else if (drugs_type == "all2"){
    drug_effect <- NULL
    for(drug in c(drugs_ordered[ceil(length(drugs_ordered)/2+1):length(drugs_ordered)])){
      drug_effect <- rbind(drug_effect, all.datasets.ci[(all.datasets.ci$Drugs==drug), ])
    }
    #drug_effect <- all.datasets.ci[(all.datasets.ci$Drugs %in% broad_effect_drugs), ]
    drug_effect$Drugs <- factor(drug_effect$Drugs, levels=c(drugs_ordered[ceil(length(drugs_ordered)/2+1):length(drugs_ordered)]), ordered = TRUE)
    w=10; h=5;t=""
  }else{
    drug_effect <- all.datasets.ci
    w=12; h=6;t="All compounds"
    drug_effect$Drugs <- factor(as.character(drug_effect$Drugs))
  }
  chemocol <- rep("black", length(levels(drug_effect$Drugs)))
  names(chemocol) <- levels(drug_effect$Drugs)
  chemocol[intersect(levels(drug_effect$Drugs), chemo)] <- RColorBrewer::brewer.pal(7,"Set1")[3]

  xx <- which(drug_effect$metric == "rCI"); yy <- which(drug_effect$metric == "CI"); zz <- which(drug_effect$metric == "kCI")
  drug_effect[xx, "drugs.value"] <- as.numeric(drug_effect$Drugs)[xx] #+ runif(length(xx), 0.1, 0.3)
  drug_effect[yy, "drugs.value"] <- as.numeric(drug_effect$Drugs)[yy] #+ runif(length(xx), -0.3, -.1)
  drug_effect[zz, "drugs.value"] <- as.numeric(drug_effect$Drugs)[zz] #+ runif(length(zz), 0.1, 0.3)

  #boxplot(cindex ~ metric + drugs, data=all.datasets.ci)
  #xx <- RColorBrewer::brewer.pal(n=12, name="Set3")[c(1, 4, 5, 7, 8, 12)]
  ll  <- length(levels(all.datasets.ci$dataset))
  xx <- c(RColorBrewer::brewer.pal(n=8, name="Set2"), RColorBrewer::brewer.pal(n=8, name="Dark2"))[1:ll]
  names(xx) <- levels(all.datasets.ci$dataset)
  library(dplyr)
  library(ggplot2)
  pp <- ggplot(drug_effect, aes(y=cindex, x=Drugs, fill=metric, color=metric)) +
    geom_boxplot(alpha=0.7,position=position_dodge(0.75), width=.6, outlier.shape=NA)+
    geom_hline(yintercept=0.75, linetype="dashed", color = "red")+
    geom_hline(yintercept=0.5, linetype="dashed", color = "lightgray")+
    #geom_point(aes(x=drugs.value, fill=dataset, group=metric, color=dataset), size=2, shape=21)+
    #geom_line(aes(x=drugs.value, group=interaction(dataset, Drugs), color=dataset, alpha=delta.alpha))+ #, linetype=delta.line.type didnt work properly
    scale_y_continuous(name="Concordance Index",
                       breaks=seq(40, 100, 20)/100,
                       limits=c(0.4, 1)) +
    # scale_x_discrete(name="Drugs") +
    ggtitle(t) +
    theme_classic() +
    theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(angle=45, hjust=1, size=11, colour=chemocol),
          legend.position="bottom") +
    scale_fill_manual(values=c(rCI=RColorBrewer::brewer.pal(name="Set1", n=9)[1],#rgb(209/255, 193/255, 219/255, alpha=0.5),#RColorBrewer::brewer.pal(name="PRGn", n=11)[4]
                               kCI=RColorBrewer::brewer.pal(name="PRGn", n=11)[4],
                               CI=RColorBrewer::brewer.pal(name="RdGy", n=11)[7], xx)) +
    scale_color_manual(values=c(rCI=RColorBrewer::brewer.pal(name="Set1", n=9)[1],#rgb(209/255, 193/255, 219/255, alpha=1),
                                kCI=RColorBrewer::brewer.pal(name="PRGn", n=11)[4],
                                CI=RColorBrewer::brewer.pal(name="RdGy", n=11)[8], xx))+
    labs(fill="metric")
  pp <- pp + geom_vline(xintercept=seq(1.5, length(unique(drug_effect$Drugs)) + .5, 1), col=RColorBrewer::brewer.pal(name="RdGy", n=11)[7])
  #pp
  ggsave(file=sprintf("%s/%s_cindex_mci.pdf", path.result, drugs_type), pp, height=h, width=w)

}

##drug taxonomy: how rank of drugs preserved using modified concordance index better than conventional approach
mci_drug_rank <- matrix(NA, nrow=length(drugs), ncol=2*choose(length(datasets), 2))
rownames(mci_drug_rank) <- drugs
kci_drug_rank <- ci_drug_rank <- mci_drug_rank
k <- 1
for(i1 in 1:(ll)){
  i2 <- i1
  ds1 <- datasets[i1]
  while(i2 < ll){
    i2 <- i2 + 1
    ds2 <- datasets[i2]
    dd <- intersectList(colnames(mci[[ds1]][[ds2]]), drugs)
    for(drug in dd){
      temp <- mci[[ds1]][[ds2]]
      o1 <- order(temp[drug, ], decreasing=T)
      o2 <- order(temp[ ,drug], decreasing=T)
      mci_drug_rank[drug, k] <- which(o1==which(colnames(temp)==drug))/ncol(temp)
      mci_drug_rank[drug, k+1] <- which(o2==which(rownames(temp)==drug))/ncol(temp)

      temp <- cindex[[ds1]][[ds2]]
      o1 <- order(temp[drug, ], decreasing=T)
      o2 <- order(temp[ ,drug], decreasing=T)
      ci_drug_rank[drug, k] <- which(o1==which(colnames(temp)==drug))/ncol(temp)
      ci_drug_rank[drug, k+1] <- which(o2==which(rownames(temp)==drug))/ncol(temp)

      temp <- wmci[[ds1]][[ds2]]
      o1 <- order(temp[drug, ], decreasing=T)
      o2 <- order(temp[ ,drug], decreasing=T)
      kci_drug_rank[drug, k] <- which(o1==which(colnames(temp)==drug))/ncol(temp)
      kci_drug_rank[drug, k+1] <- which(o2==which(rownames(temp)==drug))/ncol(temp)
    }
    k <- k + 2
  }
}
library(ggplot2)
library(reshape2)
library(plyr)
dd <- drugs[order(apply(mci_drug_rank, 1, median, na.rm=T))]
combined_data <- ldply(list(CI=t(ci_drug_rank[dd,]), mCI=t(mci_drug_rank[dd,]), kCI=t(kci_drug_rank[dd,])))
combined_data_melt <- melt(combined_data, id.vars = 1)
colnames(combined_data_melt) <- c("metric", "Drugs", "Rank")

#Plot and use facet_wrap for each data.frame
pp <- ggplot(combined_data_melt, aes(y=Rank, x=Drugs, fill=metric, color=metric)) +
  geom_boxplot(alpha=0.7,position=position_dodge(0.75), width=.6, outlier.shape=NA)+
  geom_hline(yintercept=0.8, linetype="dashed", color = "red")+
  #geom_point(aes(x=drugs.value, fill=dataset, group=metric, color=dataset), size=2, shape=21)+
  #geom_line(aes(x=drugs.value, group=interaction(dataset, Drugs), color=dataset, alpha=delta.alpha))+ #, linetype=delta.line.type didnt work properly
  scale_y_continuous(name="Rank",
                     breaks=seq(0, 1, .1),
                     limits=c(0, 1)) +
  # scale_x_discrete(name="Drugs") +
  ggtitle(t) +
  theme_classic() +
  theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11),
        legend.position="bottom") +
  scale_fill_manual(values=c(mCI=RColorBrewer::brewer.pal(name="Set1", n=9)[1],#rgb(209/255, 193/255, 219/255, alpha=0.5),#RColorBrewer::brewer.pal(name="PRGn", n=11)[4]
                             kCI=RColorBrewer::brewer.pal(name="PRGn", n=11)[4],
                             CI=RColorBrewer::brewer.pal(name="RdGy", n=11)[7], xx)) +
  scale_color_manual(values=c(mCI=RColorBrewer::brewer.pal(name="Set1", n=9)[1],#rgb(209/255, 193/255, 219/255, alpha=1),
                              kCI=RColorBrewer::brewer.pal(name="PRGn", n=11)[4],
                              CI=RColorBrewer::brewer.pal(name="RdGy", n=11)[8], xx))+
  labs(fill="metric") +
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(file=sprintf("../results/drug_taxonomy_%s.pdf", drugs_type), pp, height=h, width=w)

ranks <- cbind("mCI"=apply(mci_drug_rank, 1, median, na.rm=T), "CI"=apply(ci_drug_rank, 1, median, na.rm=T))
ranks <- ranks[drugs_ordered,]
pdf(file.path(path.result, "drug_taxonomy_barplot.pdf"), height=6, width=8)
par(mar=c(12, 4.1, 4.1, 4.1), xpd=TRUE)
rr <- as.numeric(ranks[,"CI"]) - as.numeric(ranks[,"mCI"])
rr[abs(rr) < 1e-3] <- 0
names(rr) <- rownames(ranks)
barplot(rr, ylab="delta Rank", las=2,
        col=RColorBrewer::brewer.pal(name="Set1", n=9)[3], main="Drug Taxonomy")
dev.off()
