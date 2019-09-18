options(stringsAsFactors=FALSE)
##loading pset and extracting molecular and sensitivity data
#BiocManager::install("PharmacoGx", version = "3.8")
#devtools::install_github("bhklab/mci", force=T)
#BiocManager::install("genefu", version = "3.8")
library(mCI)
library(PharmacoGx)
library(genefu)
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
gcsi_drugs <- intersect(drugs, drugNames(gCSI))
gcsi.drug.sensitivity.greff <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=gCSI, sensitivity.measure="GReff", summary.stat="median", drugs=gcsi_drugs))
gcsi.drug.sensitivity.greff.rev <- -gcsi.drug.sensitivity.greff
gcsi.drug.sensitivity.graoc <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=gCSI, sensitivity.measure="GR_AOC", summary.stat="median", drugs=gcsi_drugs))
gcsi.drug.sensitivity.aac <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=gCSI, sensitivity.measure="auc_recomputed", summary.stat="median", drugs=gcsi_drugs))
gcsi.rnaseq <- exprs(summarizeMolecularProfiles(pSet=gCSI, mDataType="rnaseq", summary.stat="median"))
gcsi.cnv <- exprs(summarizeMolecularProfiles(pSet=gCSI, mDataType="cnv", summary.stat="median"))

##Categorize drug sensitivity values
categories_extended <- c("Cytotoxic", "Constant", "Partial Growth Inhibition", "No Growth Inhibition")

categories <- c("Cytotoxic", "Constant", "Partial GI", "No GI")
gcsi.GReff.bin <- gcsi.drug.sensitivity.greff
gcsi.GReff.bin[which(gcsi.drug.sensitivity.greff <= -.1)] <- categories[1]
gcsi.GReff.bin[which(gcsi.drug.sensitivity.greff > -.1 & gcsi.drug.sensitivity.greff <= .2)] <- categories[2]
gcsi.GReff.bin[which(gcsi.drug.sensitivity.greff > .2 & gcsi.drug.sensitivity.greff <= .7)] <- categories[3]
gcsi.GReff.bin[which(gcsi.drug.sensitivity.greff > .7)] <- categories[4]


xx <- matrix(0, ncol=length(categories), nrow=length(gcsi_drugs), dimnames=list(gcsi_drugs, categories))
for(i in 1:length(gcsi_drugs)){
  tt <- table(gcsi.GReff.bin[,i])
  if(length(tt) !=0){
    xx[i, names(tt)] <- tt
  }
}
xx <- xx[-which(rowSums(xx)==0),]
xx <- apply(xx, 1, function(x){x/sum(x)})
mycol <- RColorBrewer::brewer.pal(n=4, name="Dark2")
oo <- order(xx["Cytotoxic",], decreasing=T)
xx <- xx[,oo]
nn <- colnames(xx); colnames(xx) <-NULL

pdf("../results/gCSI_categories_percentage.pdf", height=5, width=7)
par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
bp <- barplot(xx, ylab="Percentage of cells",
              col=mycol, border=F)#, main=sprintf("%s Markers", toupper(type))
#  axis(side=2, labels=seq(-.1,.4,.1), at=seq(-.1,.4,.1), las=2)#, cex.axis=1.2, las=2)
text(bp, par("usr")[3], labels=nn, srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=.9)
legend("topright", inset=c(-.35,0), legend=categories_extended, col=mycol, pch=16, bty="n", cex=.85)
dev.off()
pdf("../results/dist_gCSI_drug_response.pdf", height=24, width=8)
par(mfrow=c(6, 2))
mycol <- RColorBrewer::brewer.pal(n=4, name="Dark2")

for(drug in nn){
  if(!all(is.na(gcsi.drug.sensitivity.graoc[,drug]))){
    hist(gcsi.drug.sensitivity.aac[,drug], breaks=50, xlab="AAC", main=drug, col=mycol[1], xlim=c(0,1), border=F)
    hist(gcsi.drug.sensitivity.graoc[,drug], breaks=50, xlab=expression('GR'['AOC']), main=drug, col=mycol[2], xlim=c(-1, 2), border=F)
  }
}
dev.off()

##Acording to Supplementary Figure 11 this GR_AOC shows a stronger signal for this biomarker
cell.lines <- intersect(colnames(gcsi.rnaseq), rownames(gcsi.drug.sensitivity.greff))

predictions <- function(type, id, cell.lines){
  if(type=="cnv"){
    return(as.numeric(gcsi.cnv[id ,cell.lines]))
  }
  else if(type=="rnaseq"){
    return(as.numeric(gcsi.rnaseq[id ,cell.lines]))
  }
}

observations <- function(metric, cell.lines, drug){
  switch(metric,
         "GReff_cat"= factor(gcsi.GReff.bin[cell.lines, drug], levels=categories),
         "GReff"= gcsi.drug.sensitivity.greff.rev[cell.lines, drug],
         "GR_AOC"= gcsi.drug.sensitivity.graoc[cell.lines, drug],
         "AAC"= gcsi.drug.sensitivity.aac[cell.lines, drug]
  )
}

delta <- function(metric){
  switch(metric,
         "GReff_cat"= 0,
         "GReff"= .47,
         "GR_AOC"= .3,
         "AAC"= .2
  )
}
alternative <- function(metric){
  switch(metric,
         "GReff_cat"= "two.sided",
         "GReff"= "two.sided",
         "GR_AOC"= "two.sided",
         "AAC"= "two.sided"
  )
}

types <- c("rnaseq", "cnv")
metrics <- c("GReff_cat", "GReff", "GR_AOC", "AAC")
for(type in types){
  associations <- read.csv("../data/G2P_clarified.csv")
  associations <- associations[,-c(1, 4, 5)]
  res <- NULL
  associations <- associations[which(associations$compound %in% colnames(gcsi.drug.sensitivity.aac)),]
  ids <- rownames(featureInfo(gCSI, type))[match(associations$gene, featureInfo(gCSI, type)[, "Symbol"])]
  associations <- associations[-which(is.na(ids)),]
  ids <- ids[-which(is.na(ids))]

  for(i in 1:nrow(associations)){
    id <- ids[i]
    drug <- associations[i,"compound"]
    ci_res <- NULL
    for(metric in metrics){
      rci <- mCI::paired.concordance.index(predictions=predictions(type=type, id=id ,cell.lines=cell.lines),
                                                   observations=observations(metric, cell.lines, drug),
                                                   delta.pred=0,
                                                   delta.obs=delta(metric),
                                                   alternative=alternative(metric),
                                                   logic.operator="and", outx=T)[c("cindex", "p.value", "relevant.pairs.no", "lower", "upper", "sterr")]
      if(metric=="GReff_cat"){
        rci$cindex <- 1 - rci$cindex
      }
      names(rci) <- sprintf("%s_%s",metric, c("cindex", "p.value", "relevant.pairs.no", "lower", "upper", "sterr"))
      ci_res <- c(ci_res, rci)
    }
    res <- rbind(res, c(associations[i,c("gene", "compound")], ci_res))
  }
  colnames(res)[1:2] <- c("feature", "compound")
  associtaions_gcsi <- res[,c("feature", "compound",
                              grep("cindex", colnames(res), value=T),
                              grep("p.value", colnames(res), value=T),
                              grep("relevant.pairs.no", colnames(res), value=T))]
  associtaions_gcsi <- as.data.frame(associtaions_gcsi, stringsAsFactors=FALSE)

  associtaions_gcsi$GReff_cat_cindex <- sapply(as.numeric(associtaions_gcsi$GReff_cat_cindex), round, digits=2)
  associtaions_gcsi$GReff_cindex <- sapply(as.numeric(associtaions_gcsi$GReff_cindex), round, digits=2)
  associtaions_gcsi$GR_AOC_cindex <- sapply(as.numeric(associtaions_gcsi$GR_AOC_cindex), round, digits=2)
  associtaions_gcsi$AAC_cindex <- sapply(as.numeric(associtaions_gcsi$AAC_cindex), round, digits=2)
  associtaions_gcsi <- associtaions_gcsi[order(as.numeric(associtaions_gcsi$GReff_cat_p.value)),]

  write.csv(as.matrix(associtaions_gcsi), file=sprintf("../results/gcsi_associations_%s.csv", type))
  save(res, file=sprintf("../results/gcsi_associations_%s.RData", type))



#  library(ggrepel)
#  ggplot(results, aes(CI, -log10(pvalue))) +
#    geom_point(aes(col=sig)) +
#    scale_color_manual(values=c("black", "orange"))+
#    theme(axis.text.x=element_text(angle=90, hjust=1))+
#    geom_text_repel(data=filter(results, padj<0.05), aes(label=Gene))

}

mycol <- RColorBrewer::brewer.pal(n=4, name="Dark2")
pvl <- list()
for(type in types){
  associations <- read.csv(sprintf("../results/GReff_bonferroni/gcsi_associations_%s.csv", type))
  ids <- rownames(featureInfo(gCSI, type))[match(associations$feature, featureInfo(gCSI, type)[, "Symbol"])]
  #pdf(sprintf("../results/GReff_cat_%s.pdf", type), height=20, width=20)
  #par(mfrow=c(4, 4))
  pv <- matrix(NA, ncol=8, nrow=nrow(associations))
  colnames(pv) <- c("feature", "compound", "rCI", "rCI_pvalue", "NO_vs_ALL", "NO_vs_PAR_CONST", "NO_PAR_vs_CONST_CYTO", "PAR_CONST_vs_CYTO")
  pv[,"feature"] <- associations[,"feature"]
  pv[,"compound"] <- associations[,"compound"]
  pv[,"rCI"] <- associations[,"GReff_cat_cindex"]
  pv[,"rCI_pvalue"] <- associations[,"GReff_cat_p.value"]
  for(i in 1:nrow(associations)){
      id <- ids[i]
      drug <- associations$compound[i]

      tt <- as.data.frame(cbind(predictions(type, id ,cell.lines),
                                as.character(observations("GReff_cat", cell.lines, drug))))
      icx <- complete.cases(tt)
      tt <- tt[icx,]
      if(nrow(tt)!=0){
        m <- max(table(tt[,2]))
        xx <- matrix(NA, nrow=length(categories), ncol=m)
        rownames(xx) <- categories
        for(cat in categories){
          w=which(tt[,2]==cat)
          if(length(w)!=0){
            xx[cat, 1:length(w)] <- as.numeric(tt[w, 1])
            }
        }
        ylim <- c((min(as.numeric(xx), na.rm=T)-1), (max(as.numeric(xx), na.rm=T)+1))
        pdf(sprintf("../results/GReff/GReff_cat_%s_%s_%s.pdf", type, associations[i,"feature"], associations[i,"compound"]), height=4, width=4)
        #par(mfrow=c(4, 4))
        pv[i, c("NO_vs_ALL", "NO_vs_PAR_CONST", "NO_PAR_vs_CONST_CYTO", "PAR_CONST_vs_CYTO")] <-
          c(ifelse(!all(is.na(xx[4,])) && !all(is.na(xx[1:3,])), wilcox.test(xx[4,], xx[1:3,])$p.value, NA),
            ifelse(!all(is.na(xx[4,])) && !all(is.na(xx[2:3,])), wilcox.test(xx[4,], xx[2:3,])$p.value, NA),
            ifelse(!all(is.na(xx[3:4,]))  && !all(is.na(xx[1:2,])), wilcox.test(xx[3:4,], xx[1:2,])$p.value, NA),
            ifelse(!all(is.na(xx[2:3,])) && !all(is.na(xx[1,])), wilcox.test(xx[2:3,], xx[1,])$p.value, NA))
        invisible(genefu::boxplotplus2(xx,
                                       pt.pch = 20,
                                       .ylim=ylim,
                                       pt.col=mycol,
                                       #names=colnames(tt),#c(paste0("none (",length(res.nosignifseeds),")"), paste0("all (",length(res.allsignifseeds),")")),
                                       ylab=type,
                                       #xlab="signif seeds (# of targets)",
                                       main=sprintf("%s, %s", associations[i,"feature"], associations[i,"compound"]),
                                       .las=2,
                                       pt.cex=1))
        text(1:4, rep(ylim[2], 4), sprintf("n=%s", apply(xx, 1, function(x){length(which(!is.na(x)))})), cex=.7)
        legend("bottomleft", legend=sprintf("rCI=%.2f\np=%.2e", associations$GReff_cat_cindex[i], associations$GReff_cat_p.value[i]), bty="n", cex=.9)
        dev.off()
      }
  }
  #dev.off()
  pvl[[type]] <- pv
  write.csv(pv, file=sprintf("../results/comparison_%s.csv", type))
  pdf(sprintf("../results/volcano_%s.pdf", type), width=20, height=5)
  par(mfrow=c(1, 4))
  for(metric in metrics){
    p <- as.numeric(associations[,sprintf("%s_p.value", metric)])
    p[which(p==0)] <- 1e-24
    q <- p.adjust(p, method="bonferroni")
    cols <- rep("black", length(associations[,sprintf("%s_cindex", metric)]))
    cols[q < .05] <- "orange"
    lab <- which(-log10(q) > 5)

    plot(associations[,sprintf("%s_cindex", metric)], -log10(p),
         pch=20, col=cols, xlab="CI", main=metric, ylim=c(0,25), xlim=c(0,1))
    #abline(h=-log10(0.05), lty=2, col="gray")
    abline(v=.5, lty=2, col="gray")
    pos <- 2
    t <- -log10(p)[lab]
    # for(l in 2:length(lab)){
    #   if(t[l] < (t[l-1] +1)){
    #     if(pos[l-1]==2){
    #       pos <- c(pos, 3)
    #     }else{
    #       pos <- c(pos, 2)
    #     }
    #   }else{
    #     pos <- 2
    #   }
    # }
    pos=sample(2:4, length(lab), replace=T)
    text(associations[lab, sprintf("%s_cindex", metric)], -log10(p)[lab], labels=sprintf("%s,\n%s", associations$compound[lab], associations$feature[lab]), cex=.4, pos=pos)

  }
  dev.off()
}
# results=as.data.frame(cbind("CI"=associations[,sprintf("%s_cindex", metric)],
#                             "pvalue"=p,
#                             "padj"=q,
#                             "Gene"=sprintf("%s,\n%s", associations[,sprintf("%s_cindex", compound)], associations[,sprintf("%s_feature", metric)]),
#                             "sig"=ifelse(q<0.05, "Sig", "Not Sig")))
# results$pvalue <- as.numeric(results$pvalue)

#kruskal.test(as.numeric(gcsi.rnaseq[id ,cell.lines]),
#             factor(gcsi.GReff.bin[cell.lines, drug], levels=categories))
#boxplot(as.numeric(gcsi.rnaseq[id ,cell.lines])~factor(gcsi.GReff.bin[cell.lines, drug], levels=categories), las=2)
##delta is extraxted from biological replicates across all drugs


