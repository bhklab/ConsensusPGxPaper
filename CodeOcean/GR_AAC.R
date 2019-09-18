metric1 <- "auc_recomputed"
metric2 <- "GReff"
load(file=file.path(sprintf("../results_gcsi_%s", tolower(metric1)), "G2P_clarified_hits_raw_list.RData"))
colnames(G2P_clarified_hits) <- sprintf("%s_1", colnames(G2P_clarified_hits))
G2P_clarified_hits_temp <- G2P_clarified_hits
load(file=file.path(sprintf("../results_gcsi_%s", tolower(metric2)), "G2P_clarified_hits_raw_list.RData"))
G2P_clarified_hits <- cbind(G2P_clarified_hits, G2P_clarified_hits_temp)
##pvalue correction
G2P_clarified_hits$metamCI_fdr <- G2P_clarified_hits$metamCI_pvalue
G2P_clarified_hits$metamCI_fdr_1 <- G2P_clarified_hits$metamCI_pvalue_1
#  G2P_clarified_hits$metamCI_fdr <- G2P_clarified_hits$metamCI_fdr * 3
#  G2P_clarified_hits$metacindex_fdr <- G2P_clarified_hits$metacindex_fdr * 3
G2P_clarified_hits$metamCI_fdr <- p.adjust(G2P_clarified_hits$metamCI_fdr, method="bonferroni")
G2P_clarified_hits$metacindex_fdr <- p.adjust(G2P_clarified_hits$metamCI_fdr_1, method="bonferroni")
G2P_clarified_hits <- G2P_clarified_hits[which((G2P_clarified_hits$metamCI_fdr < 0.05 &
                                                  G2P_clarified_hits$metamCI > 0.55)|
                                                 (G2P_clarified_hits$metamCI_fdr_1 < 0.05 &
                                                    G2P_clarified_hits$metamCI_1 > 0.55)|
                                                 (G2P_clarified_hits$metamCI_fdr_1 < 0.05 &
                                                    G2P_clarified_hits$metamCI_1 < 0.45)|
                                                 (G2P_clarified_hits$metamCI_fdr < 0.05 &
                                                    G2P_clarified_hits$metamCI < 0.45)),]
G2P_clarified_hits <- G2P_clarified_hits[which((G2P_clarified_hits$metamCI < 0.5 & G2P_clarified_hits$metamCI_1 < 0.5)|
                                                 (G2P_clarified_hits$metamCI >= 0.5 & G2P_clarified_hits$metamCI_1 >= 0.5)),]
G2P_clarified_hits <- G2P_clarified_hits[order(G2P_clarified_hits$compound, G2P_clarified_hits$gene),]
save(G2P_clarified_hits, file="../results/G2P_clarified_hits.RData")
save(G2P_clarified_hits, file="../results/G2P_clarified_hits_copy.RData")
dd <- unique(G2P_clarified_hits$compound)
mycol <- c(colorRamps::primary.colors(), RColorBrewer::brewer.pal(n=9, name="Set1"));
if(length(mycol) < length(dd)){
  mycol <- rep(mycol, ceiling( length(dd)/length(mycol)))
  length(mycol)
}
mycol <- mycol[1:length(dd)]
names(mycol) <- dd
drugs_types <- list()
types <- c("cnv", "mutation", "expression")

for(ttype in types){
  load("../results/G2P_clarified_hits_copy.RData")
  temp <- G2P_clarified_hits
  temp <- temp[which(!is.na(temp$metamCI) & !is.na(temp$metamCI_1)),]
  if(ttype %in% names(table(temp$ttype))){
    temp <- temp[which(temp$ttype == ttype),]
    temp[which(temp[,"metamCI_upper_1"] == 0 & temp[,"metamCI_fdr_1"] == 1), "metamCI_upper_1"] <- 1
    temp[which(temp[,"metamCI_lower_1"] == 0 & temp[,"metamCI_fdr_1"] == 1), "metamCI_lower_1"] <- .4
    G2P_clarified_hits <-
      temp[sapply(unique(temp$compound), function(x){
          min(which(temp$compound==x & temp$metamCI==max(temp[which(temp$compound == x), "metamCI"])))
      }, simplify=T),]
    print(ttype)
    print(length(unique(G2P_clarified_hits$compound)))
    print(nrow(G2P_clarified_hits))
    gain=which((G2P_clarified_hits$metacindex <= 0.55 | G2P_clarified_hits$metacindex_fdr >= 0.05) &
                 (G2P_clarified_hits$metamCI > 0.55 & G2P_clarified_hits$metamCI_fdr < 0.05))
    print(length(gain))
    print(length(gain)/nrow(G2P_clarified_hits))
    loss=which((G2P_clarified_hits$metamCI <= 0.55 | G2P_clarified_hits$metamCI_fdr  >= 0.05) &
                 (G2P_clarified_hits$metacindex > 0.55 & G2P_clarified_hits$metacindex_fdr < 0.05))
    print(length(loss))
    print(length(loss)/nrow(G2P_clarified_hits))
    #G2P_clarified_hits[loss,]
    drugs_types[[ttype]] <- G2P_clarified_hits$compound
    save(G2P_clarified_hits, file=sprintf("../results/G2P_clarified_hits_%s.RData", ttype))
    save(G2P_clarified_hits, file="../results/G2P_clarified_hits.RData")
    pdf(sprintf("../results/%s_%s_%s.pdf", metric1, metric2, ttype), height=8, width=10)
    par(mar=c(10.1, 4.1, 5.1, 13.5), xpd=TRUE)
    m1 <- metric1; m2 <- metric2
    if(metric1=="GR_AOC"){m1 <- expression('GR'['AOC'])}
    if(metric2=="GR_AOC"){m2 <-  expression('GR'['AOC'])}
    if(metric1=="auc_recomputed"){m1 <- "AAC"}
    if(metric2=="auc_recomputed"){m2 <- "AAC"}
    plot(NA, pch='', ylab=m2, xlab=m1, xlim=c(0.35, 1), ylim=c(0.35, 1), main="")#sprintf("%s based biomarkers", ttype))
    abline(0, 1, lty=2, xpd=FALSE)
    abline(h=.5, lty=2, col="red", xpd=FALSE)
    abline(v=.5, lty=2, col="red", xpd=FALSE)

    gene_drug <- NULL
    meta_analysis <- NULL
    known.biomarkers <- G2P_clarified_hits
    for(i in 1:nrow(known.biomarkers)){#which(apply(known.biomarkers, 1, function(x){length(which(is.na(x[grep("_mCI$", names(x))])))}) != length(datasets))){
      points(known.biomarkers[i, "metamCI_1"], known.biomarkers[i, "metamCI"], pch=19, col=mycol[known.biomarkers[i, "compound"]])
      #mci_lty=ifelse(meta_mci_pvalue < 0.05, 1, 2)
      #cindex_lty=ifelse(meta_cindex_pvalue < 0.05, 1, 2)
      mci_lty <- 1
      cindex_lty <- 1

      arrows(known.biomarkers[i, "metamCI_lower_1"],
             known.biomarkers[i, "metamCI"],
             known.biomarkers[i, "metamCI_upper_1"],
             known.biomarkers[i, "metamCI"], length=0.05, angle=90, code=3, lty=mci_lty, col=mycol[known.biomarkers[i, "compound"]])
      arrows(known.biomarkers[i, "metamCI_1"],
             known.biomarkers[i, "metamCI_lower"],
             known.biomarkers[i, "metamCI_1"],
             known.biomarkers[i, "metamCI_upper"], length=0.05, angle=90, code=3, lty=cindex_lty, col=mycol[known.biomarkers[i, "compound"]])
      j <- j + 1
      gene_drug <- c(gene_drug,
                     sprintf("%s, %s",
                             known.biomarkers[i, "compound"],
                             known.biomarkers[i, "gene"]))

    }
    View(known.biomarkers[,c(1:3, grep("metamCI", colnames(known.biomarkers)))])
    legend("topright", inset=c(-.40,0), legend=gene_drug, pch=19, bty="n", col=mycol, cex=ifelse(nrow(known.biomarkers)>30,.3, .8))
    dev.off()
    pdf(sprintf("../results/pairs_%s_%s_%s.pdf", metric1, metric2, ttype), height=5, width=7)
    par(mar=c(10.1, 4.1, 5.1, 4.1))
    tt <- rbind(known.biomarkers$gCSI_V1_mCI_pairs/known.biomarkers$gCSI_V1_cindex_pairs, known.biomarkers$gCSI_V1_mCI_pairs_1/known.biomarkers$gCSI_V1_cindex_pairs_1)

    print("%s, %s keeps more pairs: %s", ttype, wilcox.test(tt[1,], tt[2,], alternative = "greater"))
    bp=barplot(tt, beside=T, col=mycol[1:2], ylab="Pairs% kept", xlab="", las=2, cex.axis=.7, cex.names=.8, ylim=c(0,1), main="")#sprintf("%s based biomarkers", ttype))
    nn <- t(cbind(sprintf("%s (%s)", known.biomarkers$compound, known.biomarkers$gene), ""))
    text(bp+.4, par("usr")[3], labels=nn, srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=1.1)
    legend("topright", legend=c(m2, m1), pch=19, bty="n", col=mycol[1:2])
    dev.off()
  }
}
xx = which(abs(as.numeric(gCSI@sensitivity$profiles[,metric1])) < 2 & !is.na(gCSI@sensitivity$profiles[,metric2]))
cells <- gCSI@sensitivity$info$cellid[xx]
tt <- as.numeric(gCSI@cell[cells,"DoublingTime"])
cols <-  colorRampPalette(c("#2b8cbe", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(140)
#x <- densCols(gCSI@sensitivity$profiles[xx, metric1], gCSI@sensitivity$profiles[xx, metric2], colramp=colorRampPalette(c("black", "white")))
#dens <- col2rgb(x)[1,] + 1L
#col <- cols[dens]

#plot(gCSI@sensitivity$profiles$auc_recomputed[xx], gCSI@sensitivity$profiles$GR_AOC[xx])
myScatterPlot(sprintf("../results/%s_%s.pdf", m1 , metric2), xlab=m1, ylab=m2,
              x=gCSI@sensitivity$profiles[xx, metric1],
              y=gCSI@sensitivity$profiles[xx, metric2], method="plain", pch=20, col=cols[tt])#, xlim=range(gCSI@sensitivity$profiles$GR_AOC[xx]))

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)

  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
lut = cols <-  colorRampPalette(c("#2b8cbe", "#00FEFF", "#45FE4F",
                                  "#FCFF00", "#FF9400", "#FF3100"))(140)
#pdf("../results/color_bar_supp_10.pdf", height=7, width=5)
color.bar(lut, 0, 140, title='Cell Division Time', nticks=5)
#dev.off()
