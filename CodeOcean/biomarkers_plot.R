#args <- commandArgs(trailingOnly=TRUE)
#datasets <- unlist(strsplit(as.character(args[1]), ","))
#library(Hmisc)
#load("../data/G2P_clarified_hits.RData")
load(file.path(path.results, "G2P_clarified_hits.RData"))
known.biomarkers <- G2P_clarified_hits
known.biomarkers[,"direction"] <- "positive"
known.biomarkers[which(known.biomarkers[, "metacindex"] < .5 | known.biomarkers[, "metamCI"] < .5), "direction"] <- "negative"
ii <- which(known.biomarkers[, "metacindex"] <= .5 &
              known.biomarkers[, "metacindex_lower"] <= .5 &
              known.biomarkers[, "metacindex_upper"] <= .5)
if(length(ii) > 0){
  known.biomarkers[ii, "metacindex"] <- abs(known.biomarkers[ii, "metacindex"] - .5) + .5
  known.biomarkers[ii, "metacindex_lower"] <- abs(known.biomarkers[ii, "metacindex_lower"] - .5) + .5
  known.biomarkers[ii, "metacindex_upper"] <- abs(known.biomarkers[ii, "metacindex_upper"] - .5) + .5

  known.biomarkers[ii, grep("cindex$", colnames(known.biomarkers))] <- abs(known.biomarkers[ii, grep("cindex$", colnames(known.biomarkers))] - .5) + .5
  tt <- abs(known.biomarkers[ii, grep("cindex_lower$", colnames(known.biomarkers))] - .5) + .5
  known.biomarkers[ii, grep("cindex_lower$", colnames(known.biomarkers))] <- abs(known.biomarkers[ii, grep("cindex_upper$", colnames(known.biomarkers))] - .5) + .5
  known.biomarkers[ii, grep("cindex_upper$", colnames(known.biomarkers))] <- tt
}
ii <- which(known.biomarkers[, "metamCI"] <= .5 &
              known.biomarkers[, "metamCI_lower"] <= .5 &
              known.biomarkers[, "metamCI_upper"] <= .5)
if(length(ii) > 0){
  known.biomarkers[ii, "metamCI"] <- abs(known.biomarkers[ii, "metamCI"] - .5) + .5
  known.biomarkers[ii, "metamCI_lower"] <- abs(known.biomarkers[ii, "metamCI_lower"] - .5) + .5
  known.biomarkers[ii, "metamCI_upper"] <- abs(known.biomarkers[ii, "metamCI_upper"] - .5) + .5

  known.biomarkers[ii, grep("mCI$", colnames(known.biomarkers))] <- abs(known.biomarkers[ii, grep("mCI$", colnames(known.biomarkers))] - .5) + .5
  tt <- abs(known.biomarkers[ii, grep("mCI_lower$", colnames(known.biomarkers))] - .5) + .5
  known.biomarkers[ii, grep("mCI_lower$", colnames(known.biomarkers))] <- abs(known.biomarkers[ii, grep("mCI_upper$", colnames(known.biomarkers))] - .5) + .5
  known.biomarkers[ii, grep("mCI_upper$", colnames(known.biomarkers))] <- tt
}
colnames(known.biomarkers)[which(colnames(known.biomarkers)=="compound")] <- "drug"
if(!"type" %in% colnames(known.biomarkers)){
  colnames(known.biomarkers)[which(colnames(known.biomarkers)=="ttype")] <- "type"
}
mycol1 <- RColorBrewer::brewer.pal(name="Dark2", n=5)
p <- c(rep(16, length(datasets)), 17)
pdf(sprintf("%s/biomarkers_CI_mCI_%s.pdf", path.results, ttype), height=4, width=8)
par(mfrow=c(2, 4))
for(i in 1:nrow(known.biomarkers)){
  biomarker <- as.vector(known.biomarkers[i, ])
  drug <- biomarker[1, "drug"]
  gene <- biomarker[1, "gene"]
  type <- biomarker[1, "type"]
  mCI <- unlist(biomarker[1, grep("mCI$", colnames(known.biomarkers))])
  CI <- unlist(biomarker[1, grep("cindex$", colnames(known.biomarkers))])
  names(mCI) <- names(CI) <- gsub("_$", "", gsub("mCI", "", names(mCI)))

  plot(CI, mCI, xlim=c(.25, 1), ylim=c(.25, 1), col=mycol1[match(names(mCI), c(datasets, "meta"))], main=sprintf("%s: %s (%s)", drug, gene, type), cex.main=.8, pch=p)
  abline(0, 1, lty=2, col="lightgray")
}
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("topright", legend=c(datasets, "Meta analysis"), pch=rep(20, length(datasets)+1), col=mycol1[1:(length(datasets)+1)], bty="n", cex=.8)
dev.off()
library(survcomp)

pdf(sprintf("%s/ci_mci_%s.pdf", path.results, ttype), height=6, width=8)
par(mar=c(4.1, 5.1, 4.1, 14), xpd=TRUE)
plot(NA, pch='', ylab='rCI', xlab='CI', xlim=c(0.4, 1), ylim=c(0.4, 1), cex.axis=1.3, cex.lab=1.8)
abline(0, 1, lty=2, xpd=FALSE)
abline(h=.5, lty=2, col="red", xpd=FALSE)
abline(v=.5, lty=2, col="red", xpd=FALSE)


gene_drug <- NULL
meta_analysis <- NULL
for(i in 1:nrow(known.biomarkers)){#which(apply(known.biomarkers, 1, function(x){length(which(is.na(x[grep("_mCI$", names(x))])))}) != length(datasets))){
  pch_i <- 19
  if(known.biomarkers[i, "direction"] == "negative"){
    pch_i <- 17
  }
  points(known.biomarkers[i, "metacindex"], known.biomarkers[i, "metamCI"], pch=pch_i, col=mycol[known.biomarkers[i, "drug"]])
  #mci_lty=ifelse(meta_mci_pvalue < 0.05, 1, 2)
  #cindex_lty=ifelse(meta_cindex_pvalue < 0.05, 1, 2)
  mci_lty <- 1
  cindex_lty <- 1

  arrows(known.biomarkers[i, "metacindex_lower"],
         known.biomarkers[i, "metamCI"],
         known.biomarkers[i, "metacindex_upper"],
         known.biomarkers[i, "metamCI"], length=0.05, angle=90, code=3, lty=mci_lty, col=mycol[known.biomarkers[i, "drug"]])
  arrows(known.biomarkers[i, "metacindex"],
         known.biomarkers[i, "metamCI_lower"],
         known.biomarkers[i, "metacindex"],
         known.biomarkers[i, "metamCI_upper"], length=0.05, angle=90, code=3, lty=cindex_lty, col=mycol[known.biomarkers[i, "drug"]])
  gene_drug <- c(gene_drug,
                sprintf("%s, %s",#: %s",
                        known.biomarkers[i, "drug"],
                        known.biomarkers[i, "gene"]
                        #capitalize(as.character(known.biomarkers[i, "type"])))
                ))

}
legend("topright", inset=c(-.65,0), legend=gene_drug, pch=19, bty="n", col=mycol[known.biomarkers[,"drug"]], cex=ifelse(nrow(known.biomarkers)>30,.3, 1.5))
dev.off()

library(forestplot)
j <- 1
for(i in 1:nrow(known.biomarkers)){#which(apply(known.biomarkers, 1, function(x){length(which(is.na(x[grep("mCI$", names(x))])))}) != length(datasets))){
  gene_drug_ss <- gene_drug[j]
  gene_drug_ss <- gsub(", ", "_", gene_drug_ss)
  pdf(sprintf("%s/%s_forest_plot_cindex_%s.pdf", path.results, ttype, gene_drug_ss), width=7, height=5)
  mci <- unlist(known.biomarkers[i, grep("mCI$", colnames(known.biomarkers))]); ii <- 1:length(mci)#which(!is.na(mci)); mci <- mci[ii]
  cindex <- unlist(known.biomarkers[i, grep("cindex$", colnames(known.biomarkers))]); cindex <- cindex[ii]
  datasets_names <- sapply(strsplit(names(mci), split="_"), function(x){x[[1]]})

  mci_lower <- unlist(known.biomarkers[i, grep("mCI_lower$", colnames(known.biomarkers))]); mci_lower <- mci_lower[ii]
  cindex_lower <- unlist(known.biomarkers[i, grep("cindex_lower$", colnames(known.biomarkers))]); cindex_lower <- cindex_lower[ii]

  mci_upper <- unlist(known.biomarkers[i, grep("mCI_upper$", colnames(known.biomarkers))]); mci_upper <- mci_upper[ii]
  cindex_upper <- unlist(known.biomarkers[i, grep("cindex_upper$", colnames(known.biomarkers))]); cindex_upper <- cindex_upper[ii]

  mci_pvalue <- unlist(known.biomarkers[i, grep("mCI_pvalue$", colnames(known.biomarkers))]); mci_pvalue <- mci_pvalue[ii]
  cindex_pvalue <- unlist(known.biomarkers[i, grep("cindex_pvalue$", colnames(known.biomarkers))]); cindex_pvalue <- cindex_pvalue[ii]

  mci_pairs <- unlist(known.biomarkers[i, grep("mCI_pairs$", colnames(known.biomarkers))]); mci_pairs <- mci_pairs[ii]
  cindex_pairs <- unlist(known.biomarkers[i, grep("cindex_pairs$", colnames(known.biomarkers))]); cindex_pairs <- cindex_pairs[ii]

  r.mean <- cbind(cindex, mci)
  r.lower <- cbind(cindex_lower, mci_lower)
  r.upper <- cbind(cindex_upper, mci_upper)
  r.pval <-  cbind(sprintf("%.1e", cindex_pvalue), sprintf("%.1e", mci_pvalue))
  r.pval[which(r.pval == "0.0e+00")] <- "<1e-16"
  r.pval[which(r.pval == "NA")] <- "_"
  r.pair <- cbind(cindex, mci)

  rownames(r.mean) <- rownames(r.lower) <- rownames(r.upper) <- rownames(r.pval) <- gsub("metamCI", "Meta analysis", datasets_names)
  colnames(r.mean) <- colnames(r.lower) <- colnames(r.upper) <- colnames(r.pval) <- c("CI (p)", "rCI (p)")

  row_names <- cbind(rownames(r.pval), r.pval)
  row_names <- rbind(c("Dataset", colnames(r.pval)), row_names)
  r.mean <- rbind(NA, r.mean)
  r.lower <- rbind(NA, r.lower)
  r.upper <- rbind(NA, r.upper)

  forestplot(row_names,
             r.mean, r.lower, r.upper,
             is.summary=c(rep(FALSE, length(datasets_names)), TRUE),
             title=gene_drug[j],
             zero=c(.49, .51),
             #grid = structure(c(2^-.5, 2^.5), gp = gpar(col = "steelblue", lty=2)),
             #boxsize=0.35,
             hrzl_lines=list("2"=gpar(lty=2, columns=1:3, col = "#000044")),
             txt_gp=fpTxtGp(label=gpar(fontfamily = "", cex = 0.8, fontface=2),
                              ticks=gpar(fontfamily = "", cex=.5, fontface=1),
                              xlab=gpar(fontfamily = "", cex=1.2, fontface=2),
                              legend=gpar(fontfamily = "", cex = 1, fontface=1)),
             col=fpColors(box=RColorBrewer::brewer.pal(n=4, name="Set2"),
                          line=RColorBrewer::brewer.pal(n=4, name="Set2"),
                          summary=RColorBrewer::brewer.pal(n=4, name="Set2")),
             xlab="Concordance Index",
             xticks= c(.3, .4, .5, .6, .7, .8, .9, 1),
             new_page = F,
             legend=c("CI", "rCI"),
             legend_args = fpLegend(pos = list("topright"),
                                    title="Method",
                                    r = unit(.1, "snpc"),
                                    gp = gpar(col="#CCCCCC", lwd=1.5)))
  j <- j + 1
  dev.off()
}

