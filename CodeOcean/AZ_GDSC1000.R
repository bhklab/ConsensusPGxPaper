library(PharmacoGx)
load("~/Google Drive/Psets/GDSC1000.RData")
load("~/Google Drive/Psets/AZ.RData")
common <- intersectPSet(pSets = list("AZ"=AZ, "GDSC"=GDSC1000), 
                        intersectOn = c("cell.lines", "drugs"))

sens_info <- common$AZ@sensitivity$info


rrs <- unique(sprintf("%s__%s", sens_info$cellid, sens_info$drugid))
all <- sprintf("%s__%s", sens_info$cellid, sens_info$drugid)
cc <- NULL
for(rr in rrs){
  cc <- c(cc, length(which(all==rr)))
}
AZ_data <- read.csv("../data/data_share_with_BHK/raw_data_for_release_280619.csv")

for( rr in rrs[which(cc > 20)]){
  rr_ <- strsplit(rr, split="__")[[1]]
  cellid <- rr_[1]
  drugid <- rr_[2]
  exp_az <- which(common$AZ@sensitivity$info$cellid == cellid &
                    common$AZ@sensitivity$info$drugid == drugid)
  exp_gdsc <- which(common$GDSC@sensitivity$info$cellid == cellid &
                      common$GDSC@sensitivity$info$drugid == drugid)
  AZ_data_sub <- AZ_data[which(AZ_data$DRUG_ID==drugInfo(AZ)[drugid, "astrazeneca.numerical_drugid"] &
                                AZ_data$CELL_LINE_NAME==cellInfo(AZ)[cellid, "raw_cell_line_name"]),]
  write.csv(AZ_data_sub, file=sprintf("../data/%s.csv", rr))
  pdf(file=sprintf("../data/%s.pdf", rr), height=8, width=8)
  mycol <- c(rep("blue", length(exp_az)), rep("red", length(exp_gdsc)))
  PharmacoGx::drugDoseResponseCurve(drug=drugid, plot.type="Fitted", ,legends.label="auc_recomputed",
                                    cellline=cellid, mycol=mycol,
                                    pSets=list("AZ"=common$AZ, "GDSC"=common$GDSC),
                                    summarize.replicates=FALSE)
  dev.off()
}

pdf(file="../data/AZ_GDSC1000_identical_exp.pdf", height=ceiling(length(uids)/4) * 5, width=4*5)
par(mfrow=c(ceiling(length(uids)/4), 4))

dd <- 0; i <- 0
for(rr in rrs){
  i <- i + 1
  print(i)
  rr <- strsplit(rr, split="__")[[1]]
  cellid <- rr[1]
  drugid <- rr[2]
  exp_az <- which(common$AZ@sensitivity$info$cellid == cellid &
                      common$AZ@sensitivity$info$drugid == drugid)
  exp_gdsc <- which(common$GDSC@sensitivity$info$cellid == cellid &
                common$GDSC@sensitivity$info$drugid == drugid)
  for(e1 in exp_az){
    d1 <- common$AZ@sensitivity$raw[e1,,"Dose"]
    d1 <- d1[which(!is.na(d1))]
    v1 <- common$AZ@sensitivity$raw[e1,1:length(d1),"Viability"]
    for(e2 in exp_gdsc){
      d2 <- common$GDSC@sensitivity$raw[e2,,"Dose"]
      d2 <- d1[which(!is.na(d2))]
      v2 <- common$GDSC@sensitivity$raw[e2,1:length(d2),"Viability"]
      if(length(d1) == length(d2)){
        if(all(d1 == d2)){# & all(v1 == v2)){
          print(sprintf("the same exp in AZ(%s), GDSC(%s), %s:%s", e1, e2, drugid, cellid))
          dd <- dd + 1
        }
      }
    }
  }
  mycol <- c(rep("blue", length(exp_az)), rep("red", length(exp_gdsc)))
  PharmacoGx::drugDoseResponseCurve(drug=drugid, plot.type="Fitted",legends.label="auc_recomputed",
                                    cellline=cellid, mycol=mycol,
                                    pSets=list("AZ"=common$AZ, "GDSC"=common$GDSC),
                                    summarize.replicates=FALSE)
}
print(sprintf("There are %s duplicated exprerimernts", dd))
dev.off()
##test intersection with the other datasets
load("~/Google Drive/Psets/gCSI_V1.RData")
intersect(rownames(gCSI@drug), rownames(AZ@drug))
load("~/Google Drive/Psets/CTRPv2.RData")
intersect(rownames(CTRPv2@drug), rownames(AZ@drug))


AZ_fitted <- read.csv("../data/data_share_with_BHK/fitted_data_for_release_280619.csv")
xx=sprintf("%s__%s",AZ_fitted$CELL_LINE_NAME, AZ_fitted$DRUG_ID)
uu <- unique(xx)
cc <- NULL
for(rr in uu){
  cc <- c(cc, length(which(xx==rr)))
}
rr=xx[which(cc==2)[1]]
rr_=strsplit(rr, split="__")[[1]]
cellid <- rr_[1]
drugid <- rr_[2]
cellid=rownames(AZ@cell)[which(AZ@cell$fitted_cell_line_name==cellid)]
drugid=rownames(AZ@drug)[which(AZ@drug$astrazeneca.numerical_drugid==drugid)]
exp_az <- which(AZ@sensitivity$info$cellid == cellid &
                  AZ@sensitivity$info$drugid == drugid)
AZ_data_sub <- AZ_data[which(AZ_data$DRUG_ID==drugInfo(AZ)[drugid, "astrazeneca.numerical_drugid"] &
                               AZ_data$CELL_LINE_NAME==cellInfo(AZ)[cellid, "raw_cell_line_name"]),]
write.csv(AZ_data_sub, file=sprintf("../data/%s.csv", rr))
PharmacoGx::drugDoseResponseCurve(drug=drugid, plot.type="Fitted", ,legends.label="auc_recomputed",
                                  cellline=cellid, mycol=mycol,
                                  pSets=list("AZ"=common$AZ, "GDSC"=common$GDSC),
                                  summarize.replicates=FALSE)
