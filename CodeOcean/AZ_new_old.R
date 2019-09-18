#badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
#dosage_no <- 9
options(stringsAsFactors = FALSE)

drugs <- read.csv("../data/data_share_with_BHK/compound_dictionary_280619.csv")
AZ_data <- read.csv("../data/data_share_with_BHK/raw_data_for_release_280619.csv")
AZ_fitted <- read.csv("../data/data_share_with_BHK/fitted_data_for_release_280619.csv")

max_conc <- unique(AZ_fitted$MAX_CONC_MICROMOLAR)
#[1] 10.00  3.00  1.00 30.00  2.50  5.00  6.00  0.05
curationCell <- readRDS("../data/curationCell.rds")
curationDrug <- readRDS("../data/curationDrug.rds")
for(col in colnames(curationCell)){
  curationCell[,col] <- as.character(curationCell[,col])
}
for(col in colnames(curationDrug)){
  curationDrug[,col] <- as.character(curationDrug[,col])
}

xx <- curationCell$unique.cellid
xx[which(is.na(xx))] <- curationCell$astrazeneca.cellid[which(is.na(xx))]
rownames(curationCell) <- xx

xx <- curationDrug$unique.drugid
xx[which(is.na(xx))] <- curationDrug$astrazeneca.drugid[which(is.na(xx))]
rownames(curationDrug) <- xx
xx <- setdiff(drugs$DRUG_ID, curationDrug$astrazeneca.numerical_drugid)
for(x in xx){
  ii <- which(drugs$DRUG_ID == x)
  curationDrug <- rbind(curationDrug, c(drugs$DRUG_NAME[ii], 
                                        drugs$DRUG_NAME[ii], 
                                        drugs$PUTATIVE_TARGET[ii],
                                        NA,
                                        NA, 
                                        NA, 
                                        NA,
                                        drugs$DRUG_ID[ii]))
}
rownames(curationDrug) <- curationDrug$unique.drugid


uids <- unique(sprintf("%s__%s", AZ_data$DRUG_ID,
                       AZ_data$CELL_LINE_NAME))
print(length(uids))
rr <- 1
myfn <- "../data/AZ_sensitivity.csv"
if(!exists(myfn)){
  AZ_data[,'uid'] <- NA
  AZ_data[,'dosage_no'] <- NA
  AZ_data[,'max_conc'] <- NA
  for(uid in uids){
    print(rr)
    rr <- rr + 1
    uid <- unlist(strsplit(uid, "__"))
    yy <- which(AZ_data$DRUG_ID==uid[1] &
                  AZ_data$CELL_LINE_NAME==uid[2])
    start_dose_i <- end_dose_i <- 1
    replicate_no <- 1
    while ( end_dose_i < length(yy)){
      if(AZ_data$CONC[yy[1]] == AZ_data$CONC[yy[2]]){
        xx <- 1
        while(AZ_data$CONC[yy[xx]] == AZ_data$CONC[yy[xx+1]]){xx <- xx + 2; if(xx > length(yy)){break}}
        temp <- AZ_data[yy[1:(xx-1)],]
        AZ_data[yy[1]:yy[(xx-1)/2],] <- temp[which(1:nrow(temp) %%2 ==0),]
        AZ_data[yy[(xx-1)/2+1]:yy[(xx-1)],] <- temp[which(1:nrow(temp) %%2 ==1),]
      }
      if(AZ_data$CONC[yy[end_dose_i + 1]] < AZ_data$CONC[yy[end_dose_i]]){
        max_conc <- AZ_data$CONC[yy[end_dose_i]]
        while(AZ_data$CONC[yy[end_dose_i]] > AZ_data$CONC[yy[end_dose_i+1]]){end_dose_i <- end_dose_i + 1; if(end_dose_i == length(yy)){break}}
      }else{
        while(AZ_data$CONC[yy[end_dose_i]] < AZ_data$CONC[yy[end_dose_i+1]]){end_dose_i <- end_dose_i + 1; if(end_dose_i == length(yy)){break}}
        max_conc <- AZ_data$CONC[yy[end_dose_i]]
      }
      uid_ <- sprintf("%s__%s__%s", uid[1], uid[2], replicate_no)
      AZ_data[yy[start_dose_i:end_dose_i], 'uid'] <- uid_
      AZ_data[yy[start_dose_i:end_dose_i], 'dosage_no'] <- end_dose_i - start_dose_i + 1
      AZ_data[yy[start_dose_i:end_dose_i],'max_conc'] <- max_conc
      replicate_no <- replicate_no + 1
      start_dose_i <- end_dose_i <- end_dose_i + 1
      print(uid_)
    }
    if(any(is.na(AZ_data[yy, 'uid']))){
      print(uid)
      break
    }
  }
  write.csv(AZ_data, myfn)
}else{
  read.csv(myfn)
}
uids <- unique(AZ_data[,'uid'])
length(uids)
##[1] 32619 Number of experiments
Az_raw_sensitivity <- array(NA, dim=c(length(uids), max(unique(AZ_data$dosage_no)), 2),
                           dimnames=list(uids, sprintf("dose%s",1:dosage_no), c("Dose", "Viability")))
sensitivity_info <- matrix(NA, nrow=length(uids), ncol=4 , 
                           dimnames=list(uids, c("cellid", "drugid", "dosage_no", "max_conc")))
pdf(file="../data/AZ_plots.pdf", height=ceiling(length(uids)/4) * 5, width=4*5)
par(mfrow=c(ceiling(length(uids)/4), 4))
rr <- 1
for(uid in uids){
  print(rr)
  rr <- rr + 1
  uid_ <- unlist(strsplit(uid, "__"))
  cell <- uid_[2]#gsub(pattern=badchars, replacement="", x=toupper(uid_[2]))
  drug <- uid_[1]
  cellid <- rownames(curationCell)[which(curationCell$raw_cell_line_name == cell)]
  drugid <- rownames(curationDrug)[which(curationDrug$astrazeneca.numerical_drugid == drug)]
  if(is.na(cellid) || is.na(drugid)){
    print(uid)
    break
  }
  xx <- AZ_data[which(AZ_data$uid==uid),]
  sensitivity_info[uid, ] <- c(cellid, drugid, xx$dosage_no[1], xx$max_conc[1])
  
  conc <- as.character(xx$CONC)
  viability <- as.character((xx$INTENSITY/max(xx$INTENSITY)) * 100)
  

  # PharmacoGx::drugDoseResponseCurve(concentrations=as.numeric(conc),
  #                                   viabilities=as.numeric(viability),
  #                                   conc_as_log=F, viability_as_pct=T, plot.type="Fitted", 
  #                                   title=sprintf("%s: %s", drugid, cellid))
  # auc <- PharmacoGx::computeAUC(concentration=conc,
  #                               viability=viability,
  #                               conc_as_log=T,
  #                               viability_as_pct=T)
  
  Az_raw_sensitivity[uid, 1:length(conc),'Dose'] <- conc
  Az_raw_sensitivity[uid, 1:length(conc),'Viability'] <- viability
}
dev.off()     
myfn <- "../data/AZ_recomputed.RData"
if(!exists(myfn)){
  recomputed <- PharmacoGx:::.calculateFromRaw(Az_raw_sensitivity, cap=100) 
  save(recomputed, file=myfn)
}else{
  load(myfn)
}

profiles <- cbind("auc_recomputed"=recomputed$AUC/100, "ic50_recomputed"=recomputed$IC50)    
rownames(profiles) <- rownames(Az_raw_sensitivity)
library(Biobase)
load("~/Google Drive/Psets/GDSC1000.RData")
pp <- pData(GDSC1000@molecularProfiles$rna)
samples <- rownames(pp)[which(pp$cellid %in% curationCell$unique.cellid)]
rna <- GDSC1000@molecularProfiles$rna[,samples]
xx <- which(is.na(curationCell$unique.tissueid))
curationCell$unique.tissueid[xx] <- curationCell$disease_type[xx]
curationCell$tissueid <- curationCell$unique.tissueid
curationTissue <- curationCell[,c("tissueid", "unique.tissueid")]

dd <- which(sensitivity_info[,"dosage_no"] != 2)
sensitivity_info <- sensitivity_info[dd,]
Az_raw_sensitivity <- Az_raw_sensitivity[dd,,]
profiles <- profiles[dd,]

AZ <- PharmacoSet(name="AZ",
                  molecularProfiles=list("rna"=rna),
                  cell=curationCell, 
                  drug=curationDrug, 
                  sensitivityInfo=sensitivity_info, 
                  sensitivityRaw=Az_raw_sensitivity, 
                  sensitivityProfiles=profiles, 
                  sensitivityN=NULL,  
                  curationCell=curationCell, 
                  curationDrug=curationDrug, 
                  curationTissue=curationTissue, datasetType="sensitivity")

save(AZ, file="~/Google Drive/Psets/AZ.RData")


