#badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
#dosage_no <- 9
options(stringsAsFactors = FALSE)

drugs <- read.csv("../data/data_share_with_BHK/compound_dictionary_280619.csv")
screens <- list.files("../data/AZ/Raw_Data", pattern=".csv")
screen <- screens[1]
AZ_data_1 <- read.csv(sprintf("../data/AZ/Raw_Data/%s", screen))
AZ_fitted <- read.csv("../data/data_share_with_BHK/fitted_data_for_release_280619.csv")
uids <- NULL
for (screen in screens){
  AZ_data_1 <- read.csv(sprintf("../data/AZ/Raw_Data/%s", screen))
  dd <- AZ_data_1[grep("L", AZ_data_1$TAG), ]
  tt <- unique(sprintf("%s__%s___%s__%s",
                       dd$COSMIC_ID,
                       dd$DRUG_ID,
                       dd$BARCODE,
                       strsplit(screen, "_")[[1]][1]))
  uids <- c(uids, tt)
}
cell_drug <- sapply(uids, function(x){strsplit(x, "___")[[1]][1]})
barcodes <- sapply(uids, function(x){strsplit(x, "___")[[1]][2]})

print(sprintf("Number of experiments: %s", length(cell_drug)))#14354
print(sprintf("Number of unique experiments: %s", length(unique(cell_drug))))#11193

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
rownames(curationCell) <- xx

xx <- curationDrug$unique.drugid
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


##[1] 14354 Number of experiments
Az_raw_sensitivity <- array(NA, dim=c(length(uids), 9, 2),
                           dimnames=list(uids, sprintf("dose%s",1:dosage_no), c("Dose", "Viability")))
sensitivity_info <- matrix(NA, nrow=length(uids), ncol=7 ,
                           dimnames=list(uids, c("cellid", "drugid", "dosage_no", "max_conc",
                                                 "screen_ID", "barcode", "scan_ID")))
#pdf(file="../data/AZ_plots.pdf", height=ceiling(length(uids)/4) * 5, width=4*5)
#par(mfrow=c(ceiling(length(uids)/4), 4))
counter <- 1
for (screen in screens){
  AZ_data_1 <- read.csv(sprintf("../data/AZ/Raw_Data/%s", screen))
  barcodes <- unique(AZ_data_1$BARCODE)
  for(i in 1:length(barcodes)){
    ff <- TRUE
    dd <- AZ_data_1[AZ_data_1$BARCODE == barcodes[i],]
    NC <- dd[grep("NC-1", dd$TAG), "INTENSITY"]
    if (length(NC) == 0){
      NC <- dd[grep("NC-0", dd$TAG), "INTENSITY"]
    }
    if (length(NC) == 0){
      print(sprintf("%s : %s, No negative controls exit for this plate!"), screen, barcodes[i])
      ff <- FALSE
    }
    if (length(unique(dd$CELL_LINE_NAME)) != 1){
      print(sprintf("%s : %s, More than one cell line on this plate!"), screen, barcodes[i])
      ff <- FALSE
    }
    B <- dd[grep("B", dd$TAG), "INTENSITY"]
    if (length(NC) == 0){
      print(sprintf("%s : %s, No blank well on this plate!"), screen, barcodes[i])
      B <- 0
    }
    PC <- dd[grep("PC", dd$TAG), "INTENSITY"]
    if(mean(PC) >= mean(NC)){
      print(sprintf("%s : %s, QC failed! Median intesity at high concentration is higher than no treatment!"), screen, barcodes[i])
      ff <- FALSE
    }
    if(length(NC) <= 5) {mNC <- median(NC)}else{mNC <- mean(NC)}
    if(length(B) <= 5) {mB <- median(B)}else{mB <- mean(B)}
    if(ff){
      ttd <- dd[grep("L", dd$TAG), ]
      drugs <- unique(ttd$DRUG_ID)
      cell <- unique(ttd$CELL_LINE_NAME)
      barcode <- unique(ttd$BARCODE)
      screen_ID <- unique(ttd$SCREEN_ID)
      for(drug in drugs){
        cellid <- rownames(curationCell)[which(curationCell$raw_cell_line_name == cell)]
        drugid <- rownames(curationDrug)[which(curationDrug$astrazeneca.numerical_drugid == as.character(drug))]

        tt <- ttd[ttd$DRUG_ID == drug, ]
        scan_ID <- paste(unique(tt$SCAN_ID), collapse=" ,")
        uid <- sprintf("%s__%s___%s__%s",
                       tt$COSMIC_ID[1],
                       tt$DRUG_ID[1],
                       tt$BARCODE[1],
                       strsplit(screen, "_")[[1]][1])
        cc <- sort(unique(tt$CONC))
        vv <- NULL
        for (conc in cc){
          vv <- c(vv, (median(tt$INTENSITY[tt$CONC==conc])-mB)/(mNC - mB))
        }
        Az_raw_sensitivity[uid, 1:length(cc),'Dose'] <- as.character(cc)
        Az_raw_sensitivity[uid, 1:length(cc),'Viability'] <- as.character(vv * 100)
        sensitivity_info[uid, ] <- c(cellid, drugid, length(cc), max(cc), screen_ID, barcode, scan_ID)
        print(counter)
        print(uid)
        print(sensitivity_info[uid, ])
        counter <- counter + 1
       }
    }
  }
}

##check replicates
xx <- sprintf("%s__%s", sensitivity_info[,"cellid"], sensitivity_info[,"drugid"])
ddd <- which(duplicated(xx))
length(ddd)
drugs <- unique(sensitivity_info[ddd, "drugid"])
for (drug in drugs){
  print(sprintf("%s: %s", drug, length(which(sensitivity_info[ddd, "drugid"] == drug))))
}

##check drugs screened
dd <- curationDrug[setdiff(curationDrug$unique.drugid, unique(sensitivity_info[, "drugid"])),"astrazeneca.numerical_drugid"]
for(drug in dd){
  for (screen in screens){
    AZ_data_1 <- read.csv(sprintf("../data/AZ/Raw_Data/%s", screen))
    ii <- which(as.character(AZ_data_1$DRUG_ID) == drug)
    print(length(ii))
  }
}
curationDrug[setdiff(curationDrug$unique.drugid, unique(sensitivity_info[, "drugid"])),
             c("astrazeneca.numerical_drugid")]


myfn <- "../data/AZ_recomputed.RData"
if(!file.exists(myfn)){
  recomputed <- PharmacoGx:::.calculateFromRaw(Az_raw_sensitivity, cap=100, nthread=4)
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


