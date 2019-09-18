#system(paste("Rscript biomarkers_assess_all_types.R", "GR_AOC", "gCSI_V1", "long_list"))
#system(paste("Rscript biomarkers_assess_all_types.R", "auc_recomputed", "gCSI_V1", "long_list"))
#system(paste("Rscript biomarkers_assess_all_types.R", "GReff", "gCSI_V1", "long_list"))
#system(paste("Rscript biomarkers_assess_all_types.R", "auc_recomputed", "CCLE,gCSI_V1,CTRPv2,GDSC1000", "long_list"))

#system(paste("Rscript biomarkers_assess_all_types.R", "auc_recomputed", "gCSI,CTRPv2,GDSC1000", "long_list"))
#system(paste("Rscript biomarkers_assess_all_types.R", "GR_AOC", "GDSC1000_updatedGRset", "long_list"))

options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
#args <- c("auc_recomputed", "CCLE,gCSI_V1,CTRPv2,GDSC1000", "long_list")
#args <- c("auc_recomputed", "CCLE,gCSI_V1,CTRPv2,GDSC1000", "ben_list")
sensitivity_metric <- as.character(args[1])
datasets <- unlist(strsplit(as.character(args[2]), ","))
biomarkers_list <- as.character(args[3])
alpha <- 0.05
delta=0.2
if(length(datasets)==1 && datasets == "gCSI_V1"){
  delta=0.17
  if(sensitivity_metric=="GR_AOC"){
    delta=0.3
  }
  if(sensitivity_metric=="GReff"){
    delta=0.47
  }
}
alternative <- "two.sided"
pairs_cutoff = 3000
path.results <- sprintf("../results_pairs_%s", pairs_cutoff)
#path.results <- sprintf("../results3")
if(length(datasets)==1 && datasets == "gCSI_V1"){
  path.results <- sprintf("%s_gcsi_%s", path.results, tolower(sensitivity_metric))
}
if(!dir.exists(path.results)){
  dir.create(path.results)
}
source("pset_loading.R")
#devtools::install_github("bhklab/mci")
library(mCI)
source("foo.R")
if(length(datasets) == 4){
  cc <- combn(4, 3)
  drugs <- NULL
  for(j in 1:ncol(cc)){
    dd <- list()
    for(i in 1:nrow(cc)){
      dd[[i]] <- colnames(drug.sensitivity(datasets[cc[i,j]]))
    }
    dd <-intersectList(dd[[1]], dd[[2]], dd[[3]])
    drugs <- union(drugs, dd)
  }
}else{
  ###Read the list of known biomarkers extracted from
  ###This list hast been already curated to only keep those markers with clear definition of the type of association between genotype and phenotype
  ##Assess the association between markers and drug response across all possible molecular data types and create seperate list for each type
  if(datasets == "gCSI_V1"){
    cells <- unique(gCSI@sensitivity$info$cellid[which(!is.na(gCSI@sensitivity$profiles$GR_AOC))])
    drugs <- unique(gCSI@sensitivity$info$drugid[which(!is.na(gCSI@sensitivity$profiles$GR_AOC))])
    gCSI <- subsetTo(pSet=gCSI, cells=cells, drugs=drugs)
  }
}

types <- c("cnv", "mutation", "expression")
##integrate the list of significant associations across molecular data types
for(ttype in types){
  G2P_clarified_hits <- NULL
  if(biomarkers_list == "long_list"){
    G2P_clarified <- read.csv(file="../data/G2P_clarified.csv")
    G2P_clarified <- G2P_clarified[which(G2P_clarified$compound %in% drugs),]
    ids <- rownames(featureInfoType("gCSI_V1", ttype))[match(G2P_clarified$gene, featureInfoType("gCSI_V1", ttype)[, "Symbol"])]
    if (any(is.na(ids))){
      G2P_clarified <- G2P_clarified[-which(is.na(ids)),]
    }
  }else if (biomarkers_list == "ben_list"){
    #G2P_clarified <- read.csv(file="../data/G2P_clarified.csv")
    G2P_clarified <- read.csv(file="../data/Ben_short_list.csv")
    #G2P_clarified <- G2P_clarified[which(G2P_clarified$compound %in% drugs),]
    ids <- rownames(featureInfoType("gCSI_V1", ttype))[match(G2P_clarified$gene, featureInfoType("gCSI_V1", ttype)[, "Symbol"])]
    if (any(is.na(ids))){
      G2P_clarified <- G2P_clarified[-which(is.na(ids)),]
    }

  }else{
    G2P_clarified <- read.csv(file="../data/known.biomarkers.csv")
    colnames(G2P_clarified)[which(colnames(G2P_clarified) == "drug")] <- "compound"
    tt <- which(colnames(G2P_clarified) == "type")
    G2P_clarified <- G2P_clarified[,-tt]
  }
  xx <- sprintf("%s_%s", G2P_clarified$compound, G2P_clarified$gene)
  ###The original list is in variant level but here for thses analyses we collapse the list to gene level
  unique_G2P_clarified <- match(unique(xx), xx)
  G2P_clarified <- G2P_clarified[unique_G2P_clarified,]
  G2P_clarified[,c(sprintf("%s_mCI", datasets), sprintf("%s_mCI_pvalue", datasets), sprintf("%s_mCI_pairs", datasets), sprintf("%s_mCI_lower", datasets), sprintf("%s_mCI_upper", datasets), sprintf("%s_mCI_se", datasets))] <- NA
  G2P_clarified[,c(sprintf("%s_cindex", datasets), sprintf("%s_cindex_pvalue", datasets), sprintf("%s_cindex_pairs", datasets), sprintf("%s_cindex_lower", datasets), sprintf("%s_cindex_upper", datasets), sprintf("%s_cindex_se", datasets))] <- NA
  ### The valididty of associations is assessed by both CI and mCI using mCI package functions
  for(i in 1:nrow(G2P_clarified)){
    drug <- G2P_clarified$compound[i]
    feature <- G2P_clarified$gene[i]
    for(dataset.name in datasets){
      id <- rownames(featureInfoType(dataset.name, ttype))[match(feature, featureInfoType(dataset.name, ttype)[, "Symbol"])]
      if(!is.null(id) && (drug %in% colnames(drug.sensitivity(dataset.name)) & !is.na(id))){
        molecular_profile <- molecular.profile(dataset.name, ttype)
        drug_sensitivity <- drug.sensitivity(dataset.name)
        cell.lines <- intersect(colnames(molecular_profile), rownames(drug_sensitivity))
        mm <- as.numeric(molecular_profile[id, which(!is.na(drug_sensitivity[cell.lines, drug]))])
        flag <- FALSE;
        if(ttype=="mutation"){
          flag <- length(which(mm == 1)) > 5 &&  (length(which(mm == 1))/length(which(!is.na(mm)))) > 0.01
        }else{
          flag <- range(mm, na.rm=T)[2]-range(mm, na.rm=T)[1] > 1
        }

        if(flag){

            G2P_clarified[i, grep(sprintf("%s_cindex", dataset.name),
                                  colnames(G2P_clarified))] <- paired.concordance.index(predictions=as.numeric(molecular_profile[id ,cell.lines]),
                                                                                        observations=drug_sensitivity[cell.lines, drug],
                                                                                        delta.pred=0,
                                                                                        delta.obs=0,
                                                                                        alternative=alternative,
                                                                                        logic.operator="and")[c("cindex", "p.value", "relevant.pairs.no", "lower", "upper", "sterr")]
            if(G2P_clarified[i, sprintf("%s_cindex_pairs", dataset.name)] < pairs_cutoff){
              G2P_clarified[i, sprintf("%s_cindex", dataset.name)] = 0.50
              G2P_clarified[i, sprintf("%s_cindex_pvalue", dataset.name)] = 1
              G2P_clarified[i, sprintf("%s_cindex_lower", dataset.name)] = 0.50-.01
              G2P_clarified[i, sprintf("%s_cindex_upper", dataset.name)] = 0.50

            }
          # G2P_clarified[i, grep(sprintf("%s_mCI", dataset.name),
          #                       colnames(G2P_clarified))] <- paired.concordance.index.weighted.version(predictions=as.numeric(molecular_profile[id ,cell.lines]),
          #                                                                             observations=as.numeric(drug_sensitivity[cell.lines, drug]),
          #                                                                             delta.pred=0,
          #                                                                             delta.obs=0,
          #                                                                             weightingFun_obs=kernel_gaussian,#weighting_gaussian,
          #                                                                             weightingFun_pred=kernel_gaussian,#weighting_gaussian,
          #                                                                             alternative="greater",
          #                                                                             logic.operator="and")[c("cindex", "p.value", "relevant.pairs.no", "lower", "upper", "sterr")]


            G2P_clarified[i, grep(sprintf("%s_mCI", dataset.name),
                                colnames(G2P_clarified))] <- paired.concordance.index(predictions=as.numeric(molecular_profile[id ,cell.lines]),
                                                                                      observations=drug_sensitivity[cell.lines, drug],
                                                                                      delta.pred=0,
                                                                                      delta.obs=delta,
                                                                                      alternative=alternative,
                                                                                      logic.operator="and")[c("cindex", "p.value", "relevant.pairs.no", "lower", "upper", "sterr")]
            if(G2P_clarified[i, sprintf("%s_mCI_pairs", dataset.name)] < pairs_cutoff){
              G2P_clarified[i, sprintf("%s_mCI", dataset.name)] = 0.50
              G2P_clarified[i, sprintf("%s_mCI_pvalue", dataset.name)] = 1
              G2P_clarified[i, sprintf("%s_mCI_lower", dataset.name)] = 0.50-.01
              G2P_clarified[i, sprintf("%s_mCI_upper", dataset.name)] = 0.50
            }
        }else{
          G2P_clarified[i, sprintf("%s_cindex", dataset.name)] = 0.50
          G2P_clarified[i, sprintf("%s_cindex_lower", dataset.name)] = 0.50-.01
          G2P_clarified[i, sprintf("%s_cindex_upper", dataset.name)] = 0.50
          G2P_clarified[i, sprintf("%s_cindex_pvalue", dataset.name)] = 1
          G2P_clarified[i, sprintf("%s_mCI", dataset.name)] = 0.50
          G2P_clarified[i, sprintf("%s_mCI_lower", dataset.name)] = 0.50-.01
          G2P_clarified[i, sprintf("%s_mCI_upper", dataset.name)] = 0.50
          G2P_clarified[i, sprintf("%s_mCI_pvalue", dataset.name)] = 1
        }
      }
    }
  }
  G2P_clarified <- cbind(G2P_clarified, cbind("metamCI"=NA,
                                              "metamCI_lower"=NA,
                                              "metamCI_upper"=NA,
                                              "metamCI_pvalue"=NA,
                                              "metacindex"=NA,
                                              "metacindex_lower"=NA,
                                              "metacindex_upper"=NA,
                                              "metacindex_pvalue"=NA))
  meta_index <- grep("meta", colnames(G2P_clarified))
  if(length(datasets) == 1){
    G2P_clarified[, meta_index] <- c(G2P_clarified[,grep("_mCI$", colnames(G2P_clarified))],
                                     G2P_clarified[,grep("_mCI_lower$", colnames(G2P_clarified))],
                                     G2P_clarified[,grep("_mCI_upper$", colnames(G2P_clarified))],
                                     G2P_clarified[,grep("_mCI_pvalue$", colnames(G2P_clarified))],
                                     G2P_clarified[,grep("_cindex$", colnames(G2P_clarified))],
                                     G2P_clarified[,grep("_cindex_lower$", colnames(G2P_clarified))],
                                     G2P_clarified[,grep("_cindex_upper$", colnames(G2P_clarified))],
                                     G2P_clarified[,grep("_cindex_pvalue$", colnames(G2P_clarified))])

  }else{
    for(i in which(apply(G2P_clarified, 1, function(x){length(which(is.na(x[grep("_mCI$", names(x))])))}) != length(datasets))){
      mci <- unlist(G2P_clarified[i, grep("_mCI$", colnames(G2P_clarified))]); ii <- which(!is.na(mci)); mci <- mci[ii]
      mci_se <- unlist(G2P_clarified[i, grep("_mCI_se$", colnames(G2P_clarified))]);mci_se <- mci_se[ii]
      mci_p <- unlist(G2P_clarified[i, grep("_mCI_pvalue$", colnames(G2P_clarified))]);mci_p <- mci_p[ii]
      mci_pair <- unlist(G2P_clarified[i, grep("_mCI_pairs$", colnames(G2P_clarified))]);mci_pair <- mci_pair[ii]
      meta_mci <- survcomp::combine.est(x=mci, x.se=mci_se, na.rm=TRUE, hetero=TRUE)

      #meta_mci_pvalue <- survcomp::combine.test(p=mci_p, hetero=FALSE,method="z.transform", weight=mci_pair)
      meta_mci_ci <- qnorm(p=alpha/2, lower.tail=FALSE) * meta_mci$se
      p <- pnorm((meta_mci$estimate - 0.5) / meta_mci$se)
      meta_mci_pvalue <- switch(alternative, less=p, greater=1 - p, two.sided=2 * min(p, 1 - p))

      cindex <- unlist(G2P_clarified[i, grep("_cindex$", colnames(G2P_clarified))]); ii <- which(!is.na(cindex)); cindex <- cindex[ii]
      cindex_se <- unlist(G2P_clarified[i, grep("_cindex_se$", colnames(G2P_clarified))]);cindex_se <- cindex_se[ii]
      cindex_p <- unlist(G2P_clarified[i, grep("_cindex_pvalue$", colnames(G2P_clarified))]);cindex_p <- cindex_p[ii]
      cindex_pair <- unlist(G2P_clarified[i, grep("_cindex_pairs$", colnames(G2P_clarified))]);cindex_pair <- cindex_pair[ii]
      meta_cindex <- survcomp::combine.est(x=cindex, x.se=cindex_se, na.rm=TRUE, hetero=TRUE)

      #meta_cindex_pvalue <- survcomp::combine.test(p=cindex_p, hetero=FALSE,method="z.transform", weight=cindex_pair)
      #meta_cindex_ci <- qnorm(p=alpha/2, lower.tail=FALSE) * meta_cindex$se
      meta_cindex_ci <- qnorm(p=alpha/2, lower.tail=FALSE) * meta_cindex$se
      p <- pnorm((meta_cindex$estimate - 0.5) / meta_cindex$se)
      meta_cindex_pvalue <- switch(alternative, less=p, greater=1 - p, two.sided=2 * min(p, 1 - p))

      G2P_clarified[i, meta_index] <- c(meta_mci$estimate,
                                        meta_mci$estimate - meta_mci_ci,
                                        meta_mci$estimate + meta_mci_ci,
                                        meta_mci_pvalue,
                                        meta_cindex$estimate,
                                        meta_cindex$estimate - meta_cindex_ci,
                                        meta_cindex$estimate + meta_cindex_ci,
                                        meta_cindex_pvalue)
    }
  }

  G2P_clarified_hits <- G2P_clarified
  G2P_clarified_hits$metamCI_fdr <- G2P_clarified_hits$metamCI_pvalue
  G2P_clarified_hits$metacindex_fdr <- G2P_clarified_hits$metacindex_pvalue

  G2P_clarified_hits$metamCI_fdr <- p.adjust(G2P_clarified_hits$metamCI_fdr, method="bonferroni")
  G2P_clarified_hits$metacindex_fdr <- p.adjust(G2P_clarified_hits$metacindex_fdr, method="bonferroni")

  all_associations <- G2P_clarified_hits
  colnames(all_associations) <- gsub("mCI", "rCI", colnames(all_associations))
  colnames(all_associations) <- gsub("fdr", "qvalue", colnames(all_associations))
  all_associations <- all_associations[order(all_associations$metarCI_qvalue),]
  all_associations <- all_associations[,-1]
  all_associations <- cbind("type"=ttype, all_associations)
  #all_associations <- all_associations[,-which(colnames(all_associations)=="ttype")]
  write.csv(all_associations, file=file.path(path.results, sprintf("all_associations_%s.csv", ttype)))
  if(biomarkers_list != "ben_list"){
    G2P_clarified_hits <- G2P_clarified_hits[which((G2P_clarified_hits$metamCI_fdr < 0.05 &
                                                      abs(G2P_clarified_hits$metamCI - .5) > 0.05)|
                                                     (G2P_clarified_hits$metacindex_fdr < 0.05 &
                                                        abs(G2P_clarified_hits$metacindex - .5) > 0.55)),]

  }
  G2P_clarified_hits <- G2P_clarified_hits[order(G2P_clarified_hits$compound, G2P_clarified_hits$gene),]
  nrow(G2P_clarified_hits)
  save(G2P_clarified_hits, file=file.path(path.results, sprintf("G2P_clarified_hits_all_%s.RData", ttype)))
}

## Create the plots and reports of the validated biomarkers
types <- c("cnv", "mutation", "expression")
for(ttype in types){
  load(file.path(path.results, sprintf("G2P_clarified_hits_all_%s.RData", ttype)))
  drugs_types <- list()
  dd <- unique(G2P_clarified_hits$compound)
  mycol <- c(c(colorRamps::primary.colors()[c(1:6, 10:15)], RColorBrewer::brewer.pal(n=8, name="Set2")[2:8]));
  if(length(mycol) < length(dd)){
    mycol <- rep(mycol, ceiling( length(dd)/length(mycol)))
    length(mycol)
  }
  mycol <- mycol[1:length(dd)]
  names(mycol) <- dd
  temp <- G2P_clarified_hits
  if(nrow(G2P_clarified_hits) > 0){
    G2P_clarified_hits <-
      temp[sapply(unique(temp$compound), function(x){
        min(which(temp$compound==x & temp$metacindex==max(temp[which(temp$compound == x), "metacindex"])))
      }),]
    print(ttype)
    print(length(unique(G2P_clarified_hits$compound)))
    print(nrow(G2P_clarified_hits))
    gain=which((abs(G2P_clarified_hits$metacindex-.5) <= 0.05 | G2P_clarified_hits$metacindex_fdr >= 0.05) &
                 (abs(G2P_clarified_hits$metamCI-.5) > 0.05 & G2P_clarified_hits$metamCI_fdr < 0.05))
    print(length(gain))
    print(length(gain)/nrow(G2P_clarified_hits))
    loss=which((abs(G2P_clarified_hits$metamCI-.5) <= 0.05| G2P_clarified_hits$metamCI_fdr  >= 0.05) &
                 (abs(G2P_clarified_hits$metacindex-.5) > 0.05 & G2P_clarified_hits$metacindex_fdr < 0.05))
    print(length(loss))
    print(length(loss)/nrow(G2P_clarified_hits))
    #G2P_clarified_hits[loss,]
    drugs_types[[ttype]] <- G2P_clarified_hits$compound
    save(G2P_clarified_hits, file=sprintf("%s/G2P_clarified_hits_%s.RData", path.results, ttype))
    save(G2P_clarified_hits, file=file.path(path.results, "G2P_clarified_hits.RData"))
    source("biomarkers_plot.R")
  }
}
    #source("Check_Biomarkers_cor.R")

types <- c("cnv", "mutation", "expression")
all_associations <- NULL
for(ttype in types){
  temp <- read.csv(file.path(path.results, sprintf("all_associations_%s.csv", ttype)))
  all_associations <- rbind(all_associations, temp)
}
all_associations <- all_associations[,-1]
all_associations <- all_associations[order(as.numeric(all_associations$metarCI_pvalue)), ]
write.csv(all_associations, file.path(path.results, "all_associations.csv"))
#intersectList(drugs_types)
#RAFA
for(ttype in types){
  temp <- read.csv(file.path(path.results, sprintf("all_associations_%s.csv", ttype)))
  temp <- temp[,c(1:5, grep('_rCI_pvalue', colnames(temp)), grep('metarCI_qvalue', colnames(temp)))]
  for(x in 6:10){
    xx <- which(temp[,x]<=0.05)
    yy <- which(temp[,x]>0.05)
    temp[xx,x] <- "significant"
    temp[yy,x] <- "not"
  }
  write.csv(temp, file.path(path.results, sprintf("summarized_associations_%s.csv", ttype)))

}


