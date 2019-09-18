library("compareDrugScreens")
data("gcsi.genomics.mask")
data("gcsi.genomics")
gcsi.genomics[which(gcsi.genomics.mask)] <- NA
table(sapply(colnames(gcsi.genomics), function(x){strsplit(x, ":")[[1]][1]}))
#loh <- grep("loh.GeneID", colnames(gcsi.genomics))
#quantile(sapply(1:length(loh), function(x){length(table(gcsi.genomics[,loh[x]]))}))

hot <- grep("hot.GeneID", colnames(gcsi.genomics))
quantile(sapply(1:length(hot), function(x){length(table(gcsi.genomics[,hot[x]]))}))

mut <- grep("mut.GeneID", colnames(gcsi.genomics))
quantile(sapply(1:length(mut), function(x){length(table(gcsi.genomics[,mut[x]]))}))

data("gcsi.line.info")


gsci.mut <- gcsi.genomics[,c(hot, mut)]
rownames(gsci.mut) <- gcsi.line.info[match(rownames(gcsi.genomics), rownames(gcsi.line.info)), "CellLineName"]
gsci.mut.bin <- apply(gsci.mut, 2, function(x){xx <- table(x); ifelse(x==names(xx)[order(as.numeric(names(xx)))[1]], 0, 1)})
#library(org.Hs.eg.db)
keys <- sapply(sapply(colnames(gsci.mut.bin), function(x){strsplit(x, ":")[[1]][2]}), function(x){unlist(strsplit(x, "_")[[1]])})
#mapping <- AnnotationDbi::select(org.Hs.eg.db,
#                                keytype="ENTREZID",
#                                keys=,
#                                columns=c("SYMBOL"))
data("gcsi.genomics.feature.info")
gene_ids <- sapply(colnames(gsci.mut.bin), function(x){strsplit(x, "_")[[1]][1]})
gene_ids <- gsub("mut.", "", gene_ids)
gene_ids <- gsub("hot.", "", gene_ids)
gene_symbols <- gcsi.genomics.feature.info[match(gene_ids, rownames(gcsi.genomics.feature.info)), "Symbol"]
gene_ids <- colnames(gsci.mut.bin)
names(gene_symbols) <- colnames(gsci.mut.bin)
gcsi_mutations <- matrix(NA, nrow=nrow(gsci.mut.bin), ncol=length(unique(gene_symbols)), dimnames=list(rownames(gsci.mut.bin), unique(gene_symbols)))
for(gene in  unique(gene_symbols)){
  ii <- which(gene_symbols == gene)
  if(length(ii) > 1){
    gcsi_mutations[,gene] <- apply(gsci.mut.bin[, names(gene_symbols)[ii]], 1, function(x){ifelse(sum(x, na.rm=T)>0,1,0)})
  }else{
    gcsi_mutations[,gene] <- gsci.mut.bin[, names(gene_symbols)[ii]]
  }
}
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
cell_annotation_all <- read.csv("~/Documents/GitHub/PharmacoGx-private/inst/extdata/cell_annotation_all.csv", na.string=c("", " "))

cellid <- cell_annotation_all$unique.cellid#rownames(gCSI@cell)
cells <- unique(rownames(gcsi_mutations))

cells[which(cells=="786-O")] <- "786-0"
cells_curated <- tolower(gsub(badchars, "", cells))
cellid_curated <- tolower(gsub(badchars, "", cellid))
min(length(cells_curated), length(cellid_curated)) - length(intersect(cells_curated, cellid_curated))
gcsi_cellid <- cellid[match(cells_curated, cellid_curated)]

load("~/Google Drive/Psets/gCSI_V1.RData", verbose=T)
options(stringsAsFactors = FALSE)

gcsi_mutations_pdata <- as.data.frame(cbind("CellLineName"=rownames(gcsi_mutations), "cellid"=NA, "batchid"=NA))

gcsi_mutations_pdata$CellLineName[which(gcsi_mutations_pdata$CellLineName == "SW 480")] <- "SW480"
gcsi_mutations_pdata$CellLineName[which(gcsi_mutations_pdata$CellLineName == "786-O")] <- "786-0"
gcsi_mutations_pdata$CellLineName[which(gcsi_mutations_pdata$CellLineName == "HCT 116")] <- "HCT-116"
gcsi_mutations_pdata$CellLineName[which(gcsi_mutations_pdata$CellLineName == "PLC/PRF/5")] <- "PLC-PRF-5"
gcsi_mutations_pdata$CellLineName[which(gcsi_mutations_pdata$CellLineName == "23132/87")] <- "23132-87"
gcsi_mutations_pdata$CellLineName[which(gcsi_mutations_pdata$CellLineName == "SW 403")] <- "SW403"
gcsi_mutations_pdata$cellid <- gCSI@curation$cell$unique.cellid[match(gcsi_mutations_pdata$CellLineName, gCSI@curation$cell$gCSI.cellid)]
rownames(gcsi_mutations) <- rownames(gcsi_mutations_pdata) <- gcsi_mutations_pdata$cellid

gcsi_mutations_fdata <- as.data.frame(cbind("Symbol"=colnames(gcsi_mutations), "BEST"=TRUE))
rownames(gcsi_mutations_fdata) <- gcsi_mutations_fdata$Symbol
gcsi_mutations_eset <- Biobase::ExpressionSet(t(gcsi_mutations))
pData(gcsi_mutations_eset) <- gcsi_mutations_pdata
fData(gcsi_mutations_eset) <- gcsi_mutations_fdata
annotation(gcsi_mutations_eset) <- "mutation"

gCSI_V1 <- PharmacoSet(molecularProfiles=list("rnaseq"=gCSI@molecularProfiles[["rnaseq"]],
                                              "cnv"=gCSI@molecularProfiles[["cnv"]],
                                              "mutation"=gcsi_mutations_eset),
                       name="gCSI",
                       cell=gCSI@cell,
                       drug=gCSI@drug,
                       sensitivityInfo=gCSI@sensitivity$info,
                       sensitivityRaw=NULL,
                       sensitivityProfiles=gCSI@sensitivity$profiles,
                       sensitivityN=NULL,
                       curationCell=gCSI@curation$cell,
                       curationDrug=gCSI@curation$drug,
                       curationTissue=gCSI@curation$tissue,
                       datasetType="sensitivity")
warnings()
gCSI <- gCSI_V1
#gCSI@molecularProfiles <- gCSI@molecularProfiles[c("rnaseq", "cnv")]
save(gCSI, file="~/Google Drive/Psets/gCSI_V1.RData")

