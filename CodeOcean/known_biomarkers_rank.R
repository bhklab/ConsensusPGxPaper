types <- c("cnv", "mutation", "expression")
names(types) <- c("amplification",  "mutation",  "expression")

type <- "expression"
drug <- "AZD0530"
feature <- "IL4R"
load(sprintf("../results/genomewide_hits_%s.RData", drug), verbose=T)
genome_wide_hits <- genome_wide_hits[which(genome_wide_hits$type==type),]
genome_wide_hits[which(genome_wide_hits$gene==feature & genome_wide_hits$compound==drug),]
mCI_rank <- genome_wide_hits[order(genome_wide_hits$metamCI, decreasing = T), ]
cindex_rank <- genome_wide_hits[order(genome_wide_hits$metacindex, decreasing = T), ]
cindex <- which(cindex_rank$gene==feature)
mCI <- which(mCI_rank$gene==feature)

dataset.name <- "CCLE"
id <- rownames(featureInfoType(dataset.name, type))[match(feature, featureInfoType(dataset.name, type)[, "Symbol"])]
drug_sensitivity <- drug.sensitivity(dataset.name)
molecular_profile <- molecular.profile(dataset.name, type)
cell.lines <- intersect(colnames(molecular_profile), rownames(drug_sensitivity))

xx <- as.numeric(molecular_profile[id ,cell.lines])
yy <- drug_sensitivity[cell.lines, drug]
cor(xx, yy, use="pairwise.complete.obs")
myScatterPlot(x=xx, y=yy, xlab=type, ylab="drug sensitivity", main=sprintf("%s, %s", drug, feature))
paired.concordance.index(predictions=xx,
                         observations=yy,
                         delta.pred=0,
                         delta.obs=0.2,
                         alternative="greater",
                         logic.operator="and")[c("cindex", "p.value", "relevant.pairs.no", "lower", "upper", "sterr")]
paired.concordance.index(predictions=xx,
                         observations=yy,
                         delta.pred=0,
                         delta.obs=0,
                         alternative="greater",
                         logic.operator="and", outx=FALSE)[c("cindex", "p.value", "relevant.pairs.no", "lower", "upper", "sterr")]

