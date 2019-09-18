library(Hmisc) #capitalize
mybarplot <- function(Filename, data, barsNo, groupNo, group.labels, ylab.label, legend.lables, main.label, yaxis = c("Regular", "Log"), cex=1.1, hh, plot.highlight){
  #compare all
  barplot.colors= RColorBrewer::brewer.pal(n=7, name="Set1")[c(1:2,4,3,5:7)]#c("steelblue3","palegreen3", "mediumpurple3","darkorange3", "indianred3")
  Sp <- c(.75,.25,.25,.25,.25)
  yr <- c(0, max(data,na.rm=TRUE)*1.02)
  idx <- NULL
  k <- 0
  for(i in 1:round(ncol(data)/barsNo))
  {
    idx <- c(idx,k + round(barsNo/2))
    k <- k + barsNo
  }

  Groups.labels <- matrix(NA, nrow = barsNo, ncol = groupNo)
  Sp <- NULL
  legend.colors <- NULL
  for(i in 1:barsNo)
  {
    Groups.labels[i,] <- ""
    Sp <- c(Sp,.85)
  }
  Groups.labels[round(barsNo/2+1),] <- group.labels
  Sp[1] <- .15

  pdf(file = sprintf("%s.pdf", Filename) , height=7, width=14)
  par(mar=c(9,5,5,8), xpd=T)
  par(oma=c(2,2,2,2))
  temp <- NULL
  for(i in 1:barsNo) {
    temp <- rbind(temp, data[seq(i, length(data), by=barsNo)])
  }
  if(yaxis == "Log") {
    if(0 %in% temp){
      temp <- temp + 1
    }
    #mp <- barplot(log10(temp), beside=TRUE, col=barplot.colors[1:barsNo], ylim=log10(c(0, max(temp)) + 1), ylab=ylab.label, space=Sp, border=NA, axes=FALSE)
    mp <- barplot(log10(temp), beside=TRUE, col=barplot.colors[1:barsNo], ylim=log10(c(0, max(temp)) + 1), ylab=ylab.label, space=c(.15,.85), border=NA, axes=FALSE, main=main.label)
    text(mp, par("usr")[3], labels=Groups.labels, srt=50, adj=c(1.1,1.1), xpd=TRUE, cex=cex)
    magicaxis::magaxis(unlog = 'y', side=2, tcl=-.3)
    #legend("topright", inset=c(-.05,0), legend=legend.lables, fill=barplot.colors[1:barsNo], bty="n")
  }else {

    #    mp <- barplot(temp, beside=TRUE, col=barplot.colors[1:barsNo], ylim=c(0, max(temp)*1.02), ylab=ylab.label, space=Sp, border=NA, axes=FALSE)
    mp <- barplot(temp, beside=TRUE, col=barplot.colors[1:barsNo], ylim=c(0, 1), ylab=ylab.label, space=c(.15,.85), border=NA, axes=FALSE, main=main.label)
    abline(h=hh, lty=2)
    text(mp, par("usr")[3], labels=Groups.labels, srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=cex)
    magicaxis::magaxis(side=2, tcl=-.3)
    legend("topright", inset=c(-.05,0), legend=legend.lables, fill=barplot.colors[1:barsNo], bty="n")
    legend("bottomright", inset=c(-.1,0), legend=plot.highlight, bty="n")
  }

  dev.off()
}
myScatterPlot <- function(Name, x, y, method=c("plain", "transparent", "smooth"), legend.label, transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), ...){
  require(grDevices) || stop("Library grDevices is not available!")
  method <- match.arg(method)
  ccix <- complete.cases(x, y)
  x <- x[ccix]
  y <- y[ccix]

  if (sum(ccix) < minp) {
    ## too few points, no transparency, no smoothing
    if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch, ...) }
  } else {
    ## enough data points
    switch(method,
           "plain"={
             if(!missing(Name)){
               pdf(file = Name, height=7, width=7)
             }
             plot(x=x, y=y, col=col, pch=pch, ...)
             if(missing(legend.label)){
               #abline(0, 1, lty=2, col="gray")
               legend("topright", legend=sprintf("r=%.2f", cor(as.numeric(x),as.numeric(y), method="spearman")), cex=1.1, bty="n")
             }else{
               legend("topright", legend=legend.label, cex=1.1, bty="n")
             }
             if(!missing(Name)){
              dev.off()
             }
           },
           "transparent"={
             if(!missing(Name)){
               pdf(file = Name, height=7, width=7)
             }
             myrgb <- grDevices::col2rgb(col, alpha=FALSE) / 255
             plot(x=x, y=y, pch=pch, col=rgb(red=myrgb[1], green=myrgb[2], blue=myrgb[3], alpha=transparency, maxColorValue=1), ...)
             if(missing(legend.label)){
               legend("topright", legend=sprintf("r=%.2f", cor(as.numeric(x),as.numeric(y), method="spearman")), cex=1.1, bty="n")
             }else{
               legend("topright", legend=legend.label, cex=1.1, bty="n")
             }
             if(!missing(Name)){
               dev.off()
             }
           },
           "smooth"={
             if(!missing(Name)){
               pdf(file = Name, height=7, width=7)
             }
             smoothScatter(x=x, y=y, col="lightgray", colramp=colorRampPalette(smooth.col), pch=smooth.pch, ...)
             #abline(0,max(y)/max(x),col="red")
             abline(lm(y~x), col = "red")
             #lines(x=c(0,0), y=c(max(x),max(y)), lty="solid", col="red")
             if(missing(legend.label)){
               legend("topright", legend=sprintf("r=%.2f", cor(as.numeric(x),as.numeric(y), method="spearman")), cex=1.1, bty="n")
             }else{
               legend("topright", legend=legend.label, cex=1.1, bty="n")
             }
             if(!missing(Name)){
               dev.off()
             }
           }
    )
  }
}

featureInfoType <- function(dataset.name, type){
  switch(dataset.name,
         "CCLE"={
           switch(type,
                  "expression"=featureInfo(CCLE, "rna")
                  #"mutation"=featureInfo(CCLE, "mutation"),
                  #"cnv"=featureInfo(CCLE, "cnv"),
                  #"fusion"=featureInfo(GDSC, "fusion")
                  )
         },
         "gCSI"={
           switch(type,
                  "expression"=featureInfo(gCSI, "rnaseq"),
                  "mutation"=featureInfo(CCLE, "mutation"),
                  "cnv"=featureInfo(gCSI, "cnv")
                  #"fusion"=featureInfo(GDSC, "fusion")
                  )
         },
         "gCSI_GR_AOC_Marc"={
           switch(type,
                  "expression"=featureInfo(gCSI, "rnaseq"),
                  #"mutation"=featureInfo(CCLE, "mutation"),
                  "cnv"=featureInfo(gCSI, "cnv"))
                  #"fusion"=featureInfo(GDSC, "fusion"))
         },
         "gCSI_V1"={
           switch(type,
                  "expression"=featureInfo(gCSI, "rnaseq"),
                  "mutation"=featureInfo(gCSI, "mutation"),
                  "cnv"=featureInfo(gCSI, "cnv"))
                  #"fusion"=featureInfo(GDSC, "fusion"))
         },
         "CTRPv2"={
           switch(type,
                  "expression"=featureInfo(CCLE, "rnaseq"),
                  "mutation"=featureInfo(CCLE, "mutation"),
                  "cnv"=featureInfo(CCLE, "cnv"))
                  #"fusion"=featureInfo(GDSC, "fusion"))
         },
         "GDSC1000"={
           switch(type,
                  "expression"=featureInfo(GDSC1000, "rna"),
                  "mutation"=gdsc1000.mutation.fData,
                # "mutation"={
                #   xx <- featureInfo(GDSC, "mutation")
                #   colnames(xx)[which(colnames(xx) == "gene_name")] <- "Symbol"
                #   xx
                # },
                  "cnv"=featureInfo(GDSC, "cnv"),
                  "fusion"=featureInfo(GDSC, "fusion"))
         },
         "GDSC1000_updatedGRset"={
           switch(type,
                  "expression"=featureInfo(GDSC1000_updatedGRset, "rna"),
                  "mutation"=gdsc1000.mutation.fData,
                  # "mutation"={
                  #   xx <- featureInfo(GDSC, "mutation")
                  #   colnames(xx)[which(colnames(xx) == "gene_name")] <- "Symbol"
                  #   xx
                  # },
                  "cnv"=featureInfo(GDSC, "cnv"),
                  "fusion"=featureInfo(GDSC, "fusion"))
         })

}
molecular.profile <- function(dataset.name, type){
  switch(dataset.name,
         "CCLE"={
           switch(type,
                  "expression"=ccle.rna)
                  #"mutation"=ccle.mutation,
                  #"cnv"=ccle.cnv,
                  #"fusion"=gdsc.fusion)
         },
         "GDSC1000"={
           switch(type,
                  "expression"=gdsc.rna,
               #   "mutation"=gdsc.mutation,
                  "mutation"=gdsc1000.mutation,
                  "cnv"=gdsc.cnv,
                  "fusion"=gdsc.fusion)
         },
         "GDSC1000_updatedGRset"={
           switch(type,
                  "expression"=gdsc_GR.rna,
                  #   "mutation"=gdsc.mutation,
                  "mutation"=gdsc1000.mutation,
                  "cnv"=gdsc.cnv,
                  "fusion"=gdsc.fusion)
         },
         "gCSI"={
           switch(type,
                  "expression"=gcsi.rnaseq,
                  #"mutation"=gcsi.mutation,
                  "cnv"=gcsi.cnv)
                  #"fusion"=gdsc.fusion)
         },
         "gCSI_GR_AOC_Marc"={
           switch(type,
                  "expression"=gcsi.rnaseq,
                  #"mutation"=ccle.mutation,
                  "cnv"=gcsi.cnv)
                  #"fusion"=gdsc.fusion)
         },
         "gCSI_V1"={
           switch(type,
                  "expression"=gcsi.rnaseq,
                  "mutation"=gcsi.mutation,
                  "cnv"=gcsi.cnv)
                  #"fusion"=gdsc.fusion)
         },
         "CTRPv2"=
           switch(type,
                  "expression"=ccle.rnaseq,
                  "mutation"=ccle.mutation,
                  "cnv"=ccle.cnv))
                  #"fusion"=gdsc.fusion))
}
drug.sensitivity <- function(dataset.name){
  switch(dataset.name,
         "AZ"=AZ.drug.sensitivity,
         "GDSC1"=GDSCv1.drug.sensitivity,
         "GDSC2"=GDSCv2.drug.sensitivity,
         "CCLE"=ccle.drug.sensitivity,
         "gCSI"=gcsi.drug.sensitivity,
         "GRAY"=gray.drug.sensitivity,
         "gCSI_GR_AOC_Marc"=gcsi.drug.sensitivity,
         "gCSI_V1"=gcsi.drug.sensitivity,
         "CTRPv2"=ctrpv2.drug.sensitivity,
         "GDSC1000"=gdsc1000.drug.sensitivity,
         "GDSC1000_updatedGRset"=gdsc1000_GR.drug.sensitivity)
}
pset <- function(dataset.name){
  switch(dataset.name,
         "AZ"=AZ,
         "GDSC1"=GDSCv1,
         "GDSC2"=GDSCv2,
         "CCLE"=CCLE,
         "gCSI"=gCSI,
         "GRAY"=GRAY,
         "GRAY2013_updated"=GRAY2013_updated,
         "GRAY2017_updated"=GRAY2017_updated,
         "gCSI_GR_AOC_Marc"=gCSI,
         "gCSI_V1"=gCSI,
         "CTRPv2"=CTRPv2,
         "GDSC1000"=GDSC1000,
         "GDSC1000_updatedGRset"=GDSC1000_updatedGRset)
}
plot.multi.dens <- function(s, t, cols, xlim, xlab)
{
  if(missing(cols)){
    cols <- l:length(s)
  }
  junk.x = NULL
  junk.y = NULL
  for(i in 1:length(s)) {
    junk.x = c(junk.x, density(s[[i]])$x)
    junk.y = c(junk.y, density(s[[i]])$y)
  }
  if(missing(xlim)){
    xr <- range(junk.x)
  }else{
    xr <- xlim
  }
  yr <- range(junk.y)
  plot(density(s[[1]]), xlim = xr, ylim = yr, xlab=xlab, main = t)
  for(i in 1:length(s)) {
    lines(density(s[[i]]), xlim = xr, ylim = yr, col = cols[i])
  }
}
### Adapted from:
### https://gist.github.com/webbedfeet/7031404fc3f500f6258e
##Forestplot
credplot.gg <- function(d, xlab='Variable', ylab="Y", col="black"){
  # d is a data frame
  # d$x gives variable names
  # d$y gives center point
  # d$ylo gives lower limits
  # d$yhi gives upper limits
  # d$annot labels the graph
  require(ggplot2)
  p <- ggplot(d, aes(x=dataset, y=y, ymin=ymin, ymax=ymax, label=label, color=metric))+
    geom_pointrange(fatten=d$size, colour=col, position=position_dodge(width=0.5))+
    #geom_errorbar(data=d, aes(ymin=cindex, ymax=cindex), width=d$cindex_upper-d$cindex_lower) +
    geom_hline(yintercept = 0.5, linetype=2)+
    #geom_vline()
    coord_flip(ylim=c(0.45, max(d$ymax) + 0.1))+
    xlab(ylab) +
    ylab(xlab) + theme_bw() +geom_text(vjust = 0, nudge_x = 0.3, size=5)
  return(p)
}

plot_ci_heatmap <- function(file_name, data){
  if(nrow(data) < 2){
    print("At least two rows and columns required to be present in data for plotting heatmap.")
    return()
  }
  pdf(file=file_name, width=5, height=5)
  cex <- switch(as.character(nrow(data)),
                "69"=0.3,
                "131"=0.2,
                "27"=0.7,
                "26"=0.7,
                "51"=0.4,
                1)

  cc <- rep("black", nrow(data))
  xx <- which(!(1:nrow(data) == apply(data, 1, which.max)))
  if(length(xx) > 0){cc[xx] <- "red"}
  colRow <- cc

  cc <- rep("black", ncol(data))
  yy <- which(!(1:ncol(data) == apply(data, 2, which.max)))
  if(length(yy) > 0){cc[yy] <- "red"}
  colCol <- cc
  gplots::heatmap.2(data,
                    Colv=NA,
                    Rowv=NA,
                    dendrogram="none",
                    trace="none",
                    key=FALSE,
                    main="MCI",
                    col=colorRampPalette(c("green","red"))(10),
                    colRow=colRow,
                    colCol=colCol,
                    cexRow=cex,
                    cexCol=cex,
                    margins=c(8,8),
                    ylab=sprintf("%s, %s out of %s", ds1, length(xx), nrow(data)),
                    xlab=sprintf("%s, %s out of %s", ds2, length(yy), ncol(data)))

  #heatmap(mci.temp, Rowv=NA, Colv=NA,col = colorRampPalette(c("green","red"))(10), margins=c(8,8), cexRow=cex, cexCol=cex, colR)
  dev.off()
}
