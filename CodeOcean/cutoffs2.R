set.seed(42)
library(wCI)
if (!exists("data_path")){datapath <- "data"}
if (!exists("results_path")){datapath <- "results"}
if (!file.exists(data_path)){dir.create(data_path)}
if (!file.exists(results_path)){dir.create(results_path)}
replicates <- function(pset, pset_name, sens_metric, sens_factor=1, max_reps=10000, save_=TRUE, force_=FALSE, unique_replicates=FALSE){
  all.rep.pairs <- NULL
  reps_path <- sprintf("%s/%s_reps.RData", data_path, pset_name)
  if (!force_ & file.exists(reps_path)){
    load(reps_path)
    return(all.rep.pairs)
  }
  dd <- apply(pset@sensitivity$n, 2, function(x){length(which(x > 1))})
  dd.rep <- intersect(names(which(dd > 0)), drugNames(pset))
  all.rep.pairs <- NULL
  jj <- 1
  for(drug in dd.rep){
    kk <- 1
    cells <- names(which(pset@sensitivity$n[,drug] > 1))
    ii <- sensitivityInfo(pset)[which(sensitivityInfo(pset)$drugid==drug & sensitivityInfo(pset)$cellid %in% cells),]
    for(cell in cells){
      print(sprintf("%s (%s out of %s cells) out of %s", jj, kk, length(cells), length(dd.rep)))
      xx <- sensitivityProfiles(pset)[rownames(ii)[which(ii$cellid == cell)], sens_metric, drop=T]/sens_factor
      cell_specific <- NULL
      for(i in 1:(length(xx)-1)){
        if(pset_name %in% c("AZ", "GDSC_v1", "GDSC_v2")){
          cell_specific <- rbind(cell_specific, c(xx[i], xx[i+1]))
        }else{
          for(j in (i+1):length(xx)){
            cell_specific <- rbind(cell_specific, c(xx[i], xx[j]))
          }
        }
      }
      if(unique_replicates){
        rand_id <- sample(1:nrow(cell_specific), size=1)
        cell_specific <- cell_specific[rand_id, ]
      }
      all.rep.pairs <- rbind(all.rep.pairs, cell_specific)
      kk <- kk + 1
    }
    jj <- jj + 1
  }
  if(save_){
    save(all.rep.pairs, file=sprintf("%s/%s_all_reps.RData", data_path, pset_name))
  }
  if (nrow(all.rep.pairs) > max_reps){
    rand_id <- sample(1:nrow(all.rep.pairs), size=max_reps)
    all.rep.pairs <- all.rep.pairs[rand_id,]
  }
  if(save_){
    save(all.rep.pairs, file=reps_path)
  }
  return(all.rep.pairs)
}
ci_compute <- function(pset_name, predictions, observations, metrics, delta=0.2, save_=TRUE, force_=FALSE){
  ci_values <- list()
  ci_values_path <- sprintf("%s/%s_ci_values.RData", data_path, pset_name)
  if (!force_ & file.exists(ci_values_path)){
    load(ci_values_path)
    return(ci_values)
  }
  if ("mci" %in% metrics){
    mci <- wCI::paired.concordance.index(predictions, observations, delta.pred=delta, delta.obs=delta, logic.operator="or")
    ci_values[["mci"]] <- mci
  }
  if ("ci" %in% metrics){
    ci <- wCI::paired.concordance.index(predictions, observations, delta.pred=0, delta.obs=0, logic.operator="or")
    ci_values[["ci"]] <- ci
  }
  if ("kci" %in% metrics){
    kci <- wCI::paired.concordance.index.weighted.version(predictions, observations, delta.pred=0, delta.obs=0, logic.operator="or",
                                                          weightingFun_obs=wCI:::kernel_gaussian,
                                                          weightingFun_pred=wCI:::kernel_gaussian,
                                                          alternative="greater", CPP=TRUE)
    ci_values[["kci"]] <- kci
  }
  if ("kcil" %in% metrics){
    kcil <- wCI::paired.concordance.index.weighted.version(predictions, observations, delta.pred=0, delta.obs=0, logic.operator="or",
                                                           weightingFun_obs=wCI:::kernel_laplace,
                                                           weightingFun_pred=wCI:::kernel_laplace,
                                                           alternative="greater", CPP=TRUE)
    ci_values[["kcil"]] <- kcil
  }
  if(save_){
    save(ci_values, file=ci_values_path)
  }
  return(ci_values)
}
replicates_plot <- function(pset_name, reps, ci_list){
  x <- densCols(reps[,1], reps[,2], colramp=colorRampPalette(c("black", "white")))
  dens <- col2rgb(x)[1,] + 1L
  col <- cols[dens]
  pdf(sprintf("%s/%s_reps.pdf", results_path, pset_name), height=5, width=5)
  myScatterPlot(x=reps[,1], y=reps[,2],
                xlab="Replicate 1 (AAC)", ylab="Replicate 2 (AAC)",
                main=sprintf("All replicates in %s", pset_name),
                pch=20, method="plain", legend.label="", xlim=c(0,1), ylim=c(0,1), col=col)
  legend_label <- NULL
  for(i in 1:length(ci_list)){
    legend_label <- c(legend_label, sprintf("%s=%.2f", names(ci_list)[i], ci_list[i]))
  }
  legend_label <- paste(legend_label, collapse=", ")
  legend("topleft", legend=legend_label, bty="n")
  dev.off()
}
cutoff_estimate <- function(pset_name, reps, xlim=c(-.8, .8)){
  ## Distribution of delta AAC values across all replicates
  pdf(sprintf("%s/%s_cutoff.pdf", results_path, pset_name), height=5, width=5)
  par(mar=c(4.1, 4.1, 4.1, 5.1), xpd=TRUE)
  delta_rep_aac <- reps[,1] - reps[,2]
  hist(delta_rep_aac, breaks=100, col="gray",
       main=sprintf("All replicates in %s", pset_name), xlab="Delta AAC", xlim=xlim, axes=FALSE, tck=-.01)
  axis(1, at=c(-.6, -.4, -.2, 0, .2 , .4, .6), labels=c(-.6, -.4, -.2, 0, .2 , .4, .6), las=2)
  axis(2)
  CI_95 <- quantile(abs(delta_rep_aac), probs = seq(0, 1, 0.05))["95%"]
  abline(v=CI_95, col="red", lty=2,  xpd=FALSE)
  if(pset_name == "gCSI" | length(grep("GRAY", pset_name)) > 0){
    y = 300
  }else{
    y = 1000
  }
  text(x=CI_95+.1, y=y, labels=round(CI_95, digits=2), col="red", cex=.9)
  abline(v=-CI_95, col="red", lty=2,  xpd=FALSE)
  #abline(v=-CI_95, col="red", lty=2)
  legend("topright", inset=c(-.30,0), legend=sprintf("%s replicates", nrow(reps)), bty="n")
  dev.off()
  return(CI_95)
}
#'
#' #' ## Threshold distribution across drugs in AZ
#' ## ----AZ_drugs----------------------------------------------------------
#'
#' qq.AZ.avg <- apply(qq.AZ, 2, mean)
#' qq.AZ.max <- apply(qq.AZ, 2, max)
#' AZ.drug.max.delta.auc <- apply(qq.AZ, 2, which.max)
#' AZ.quantile.over.combined.drugs.delta.auc <- quantile(AZ.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]
#'
#' hist(qq.AZ, col="AZ", breaks=100, xlab="delta AAC", main="AZ")
#' legend("topright", legend=paste(c("90%","95%"), round(AZ.quantile.over.combined.drugs.delta.auc, digits=2), sep=":"), bty="n")
#'
