datasets <- c("GRAY")
sensitivity_metric <- c("auc_recomputed", "AAC", "auc_recomputed", "auc_recomputed", "auc_recomputed")
sensitivity_factor <- c(1, 100, 1, 1, 100)
names(sensitivity_metric) <- names(sensitivity_factor) <- datasets
source("pset_loading.R")
source("foo.R")
data_path <- "../data"
results_path <- "../results_Cello"
source("cutoffs2.R")
ci_values <- NULL
all1 <- all2 <- d_name <- NULL
metrics <- c("ci", "kci", "mci")
names(metrics) <- c("CI", "kCI", "rCI")
cols <-  colorRampPalette(c("#2b8cbe", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100"))(256)

datasets <- c("GRAY", "GRAY2013_updated", "GRAY2017_updated")
sensitivity_metric <- c("auc_recomputed","auc_recomputed", "auc_recomputed")
sensitivity_factor <- c(1, 100, 100)
names(sensitivity_metric) <- names(sensitivity_factor) <- datasets
for(dataset in datasets){
  ps <- pset(dataset)
  pdf(sprintf("%s/%s_hist_aac.pdf", results_path, dataset) , height=5, width=5)
  hist(ps@sensitivity$profiles$auc_recomputed/sensitivity_factor[dataset], main=dataset, xlab="AAC")
  dev.off()
   reps <- replicates(pset=ps,
                      pset_name=dataset,
                     sens_metric=sensitivity_metric[dataset],
                     sens_factor=sensitivity_factor[dataset])
  ci_values_temp <- ci_compute(dataset,
                               predictions=reps[,1],
                               observations=reps[,2],
                               metrics=metrics)
  xx <- dataset
  ci_list <- list()
  for(i in 1:length(metrics)){
    xx <- c(xx, ci_values_temp[[metrics[i]]]$cindex)
    ci_list[[names(metrics)[i]]] <- ci_values_temp[[metrics[i]]]$cindex
  }
  ci_values <- rbind(ci_values, xx)
  replicates_plot(pset_name=dataset, reps=reps, ci_list=ci_list)
  delta <- cutoff_estimate(pset_name=dataset, reps=reps)

  all1 <- c(all1, reps[,1])
  all2 <- c(all2, reps[,2])
  d_name <- c(d_name, rep(pSetName(ps), nrow(reps)))
}
all_replicates <- cbind("rep1"=all1, "rep2"=all2, "dataset"=d_name)
ci_values_temp <- ci_compute("ALL",
                             predictions=all1,
                             observations=all2,
                             metrics=metrics[c("CI", "rCI")])
xx <- "ALL"
ci_list <- list()
for(i in 1:length(metrics)){

  if(names(metrics)[i] == "kCI"){
    ci_list[[names(metrics)[i]]] <- 0.89
    xx <- c(xx, 0.89)
  }else{
    ci_list[[names(metrics)[i]]] <- ci_values_temp[[metrics[i]]]$cindex
    xx <- c(xx, ci_values_temp[[metrics[i]]]$cindex)
  }
}
ci_values <- rbind(ci_values, xx)
replicates_plot(pset_name="ALL", reps=cbind(all1, all2), ci_list=ci_list)
delta <- cutoff_estimate(pset_name="ALL", reps=cbind(all1, all2))



