#How to create figures
1-Figure 1
The distribution of biological replicates AAC values is created with all_replicates.kci




# How to add other sensitivity metrics to PSets



Lets say you want to compute GR_AOC values for all the experiments in GDSC1000 (updated version of GDSC) you need to take the following steps:

1. Load the PSet (or download PSet if you want to do it outside of this environment)
2. Access to raw sensitivity data which is a three dimensional array
3. Compute GR_AOC values for all the experiments (Note that viability values ranges might exceed [0, 100] so you might need to truncate them to sit between [0,100])
4. Add the recomputed values to sensitivity$profiles as a new column.


```{r pseudo_code}
library(PharmacoGx)
load("data/GDSC1000.RData")
raw.sensitivity <- GDSC1000@sensitivity$raw
experiments <- vector(length=dim(raw.sensitivity)[1])
names(experiments) <- dimnames(raw.sensitivity)[[1]] 
GR_AOC <- sapply(names(experiments), function(exp, raw.sensitivity) {compute_GR_AOC(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"]))
GDSC1000@sensitivity$profiles [,"GR_AOC"] <- GR_AOC
```
After these steps you would be able to get a matrix of all GR_AOC values in GDSC10000,
in which compounds are in rows and cells are in columns using summarizeSensitivityProfiles function.
# How to use this capsule to validate known biomarkers

In run.sh script 
1. Assign any alternative metric (e.g. "GR_AOC") to sensitivity_metric 
2. Subset datasets to the datasets you have added this new metric to their profiles using the approach in previous section
3. Run the capsule
