#!/bin/bash
if [ $# -eq 2 ]; then
    # assign the provided arguments to variables
    sensitivity_metric=$1
    datasets=$2
else
    # assign the default values to variables
    sensitivity_metric="GR_AOC"
    datasets="GDSC1000_updatedGRset"
fi
Rscript biomarkers_clarify_list.R 
Rscript biomarkers_assess.R $sensitivity_metric $datasets

#Rscript biomarkers_plot.R $datasets
