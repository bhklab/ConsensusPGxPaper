#' Generate drug response prediction models and use them to predict external test sets
#'
#'Function to generate drug response prediction models and predict external test sets
#'
#'@param LigandDesc A dataframe with ECFP4 descriptors for drugs and LNIC50 responses corresponding to different cell lines
#'@param MutCNV A dataframe with principal component scores of the mutation / CNV status of the cell lines
#'@return A dataframe with the various test set performances
#'@export


DR_Models<-function(LigandDesc_LNIC50,MutCNV)
{

#Write Test set performances to a file
Predictions_colnames<-c()
Predictions_colnames<-cbind(Predictions_colnames,"Test_RMSE","Test_Pearson_Correlation","Test_RMSE_cell","Test_Rsquared_cell","Test_RMSE_Drug","Test_Rsquared_Drug")
write.table(Predictions_colnames,file="Drugs_Testset_predictions_rf_70_30_splitting_Mut_CNV.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)


Drugs_list<-unique(as.vector(LigandDesc_LNIC50$Drug_Name))
registerDoMC(cores=20)

 set.seed(2)
 Predictions_all<-c()

 #Split the total number of drugs into 2 (70% Training; 30% Test)
 Training_drugs<-sample(Drugs_list,size=0.7*length(Drugs_list))
 Training_index<-which(LigandDesc_LNIC50$Drug_Name %in% Training_drugs)
 Training_matrix<-LigandDesc_LNIC50[Training_index,]
 Test_matrix<-LigandDesc_LNIC50[-Training_index,]


  #Combine Ligand LigandDesc_LNIC50riptors and Gene Expression
  Training_set<-merge(Training_matrix,MutCNV,by='CellLine')
  Test_set<-merge( Test_matrix,MutCNV,by='CellLine')

  Training_data_extract<-as.data.frame(Training_set[,!(colnames(Training_set)) %in% c("Drug_Name","CellLine","LN_IC50")])
  Test_data_extract<-as.data.frame(Test_set[,!(colnames(Test_set)) %in% c("Drug_Name","CellLine","LN_IC50")])
  Training_sens<-as.vector(Training_set$LN_IC50)
  Test_sens<-as.vector(Test_set$LN_IC50)


  Drugs_orig_Pred<-c()
  Predictions<-c()

  #Build random forest models on 70% of the data and use it to test the remaining 30%; Splits were done 5 times and run parallely in the cluster
  fullModel<-randomForest(as.matrix(Training_data_extract),Training_sens,importance=TRUE)
  #save(fullModel, file = "mymodel_1.rda")
  Predicted_values_test<-as.vector(predict(fullModel,Test_data_extract))
  Predictions<-as.data.frame(cbind(Drug=Test_matrix$Drug_Name,Cell=Test_matrix$CellLine,Test_sens,Predicted_values_test))

  #Overall R2 and RMSEs; All observations in the test set considered
  Test_Rsquared<-(cor(Test_sens,Predicted_values_test,method="pearson"))^2
  Test_RMSE<-round(rmse(Test_sens,Predicted_values_test),digits=3)
  CellLines_unique<-as.vector(unique(Predictions$Cell))
  Drugs_unique<-as.vector(unique(Predictions$Drug))
  Cell_Rsquared_all<-c()
  Drug_Rsquared_all<-c()

  #Test set Predictions for each cell line; R2 and RMSE are computed, considering the performance of each cell line
  for(j in 1:length(CellLines_unique))
  {
      cell_pred<-as.vector(CellLines_unique[j])
      Predictions_cell<-subset(Predictions,Predictions$Cell==cell_pred)
      Test_Rsquared_cell<-(cor(as.numeric(as.vector(Predictions_cell$Test_sens)),as.numeric(as.vector(Predictions_cell$Predicted_values_test)),method="pearson"))^2
      Test_RMSE_cell<-round(rmse(as.numeric(as.vector(Predictions_cell$Test_sens)),as.numeric(as.vector(Predictions_cell$Predicted_values_test))),digits=3)
      Cell_Rsquared<-as.data.frame(cbind(cell_pred,Test_Rsquared_cell,Test_RMSE_cell))
      Cell_Rsquared_all<-rbind(Cell_Rsquared_all,Cell_Rsquared)
  }

  #Test set Predictions for each drug; R2 and RMSE are computed, considering the performance of each drug
  for(k in 1:length(Drugs_unique))
  {
      Drug_pred<-as.vector(Drugs_unique[k])
      Predictions_drug<-subset(Predictions,Predictions$Drug==Drug_pred)
      Test_Rsquared_Drug<-(cor(as.numeric(as.vector(Predictions_drug$Test_sens)),as.numeric(as.vector(Predictions_drug$Predicted_values_test)),method="pearson"))^2
      Test_RMSE_Drug<-round(rmse(as.numeric(as.vector(Predictions_drug$Test_sens)),as.numeric(as.vector(Predictions_drug$Predicted_values_test))),digits=3)
      Drug_Rsquared<-as.data.frame(cbind(Drug_pred,Test_Rsquared_Drug,Test_RMSE_Drug))
      Drug_Rsquared_all<-rbind(Drug_Rsquared_all,Drug_Rsquared)
  }

  Predictions_all<-data.frame(RMSE_overall=Test_RMSE,Pearson_correlation_overall=Test_Rsquared,RMSE_cell=mean(as.numeric(as.vector(Cell_Rsquared_all$Test_RMSE_cell))),Pearson_correlation_cell=mean(as.numeric(as.vector(Cell_Rsquared_all$Test_Rsquared_cell))),RMSE_Drug=mean(as.numeric(as.vector(Drug_Rsquared_all$Test_RMSE_Drug))),Pearson_correlation_Drug=mean(as.numeric(as.vector(Drug_Rsquared_all$Test_Rsquared_Drug))))
  write.table(Predictions_all,file="Drugs_Testset_predictions_rf_70_30_splitting_Mut_CNV.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  return(Predictions_all)
}
