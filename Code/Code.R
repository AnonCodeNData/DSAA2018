#setwd("~/Desktop/Academia/Manuscripts/2018 - SMOTEBoost (ECML)/")
#load("~/Desktop/Academia/Datasets/RegressionDatasets_FULL.RData")

setwd("~/Boosting")
load("RegressionDatasets_FULL.RData")

library(Metrics)
library(rpart)
library(UBL)
library(performanceEstimation)

source("Auxiliary.R")
source("Functions.R")

exp1 <- performanceEstimation(PredTask(DSs[[i]]@formula,DSs[[i]]@data),
  c(workflowVariants("wf.bagging",nbags=c(10,25,50)),
    workflowVariants("wf.rf",ntrees=c(100,250,500)),
    workflowVariants("AdaBoost.RT",t_final=c(50,100,200)),
    workflowVariants("AdaBoost.R2",t_final=c(50,100,200)),
    workflowVariants("AdaBoost.RTPlus",t_final=c(50,100,200)),
    workflowVariants("BEMBoost",t_final=c(50,100,200)),
    workflowVariants("wf.bagging_SM",perc.O=c(1.5,2,3),nbags=c(10,25,50)),
    workflowVariants("wf.rf_SM",perc.O=c(1.5,2,3),ntrees=c(100,250,500)),
    workflowVariants("SMOTEBoost.RT",perc.O=c(1.5,2,3),t_final=c(50,100,200)),
    workflowVariants("SMOTEBoost.R2",perc.O=c(1.5,2,3),t_final=c(50,100,200)),
    workflowVariants("SMOTEBoost.RTPlus",perc.O=c(1.5,2,3),t_final=c(50,100,200)),
    workflowVariants("SMOTEBoost.BEM",perc.O=c(1.5,2,3),t_final=c(50,100,200))),
  EstimationTask("totTime",method=CV(nReps = 2, nFolds=5))
)


res <- c()

for(wf in 1:72) {
  res_aux <- c()
  for(it in 1:10) {
    res_aux <- rbind(res_aux,getIterationsInfo(exp1,workflow=wf,task=1,it=it)$eval)
  }
  res_aux[is.na(res_aux)] <- 0
  res_aux <- as.data.frame(res_aux)
  res_aux["wf"] <- names(exp1[[1]])[wf]
  
  res <- rbind(res,res_aux)
}

write.csv(res,file=paste0("res",i,".csv"),row.names=FALSE)

agg <- aggregate(res,by=list(res$wf),FUN=mean)
head(agg[order(agg$F1.u,decreasing=TRUE),],8)
head(agg[order(agg$mu,decreasing=TRUE),],8)
