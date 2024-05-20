## code to prepare `DATASET_RV` dataset goes here
if (!requireNamespace("HARModel", quietly = TRUE)) {
  stop("Cannot simulate as package HARModel not installed.")
}
sample_RV<-matrix(NA,200,30)
for(i in 1:30){
  sim<-HARModel::HARSimulate(len=200, periods = c(1, 5, 22),
                   coef = c(0.01, 0.36 ,0.28 , 0.28), errorTermSD = 0.001)
  sample_RV[,i]<-sim@simulation
}
colnames(sample_RV)<-c(paste(rep("Var",30),1:30))

usethis::use_data(sample_RV, overwrite = TRUE)
