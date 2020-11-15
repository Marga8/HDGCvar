## code to prepare `DATASET` dataset goes here
SimulVAR=function(T_,g){
  coef1<-matrix(NA,nrow=g,ncol=g)
  for (i in 1:g) {
    for (j in 1:g) {
      coef1[i,j]<-((-1)^(abs(i-j)))*(0.4^(abs(i-j)+1))
    }
  }
  presample<-1
  T_new<-T_+presample
  eps1<-rnorm(ncol(coef1)*T_new,0,1)
  eps<-matrix(eps1,nrow=ncol(coef1))
  X <- matrix(nrow=ncol(coef1),ncol=T_new)
  X[,1] <- eps[,1]
  for (t in 2:T_new) {
    X[,t] <- (coef1)%*%X[,t-1]+eps[,t]
  }
  finseries<- X[,(1+presample):T_new]
  return(t(finseries))
}

set.seed(123)
dataset_I0<-as.matrix(SimulVAR(200,30))
sample_dataset_I1<-as.matrix(diffinv(dataset_I0))
colnames(sample_dataset_I1)<-c(paste(rep("Var",30),1:30))

usethis::use_data(sample_dataset_I1, compress = "xz", overwrite = TRUE)
