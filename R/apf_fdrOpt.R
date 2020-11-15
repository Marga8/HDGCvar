#' @title Wrapper around APFr
#' @keywords internal
#' @description Wrapper around APFr \code{ apf_plot} to get the desired threshold on the p values such that the FDR is at most \code{ FDR_max}
#' and APF is closest to \code{ APF_min}
#' @param mat_input the input matrix
#' @param FDR_max the maximum false discovery rate desired
#' @param APF_min the minimum average power function desired
#' @param verbose logical for verbose console
#' @return the desired constrained opt result
apf_fdrOpt<-function(mat_input,FDR_max,APF_min,verbose=F){
  if(is.null(FDR_max)){
    index2<-which.min(abs(mat_input[,4]-APF_min))
    result<-mat_input$Gamma[index2]
  }
  else if(is.null(APF_min)){
    index1<-which.min(abs(mat_input[,3]-FDR_max))
    result<-mat_input$Gamma[index1]
  }
  else{
    index1<-which.min(abs(mat_input[,3]-FDR_max))
    index2<-which.min(abs(mat_input[,4]-APF_min))
    if(index1==index2){
      result<-mat_input$Gamma[index1]
    }
    if(index1!=index2){
      index3<-ceiling((index1+index2)/2)
      result<-mat_input$Gamma[index3]
    }
  }
  if(verbose==T){
    Vresult<-list(result,paste("FDR=",mat_input[,3][index1]),paste("APF=",mat_input[,4][index2]))
    return(Vresult)
  }
  if(verbose==F){
    return(result)
  }
}
