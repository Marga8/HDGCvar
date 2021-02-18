#' @title Plot High-Dimensional Granger causality Networks
#'
#' @param Comb  output from: \code{\link{HDGC_VAR_all_I0}}, \code{\link{HDGC_VAR_multiple_pairs_I0}}, \code{\link{HDGC_VAR_all}}, \code{\link{HDGC_VAR_multiple_pairs}}, \code{\link{HDGC_HVAR_all}} or \code{\link{HDGC_HVAR_multiple_pairs}}
#' @param Stat_type either FS_cor (default), Asymp or Asymp_Robust respectively for F-stat small sample correction, standard Chi square test, standard chi square test with heteroscedasticity correction
#' @param alpha the desired probability of type one error, default is 0.01.
#' @param multip_corr A list: first element is logical, if TRUE a multiple testing correction using \code{\link[stats]{p.adjust}} is used. The second
#'                    element of the list define the p.adjust.method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'                    "fdr", "none")). If the second element gets the name "APF_FDR" then \code{\link[APFr]{apf_fdr}} is called which uses empirical Bayes is called and a third and fourth
#'                    elements in the mutip_corr list are required: gamm=c(a,b,c) requires a min (a), max (b) and step length (c) values to be set for the threshold on the p_values,
#'                    \code{\link[APFr]{apf_fdr}} requires one or two values: either (NULL,value) or (value,NULL) if one wants to have specified amount of average power (fdr) no matter fdr (average power).
#'                    If both (value,value) are given, the calculated threshold will find the closest combination to both apf and fdr desired. The last element
#'                    of the list is logical: verbose=TRUE if one wants to know how much apf/fdr the testing has.
#' @param ... all parameters for the network plot: see example and \code{\link[igraph]{graph_from_adjacency_matrix}} documentation.
#' @param cluster A list: first element is logical, if TRUE a cluster plot using \code{\link[igraph]{cluster_edge_betweenness}} is plotted.
#'                Other elements are respectively: vertex.size, vertex.label.color,vertex.label.cex, vertex.label.dist, edge.curved (see \code{\link[igraph]{graph_from_adjacency_matrix}} for details).
#' @return a \code{\link[igraph]{graph_from_adjacency_matrix}} network
#' @export
#' @importFrom igraph graph_from_adjacency_matrix as.undirected cluster_edge_betweenness V V<-
#' @importFrom APFr apf_fdr apf_plot
#' @importFrom stats p.adjust
#' @importFrom graphics par
#' @examples \dontrun{Plot_GC_all(Comb, "FS_cor",alpha=0.01,multip_corr=list(F), directed=T, layout.circle}
#' @references Hecq, A., Margaritella, L., Smeekes, S., "Inference in Non Stationary High Dimensional VARs" (2020, check the latest version at https://sites.google.com/view/luca-margaritella )
#' @references Hecq, A., Margaritella, L., Smeekes, S., "Granger Causality Testing in High-Dimensional VARs: a Post-Double-Selection Procedure." arXiv preprint arXiv:1902.10991 (2019).
#' @references Newman, Mark EJ, and Michelle Girvan. "Finding and evaluating community structure in networks." Physical review E 69.2 (2004): 026113.
#' @references Quatto, Piero, et al. "Brain networks construction using Bayes FDR and average power function." Statistical Methods in Medical Research 29.3 (2020): 866-878.
Plot_GC_all<-function(Comb,Stat_type="FS_cor", alpha=0.01, multip_corr=list(F,"bonferroni",gamm = c(1e-04, 0.1, 0.001),fdr.apf=c(0.05,0.6),verb=F),...,
                      cluster=list(F,10,"black",0.51, 1, 0)){
  if(Stat_type=="FS_cor"){
    #build adjaciency matrix
    input<-as.matrix(Comb[["tests"]][,,2,2])}#p_values
    input<-t(input)
  if(Stat_type=="Asymp"){
    #build adjaciency matrix
    input<-as.matrix(Comb[["tests"]][,,2,1])}#p_values
    input<-t(input)
  if(Stat_type=="Asymp_Robust"){
    #build adjaciency matrix
    input<-as.matrix(Comb[["tests"]][,,2,3])}#p_values
    input<-t(input)
  if(multip_corr[[1]]==F){
    input[input < alpha] <- 1 #put =1 values < alpha
    input[is.na(input)] <- 0 #put =0 the diagonal
    input[input != 1] <- 0 #put =0 values > alpha
    network=graph_from_adjacency_matrix(input, mode='directed',diag=F,add.rownames = TRUE )
    V(network)$label = rownames(input)
  }
  if(multip_corr[[1]]==T){
    if(multip_corr[[2]]=="APF_FDR"){
      input_without_NA<-as.vector(input)
      input_without_NA<-input_without_NA[!is.na(input_without_NA)] #take away NAs
      APF_lst<-apf_fdr(data=input_without_NA, type = "pvl", lobs =300 ,
                       seed = 123, gamm = c(multip_corr[[3]]))
      plotAPF_FDR<-apf_plot(APF_lst, tab = TRUE, APF_inf = 0.5, FDR_sup = 0.1)
      alpha_threshold<-apf_fdrOpt(plotAPF_FDR,FDR_max=multip_corr$fdr.apf[1],APF_min=multip_corr$fdr.apf[2],verbose=F)
      if(multip_corr$verb==T){
        apf_fdrOpt(plotAPF_FDR,FDR_max=multip_corr$fdr.apf[1],APF_min=multip_corr$fdr.apf[2],verbose=T)
      }
      input[input < alpha_threshold] <- 1 #put =1 values < alpha
      input[is.na(input)] <- 0 #put =0 the diagonal
      input[input != 1] <- 0 #put =0 values > alpha
      network=graph_from_adjacency_matrix(input, mode='directed',diag=F,add.rownames = TRUE )
      V(network)$label = rownames(input)
    }
    else if (multip_corr[[2]]!="APF_FDR"){
      adj_input<-p.adjust(as.vector(input),method =multip_corr[[2]]) #adjust p-values for multiple testing
      adj_pval_mat<-matrix(adj_input,nrow =nrow(input) ,ncol =ncol(input),byrow = F) #put them back in a matrix
      adj_pval_mat[adj_pval_mat < alpha] <- 1 #put =1 values < alpha
      adj_pval_mat[is.na(adj_pval_mat)] <- 0 #put =0 the diagonal
      adj_pval_mat[adj_pval_mat != 1] <- 0 #put =0 values > alpha
      rownames(adj_pval_mat)<-rownames(input)
      colnames(adj_pval_mat)<-colnames(input)
      network=graph_from_adjacency_matrix(adj_pval_mat, mode='directed',diag=F,add.rownames = TRUE )
      V(network)$label = rownames(adj_pval_mat)
    }
  }
  par(mfrow=c(1,1))
  plot(network,...)
  if(cluster[[1]]==TRUE){
    ## make graph undirected ##
    net.sym <- as.undirected(network, mode= "collapse")
    ceb <- cluster_edge_betweenness(net.sym,directed=T)
    plot(ceb, net.sym, vertex.size=cluster[[2]],vertex.label.color=cluster[[3]],
         vertex.label.cex=cluster[[4]], vertex.label.dist=cluster[[5]], edge.curved=cluster[[6]] )
    #dendPlot(ceb, mode="hclust")
  }
}

