#'  Principle Component Analysis
#'
#'
#' @description Removes PCs that are linearly dependent with the given covariates and also user can specify how many PCs to chooose as co-factors
#'
#'
#' @param npc Number of Principle Components (PCs) that are specified by the user 
#' @param X     Markers data in the form n by m with n number of individuals and m number of markers
#' @param CV     Covariates matrix in the form n by t with n number of individuals and t number of co-variates
#'
#' @export
#' @return Principle Component Analysis
#' 

JKPCA = function(X,CV = NULL, npc = 5){
  X = dplyr::select_if(X, is.numeric);
  jkpca = prcomp(X); #perform pca with prcomp() defaults
  pc = jkpca$x; #save principal component information
  pc = pc[,c(1:npc)]; #filter for number of desired pc
  pc.out = pc; #if no CV, pca is computed and returned
  
  if(!is.null(CV)){
    
  CV = dplyr::select_if(CV, is.numeric);
  JKmatrix = as.matrix(cbind(CV, pc));
  nc = ncol(CV)+1;
  #determine linear dependence by matrix rank with successive removal of PCs
  rankByRemoved <- sapply(nc:ncol(JKmatrix), function (x) qr(JKmatrix[,-x])$rank);
  #identify columns which, if removed, do not change the rank and are therefore dependent
  removeIDs = which(rankByRemoved == max(rankByRemoved));
  
  #filter these columns from the PC matrix
  pc.out = pc[,-c(removeIDs)];
  
  }
  
  return(pc.out);
}

#https://stats.stackexchange.com/questions/16327/testing-for-linear-dependence-among-the-columns-of-a-matrix