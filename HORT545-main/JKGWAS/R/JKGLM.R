#' Compute SNP p-values for GWAS by GLM 
#' 
#' @description Compute p-values for the association tests between phenotype and SNPs by Generalized Linear Model fitted for the fixed effects of SNP, with option to include Principal Compoenents and Covariates.
#' 
#' @param X Matrix containing genotype data organized with SNPs as columns, individuals as rows
#' @param y Matrix containing SNP positions and chromosomes corresponding to X
#' @param CV Optional argument to include matrix containing covariates. Default is NULL.
#' @param PC Optional argument to include matrix containing principal components such as an object of JKPCA. Default is NULL.
#' 
#' @return Pvals array of length ncol of X for SNPs from fitted GLM
#' @export
#'
## GWAS by GLM Function

JKGLM = function(X, y, CV=NULL, PC=NULL){
  
  X = dplyr::select_if(X, is.numeric); 
  y = dplyr::select_if(y, is.numeric);

  #get dimensions of X
  ncol.X = ncol(X);
  nrow.X = nrow(X);
  
  #get dimensions of y
  n.y=nrow(y);
  
  save.Pvals = matrix(); #create empty matrix to save pvalues
  
  #iterate through all snps
  for(i in 1:ncol.X){
    snp = X[,i];
    if(max(snp) == min(snp)){
      p = 1;
    } else {
      ## No user input for Covariate and Principal Components
      ## only snps are modeled
      if(is.null(CV) && is.null(PC)){
        JK = as.matrix(cbind(1, snp));
      ## User input for Covariates only
      ## snps and covariates are modeled
      } else if(!is.null(CV) && is.null(PC)){
        CV = dplyr::select_if(CV, is.numeric);
        JK = as.matrix(cbind(1,CV,snp));
      ## User input for  PCs only
      ## pcs and snps are modeled
      } else if(is.null(CV) && !is.null(PC)){
        JK = as.matrix(cbind(1,PC,snp));
      ## User input for both PCs and Covariates
      ## snps, PCs, CVs are all modeled
      } else if(!is.null(CV) && !is.null(PC)){
        CV = dplyr::select_if(CV, is.numeric);
        JK = as.matrix(cbind(1,CV,PC,snp));
      }
    
      ## construct GLM
      LHS=as.matrix(t(JK))%*%as.matrix(JK);
      # compute the inverse of the LHS matrix
      inv.LHS=base::solve(LHS);
      RHS=t(JK)%*%as.matrix(y);
      b=inv.LHS%*%RHS;
      yb=JK%*%b;
      e=y-yb;
      ve=sum(e^2)/(n.y-1);
      vt=inv.LHS*ve;
      # compute t test
      t=b/sqrt(diag(vt));
      # get p value from t test
      p=2*(1-pt(abs(t),n.y-2));
    }
      save.Pvals[i] = p[length(p)];
    
  }
  return(save.Pvals);
}


