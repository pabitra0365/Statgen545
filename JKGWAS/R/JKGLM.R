#GWAS by GLM function

JKGLM = function(X,y,CV=Null,PC=Null){
  ##First build our regression matrix based on inputs
  ##Genomic data = X,which is mandatory

  ncol.X=ncol(X):
  nrow.X= nrow(X):
    save.Pvals =matrix():
    for (i in 1:ncol.X){
      snp=X[,i]
      if(max(snp)==min(snp)){
        p=1}else{
          if(is.null(CV) & is.null(PC)){
            JK=matrix(cbind(1,snp))}

         else if(!is.null(CV) & is.null(PC)){
            JK=as.matrix(cbind(1,CV,snp))}

         else if(is.null(CV) & !is.null(PC)){
            JK=as.matrix(cbind(1,PC,snp))}

         else if(!is.null(CV) & !is.null(PC)){
           ##exclude dependency
           for(i in 1:ncol.PC){
             current.pc=PC[,i]
             df=as.data.frame(cbind(y,CV,current.pc))
             names(df)=
               df.lm=lm(y~CV+)
           }
            JK=as.matrix(cbind(1,CV,PC,snp))}

          LHS=t(JK)%*%JK
          inv.LHS=solve(LHS)
          RHS=t(JK)%*%y
          b=inv.LHS%*%RHS
          yb=JK%*%b
          e=y-yb
          n.y=length(y)
          ve=sum(e^2)/(n.y-1)
          vt=inv*ve
          t=b/sqrt(diag(vt))
          p=2*(1-pt(abs(t),n.y-2))
        } #end of testing variation
      save.Pvals[i]=p[length(p)]
    } #end of looping for markers

    return(save.Pvals)
  }
