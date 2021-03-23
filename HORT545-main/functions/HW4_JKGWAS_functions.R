# Function for StaGen HW 4 Question 5
compareGWASnTimes = function(n = 10, X, qtn = 10, CV = CV, PC = PC) {
  X = dplyr::select_if(X, is.numeric);
  # create array to store average for each iteration
  store.truepos = data.frame();
  store.truerate = data.frame();
  store.falsepos = data.frame();
  store.falserate = data.frame();
  
  for(i in 1:n){
    # first, simulate phenotype
    # same parameters as Q2-4
    G2P.sim = G2P(X= X,
                  h2= 0.75,
                  alpha=1,
                  NQTN=qtn,
                  distribution="norm");
    
    # get QTN positions in GD data
    G2P.sim.qtn = G2P.sim$QTN.position;
    
    #get vector of pvals
    #by correlation
    p.list.cor = GWASbyCor(X, G2P.sim$y);
    p.list.cor = as.matrix(t(p.list.cor));
    
    #change NAs to 1s for correct Bonferroni correction
    if(sum(is.na(p.list.cor)!=0)){p.list.cor[is.na(p.list.cor)]=1}
    
    y = as.data.frame(G2P.sim$y); #input needed for JKGLM
    
    #by glm
    p.list.glm = JKGLM(X, y, CV = CV, PC = PC);
    p.list.glm = as.matrix(p.list.glm)
    
    #put results in list
    p.df = list(p.list.cor, p.list.glm);
    names(p.df) = c("COR", "GLM");
    
    #perform test on each list item
    nc = length(p.df);
    for(j in 1:nc){
      p.index = p.df[[j]];
      #identify qtn pvals
      qtn.p = p.index[G2P.sim.qtn];
      #sort vector of pvals as increasing
      #p.index.sort =sort(p.index, decreasing = FALSE, na.last = TRUE);
      p.sig.cor = p.index[p.index < (0.05/length(p.index))];
      
      #find all qtn pvals in total list
      truePos.cor=intersect(p.sig.cor, qtn.p);
      #find all significant pvals not in qtn list
      falsePos.cor=setdiff(p.sig.cor, qtn.p);
      #save the number of true positives for this iteration
      
      store.truepos[i,j] = length(truePos.cor);
      store.truerate[i,j] = length(truePos.cor)/length(p.sig.cor);
      store.falsepos[i,j]=length(falsePos.cor);
      store.falserate[i,j] = length(falsePos.cor)/length(p.sig.cor);
    }
    
  }
  names(store.truepos) = c("COR", "GLM");
  names(store.truerate) = c("COR", "GLM");
  names(store.falsepos) = c("COR", "GLM");
  names(store.falserate) = c("COR", "GLM");
  
  #return a list with 4 lists comparing Cor and GLM
  output = list(store.truepos, store.truerate, store.falsepos, store.falserate);
  names(output) = c("TruePos", "TPR", "FalsePos", "FPR");
  return(output);
  
  #summary(store.truepos);
}


# Function for StaGen HW 4 Extra Credit
simBLINKnTimes = function(n, X, SNP){
  #create arrays to store true positive and time info
  store.truepos = array();
  store.truerate = array();
  store.falsepos = array();
  store.falserate = array();
  BLINK_time= array();
  
  for(i in 1:n){
    #simulate phenotype
    mySim=GAPIT.Phenotype.Simulation(GD=X,
                                     GM=SNP,
                                     h2=.75,
                                     NQTN=10, 
                                     effectunit=.95,
                                     QTNDist="normal",
                                     CV=CV);
    gapit.simY = mySim$Y;
    gapit.qtn = mySim$QTN.position;
    BLINK_start=Sys.time(); #start timing BLINK
    #run BLINK
    myGAPIT_BLINK = GAPIT(Y= gapit.simY,
                          GD=X,
                          GM=SNP, 
                          CV=CV, 
                          PCA.total = 3, 
                          QTN.position = gapit.qtn, 
                          model = "Blink",
                          file.output = F)## GWAS with Blink
    BLINK_end=Sys.time(); # end timing BLINK
    GWAS_result = myGAPIT_BLINK$GWAS ##Extracting the GWAS result from Blink
    GWAS_p=GWAS_result$P.value
    p.index = GWAS_p;
    #identify qtn pvals
    qtn.p = p.index[gapit.qtn];
    #sort vector of pvals as increasing
    #p.index.sort =sort(p.index, decreasing = FALSE, na.last = TRUE);
    p.sig.cor = p.index[p.index < (0.05/length(p.index))];
    
    #find all qtn pvals in total list
    truePos.cor=intersect(p.sig.cor, qtn.p);
    #find all significant pvals not in qtn list
    falsePos.cor=setdiff(p.sig.cor, qtn.p);
    
    #save the number of true positives for this iteration
    
    store.truepos[i] = length(truePos.cor);
    store.truerate[i] = length(truePos.cor)/length(p.sig.cor);
    store.falsepos[i]=length(falsePos.cor);
    store.falserate[i] = length(falsePos.cor)/length(p.sig.cor);
    BLINK_time[i]= BLINK_start- BLINK_end; #get time difference
    
  }
  #output arrays in list
  output = list(store.truepos, store.truerate, store.falsepos, store.falserate, BLINK_time);
  names(output) = c("TruePos", "TPR", "FalsePos", "FPR", "BLINK_Time");
  return(output);
}


# Function for StaGen HW 4 Extra Credit

simJKGWASnTimes = function(n, X, SNP) {
  #create arrays to store true positive and time info
  store.truepos = array();
  store.truerate = array();
  store.falsepos = array();
  store.falserate = array();
  GLM_time= array();
  
  for(i in 1:n){
    #simulate phenotype
    mySim=GAPIT.Phenotype.Simulation(GD=X,
                                     GM=SNP,
                                     h2=.75,
                                     NQTN=10, 
                                     effectunit=.95,
                                     QTNDist="normal",
                                     CV=CV);
    gapit.simY = mySim$Y
    gapit.qtn = mySim$QTN.position
    
    y = as.data.frame(gapit.simY); #needed input for JKGLM
    
    GLM_start=Sys.time(); # start timing GLM
    p.list.glm = JKGLM(X, y, CV = CV, PC = PC);
    GLM_end=Sys.time(); # end timing GLM
    
    p.index = p.list.glm;
    #identify qtn pvals
    qtn.p = p.index[gapit.qtn];
    
    #p.index.sort =sort(p.index, decreasing = FALSE, na.last = TRUE);
    p.sig.cor = p.index[p.index < (0.05/length(p.index))];
    
    #find all qtn pvals in total list
    truePos.cor=intersect(p.sig.cor, qtn.p);
    #find all significant pvals not in qtn list
    falsePos.cor=setdiff(p.sig.cor, qtn.p);
    #save the number of true positives for this iteration
    
    store.truepos[i] = length(truePos.cor);
    store.truerate[i] = length(truePos.cor)/length(p.sig.cor);
    store.falsepos[i]=length(falsePos.cor);
    store.falserate[i] = length(falsePos.cor)/length(p.sig.cor);
    GLM_time[i]= GLM_start-GLM_end
  }
  output = list(store.truepos, store.truerate, store.falsepos, store.falserate, GLM_time);
  names(output) = c("TruePos", "TPR", "FalsePos", "FPR", "GLM_Time");
  return(output);
}