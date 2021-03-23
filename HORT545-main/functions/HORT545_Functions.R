## function library for HORT 545
## All Credit goes to the lab of Dr. Zhiwu Zhao - http://www.zzlab.net/StaGen/2021/R/

getStaGen545AbbreviationMeaning = function(abbrev){
  if(abbrev == "GBLUP"){
    print("Genomic Best Linear Unbiased Prediction");
  } else if(abbrev == "GBS") {
    print("Genotyping by Sequencing");
  } else if(abbrev == "GWAS"){
    print("Genome-Wide Association Studies");
  } else if(abbrev == "KNN"){
    print("K Nearest Neighbors");
  } else if(abbrev == "LD"){
    print("Linkage Disequilibrium");
  } else if(abbrev == "MAF"){
    print("Minor Allele Frequency");
  } else if(abbrev == "QTL") {
    print("Quantitative Trait Loci");
  } else if(abbrev == "QTN"){
    print("Quantitative Trait Nucleotide");
  } else if(abbrev == "RFLP"){
    print("Restriction Fragment Length Polymorphism");
  } else if(abbrev == "RIL"){
    print("Recombinant Inbred Lines");
  } else if(abbrev == "SNP") {
    print("Single Nucleotide Polymorphism");
  } else if(abbrev == "SSR") {
    print("Simple Sequence Repeats");
  } else {
    print("Abbreviation input is not recognized. Make sure to use all-caps.");
  }
}

StochasticImpute = function(X){
  n=nrow(X)
  m=ncol(X)
  
  ## AA is coded as 2
  ## AB is coded as 1 - 1 is always heterozygous
  ## BB is coded as 0, so when we are looking for frequency of A, we are looking only for A not AA
  
  fn=colSums(X, na.rm=T) # sum of genotypes for all individuals 
  fc=colSums(floor(X/3+1),na.rm=T) #count number of non missing individuals - basically just assigns a 1 to every existing observation
  fa=fn/(2*fc) #Frequency of allele A - there are 2fc total numbers of positions - aka 2 possible allelles per ob
  for(i in 1:m){
    index.a=runif(n)<fa[i] # puts true for values that are less than the frequency of allelle A in column i
    index.na=is.na(X[,i]) # puts true for all missing values
    index.m2=index.a  &  index.na # cumulative of missing values or lower than freq of A
    index.m0=!index.a  &  index.na # values that are greater than the frequency of A and missing
    X[index.m2,i]=2 #assign 2 to all missing values lower than freq A, bc we assume they are A
    X[index.m0,i]=0 # assign 0 to all missing values greater than freq A, bc we assume they are B
  }
  return(X)} 

G2P=function(X,h2,alpha,NQTN,distribution){
  n=nrow(X)
  m=ncol(X)
  #Sampling QTN
  QTN.position=sample(m,NQTN,replace=F)
  SNPQ=as.matrix(X[,QTN.position])
  QTN.position
  #QTN effects
  if(distribution=="norm")
  {addeffect=rnorm(NQTN,0,1)
  }else
  {addeffect=alpha^(1:NQTN)}
  #Simulate phenotype
  effect=SNPQ%*%addeffect
  effectvar=var(effect)
  residualvar=(effectvar-h2*effectvar)/h2
  residual=rnorm(n,0,sqrt(residualvar))
  y=effect+residual
  return(list(addeffect = addeffect, y=y, add = effect, residual = residual, QTN.position=QTN.position, SNPQ=SNPQ))
}


#To perform GWAS with correlation
#Input: y-vector of phenotype, X matrix of numeric genotype with rows as individuals and SNPs as column
#Output: vector of probability for SNPs in same order
GWASbyCor=function(X,y){
  n=nrow(X)
  r=cor(y,X) #calculate correlation between y phenotype and each snp
  n=nrow(X)
  t=r/sqrt((1-r^2)/(n-2)) #convert the r to the t statistic
  p=2*(1-pt(abs(t),n-2))
  zeros=p==0
  p[zeros]=1e-10
  return(p)}

cort=function(n=10000,df=100){
  z=replicate(n,{
    x=rnorm(df+2)
    y=rnorm(df+2)
    r=cor(x,y)
    t=r/sqrt((1-r^2)/(df))
  })
  return(z)}

