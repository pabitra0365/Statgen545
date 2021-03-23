#' Manhattan Plot for GWAS Visualization
#'
#'
#' @description Visualize the GWAS by GLM results by Manhattan plot. User can also specify QTN.
#'
#'
#' @param Pvals Input vector of Pvals such as an object returned by the JKGLM function
#' @param SNP Input matrix containing SNP location. Must include columns Position and Chromosome
#' @param sig.cutoff Significance threshold for visualization. If NULL, then uses bonferroni correction with alpha = 0.05
#' @param QTN Vector of QTN positions that is provided by the user to highlight position in the Manhattan plot. If NULL, then QTN will not be identified.
#'
#' @return Manhattan plot with user inputs.
#' @export
#'

JKManhattan = function(Pvals, SNP, sigcutoff = NULL, QTN = NULL){
  if(is.null(sigcutoff)){
    sigcutoff = 0.05/length(Pvals); #default Bonferonni correction
  }
  
  #colors to differentiate chromosomes
  nChrom =length(unique(SNP$Chromosome));
  color.vector <- rep(c("#3c1642","#086375","#1DD3B0","#e01a4f"),nChrom);
  m=length(Pvals);
  #plot the positions with -log of the pvalues for better visualization
  plot(seq(1:m),t(-log10(Pvals)),
       col=color.vector[SNP[,2]], #2 is the chromosome column
       xlab="SNP Positions",
       ylab = "-log(Pvalues)",
       main ="Manhattan Plot")
  #plot a line at the significance cutoff
  abline(h=(-log10(sigcutoff)), lwd=2, col = "black")
  # if user provides QTN, identify positions with black points and vertical lines
  if(!is.null(QTN)){
    points(QTN,-log10(Pvals[QTN]), pch=19)
    abline(v = QTN, lty = 2, lwd = 1, col = "black")
  }
}



