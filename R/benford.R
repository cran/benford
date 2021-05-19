####################################################
benford <- function(x, plot=FALSE) {
  x <- x[which(!is.na(x))]
  dist <- abs(x)
  if(max(dist)/min(dist) < 1000) {
    print("Benford's analysis perform better when the analyzed data set ranges for more than
three orders of magnitude than when it ranges for just one or two orders of magnitude (Fewster, 2009).")
  }
  first <- as.integer(sapply(strsplit(as.character(dist),""), `[[`, 1))
  first <- first[which(first!=0)]

  ### Compute RMSD, Likelihood and Chi-square for each country
  benford_dist <- c(.301,.176,.125,.097,.079,.067,.058,.051,.046)
  benford.data <- NULL
  benford_hist <- hist(first, breaks=c(.5,1.5,2.3,3.5,4.5,5.5,6.5,7.5,8.5,9.5),
                       right=F, plot=F)
  benford_counts <- benford_hist$counts / length(first)
  rmsd <- sqrt(mean((benford_counts - benford_dist)^2, na.mr=T))
  chisq <- chisq.test(matrix(c(benford_hist$counts, benford_dist * length(first)),
                             ncol=2, byrow=F))
  likelihood <- sum(dnorm(benford_counts, mean=benford_dist, sd=sd(benford_counts), log=T))

  chart <- matrix(c(n=1:9, data=benford_hist$counts, benford=benford_dist*length(first)),
                  ncol=3)

  if(plot == T) {
    plot(chart[,1], chart[,2], type="h", ylim=c(0,length(first)/2.5),
         main=paste("p =", round(chisq[[3]],4)), xlab="First digit",
         ylab="counts")
    points(chart[,1], chart[,3], col="red")
    lines(chart[,1], chart[,3], col="red")
  }

  exit <- list(c(p=chisq[[3]], RMSD=rmsd, LogLikelihood=likelihood), chart)
  return(exit)
}



