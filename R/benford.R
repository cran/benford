####################################################
benford <- function(x, plot=FALSE, mode=1) {
  x <- x[which(!is.na(x))]
  dist <- abs(x)
  if(max(dist)/min(dist) < 1000) {
    print("Benford's analysis perform better when the analyzed data set ranges for more than
three orders of magnitude than when it ranges for just one or two orders of magnitude (Fewster, 2009).")
  }

  if(mode == 1) {
  first <- as.integer(sapply(strsplit(as.character(dist),""), `[[`, 1))
  n <- c(1:9)
  benford_dist <- log10(1+1/n)
  first <- first[which(first!=0)]
  xlab <- "First digit"
  } else if(mode == 2) {
    dist <- dist[dist>=10]
    first <- as.integer(sapply(strsplit(as.character(dist),""), `[[`, 2))
    n <- c(0:9)
    benford_dist <- NULL
    for(d in c(0:9)){
      benford_dist[d+1] <- sum(log10(1+1/(10*c(1:9)+d)))
    }
    xlab <- "Second digit"
    } else if(mode == 12) {
        dist <- dist[dist>=10]
        first <- as.integer(paste(sapply(strsplit(as.character(dist),""), `[[`, 1),
                              sapply(strsplit(as.character(dist),""), `[[`, 2),
                              sep=""))
        n <- c(10:99)
        benford_dist <- log10(1+1/n)
        first <- first[which(first!=0)]
        xlab <- "first-two digits"
        } else if(mode == 123){
      dist <- dist[dist>=100]
      first <- as.integer(paste(sapply(strsplit(as.character(dist),""), `[[`, 1),
                              sapply(strsplit(as.character(dist),""), `[[`, 2),
                              sapply(strsplit(as.character(dist),""), `[[`, 3),
                              sep=""))
    n <- c(100:999)
    benford_dist <- log10(1+1/n)
    first <- first[which(first!=0)]
    xlab <- "first-three digits"
    } else {
    print("Mode must be 1 for first digit analysis, 2 for second digit, 3 for
    the third digit or 12 for first-two and 123 for first-three digits analysis")
    }

  ### Compute RMSD, Likelihood and Chi-square for each country
  benford.data <- NULL
  benford_hist <- hist(first, breaks=c(n-.5,max(n+.5)),
                       right=F, plot=F)
  benford_counts <- benford_hist$counts / length(first)
  rmsd <- sqrt(mean((benford_counts - benford_dist)^2, na.mr=T))
  chisq <- chisq.test(matrix(c(benford_hist$counts, benford_dist * length(first)),
                             ncol=2, byrow=F))
  likelihood <- sum(dnorm(benford_counts, mean=benford_dist, sd=sd(benford_counts), log=T))

  chart <- matrix(c(n=1:length(benford_dist), data=benford_hist$counts, benford=benford_dist*length(first)),
                  ncol=3)

  if(plot == T) {
    if(max(chart[,3]) >= max(chart[,2])){
      ymax <- max(chart[,3])*1.05 } else {
        ymax <- max(chart[,2])*1.05
      }
    plot(chart[,1], chart[,2], type="h", ylim=c(0,ymax),
         main=paste("p =", round(chisq[[3]],4)), xlab=xlab,
         ylab="counts")
    points(chart[,1], chart[,3], col="red")
    lines(chart[,1], chart[,3], col="red")
  }

  exit <- list(c(p=chisq[[3]], RMSD=rmsd, LogLikelihood=likelihood), chart)
  return(exit)
}

