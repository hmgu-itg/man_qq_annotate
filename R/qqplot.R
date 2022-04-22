save_qqplot = function(outfile, pvalue, image.type = 'png', title="") {
  qqfile = paste0(outfile, ".qq.", image.type)
  if(image.type=="pdf") {
    pdf(qqfile)
  } else if(image.type=="png") {
    png(qqfile)
  }
  lambdavalue = fastqq(pvalue, title=title)
  # Save lambda value to a separate file
  lambdafile=paste0(outfile, ".lambda.txt")
  cat(lambdavalue), file=lambdafile, sep="\n")

  dev.off()
  return(NULL)
}

compute_qqplot = function(data, X_GRID=800, Y_GRID=800){
  ### QQ plot
  ## Expects data to be a vector of p values
  print("entering qq function")
  obspval <- sort(data)
  nrows=length(data)
  logobspval <- -(log10(obspval))
  exppval <- c(1:length(obspval))
  logexppval <- -(log10( (exppval-0.5)/length(exppval)))
  obsmax <- trunc(max(logobspval))+1
  expmax <- trunc(max(logexppval))+1

  yres=(max(logobspval)-min(logobspval))/Y_GRID
  xres=(max(logexppval)-min(logexppval))/X_GRID
  ymax=max(logobspval)
  xmax=max(logexppval)
  index=1
  indey=1
  newx=rep(NA, X_GRID)
  newy=rep(NA, X_GRID)
  ord=rep(NA, X_GRID)
  lowx=xmax
  lowy=ymax
  i=1
  while(index<nrows){
    lowx=lowx-xres;
    lowy=lowy-yres;
    before=index;
    while(logexppval[index]>=lowx & logobspval[index]>=lowy){
      index=index+1
      if (index>nrows){break;}
    }
    lowx=logexppval[index]
    lowy=logobspval[index]

    if(before==index){next;}
    newx[i]=logexppval[before]-0.5*xres
    newy[i]=logobspval[before]-0.5*yres
    ord[i]=before
    i=i+1;
  }
  newx=zoo::na.trim(newx)
  newy=zoo::na.trim(newy)
  ord=zoo::na.trim(ord)
  return(data.frame(x=newx, y=newy, order=ord))
}

#' Display a QQ-plot of a p-value vector from a GWAS
#'
#' @param pvalue A vector of p-values
#' @param X_GRID The horizontal resolution of the grid used to simplify the plot. 800 is good for most cases. 
#' @param Y_GRID The vertical resolution of the grid used to simplify the plot. 800 is good for most cases. 
#' @return NULL
#' @examples
#' \donttest{
#' library(data.table)
#' mygwas=read.table("GWAS.GEMMA.assoc.txt.gz")
#' qqplot(mygwas$P_SCORE)
#' }
#' @export
fastqq = function(pvalue, X_GRID=800, Y_GRID=800, title="") {
  ret = compute_qqplot(pvalue)
  nn = length(pvalue)
  upper=rep(NA, nrow(ret))
  lower=rep(NA, nrow(ret))
  k=0
  for (i in ret$order){
    k=k+1
    upper[k]=qbeta(0.95, i, nn-i+1)
    lower[k]=qbeta(0.05, i, nn-i+1)
  }

  plot(ret$x, ret$y, pch=20, col="darkslategray", type="n", xlab="Expected quantiles", ylab="Observed quantiles", main=title)
  xx = -log10((ret$order)/(nn+1))
  polygon(c(xx, rev(xx)), c(-log10(upper), -log10(rev(lower))), border=NA, col="gray80")
  lines(xx, -log10(upper), col="gray", lty=2, lwd=2)
  lines(xx, -log10(lower), col="gray", lty=2, lwd=2)
  lambdavalue = lambdaCalc(pvalue)


  text(substitute(paste(lambda, "=", lambdaval), list(lambdaval=lambdavalue)), x=1, y=max(ret$y)-1, cex=1.5)
  abline(a=0, b=1, col="firebrick", lwd=2)
  points(ret$x, ret$y, pch=20, col="dodgerblue4")
  return(lambdavalue)
}

lambdaCalc = function(pval, round=NULL) {
  lambda = median(qchisq(pval, 1, lower.tail = F), na.rm=T) / qchisq(0.5, 1)

  if (is.null(round)==FALSE) {
    lambda=round(lambda,round)
  }
  print(paste("lambda is", lambda))
  return(lambda)
}
