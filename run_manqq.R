library(argparse)
library(zoo)

source("~sh29/repos/man_qq/manqq_functions.R")

# create parser object
parser <- ArgumentParser(description="A program to plot Manhattan and QQ plots")


parser$add_argument("--chr-col", 
                    type="character",
                    default="chr",
                    help="The column NAME for the chromosome column, default chr",
                    metavar="[character]")

parser$add_argument("--pval-col", 
                    type="character",
                    default="pval",
                    help="The column NAME for the chromosome column, default pval",
                    metavar="[character]")

parser$add_argument("--pos-col", 
                    type="character",
                    default="pos",
                    help="The column NAME for the chromosome column, default pos",
                    metavar="[character]")

parser$add_argument("--a1", 
                    type="character",
                    default="pos",
                    help="The column NAME for the effect allele column, default a1",
                    metavar="[character]")

parser$add_argument("--a2", 
                    type="character",
                    default="pos",
                    help="The column NAME for the non-effect column, default a2",
                    metavar="[character]")

#parser$add_argument("--sig-thresh-line", 
#                    type="double",
#                    default=-1.0,
#                    help="The significance threshold for the line",
#                    metavar="[DOUBLE]")
#
#parser$add_argument("--title", 
#                    type="character",
#                    default="",
#                    help="An optional text title to add to each plot",
#                    metavar="[character]")

parser$add_argument("--build", type="integer",default=38,
                    help="The genome build the positions refer to",
                    metavar="[integer]")

# The input file used to create the plots
parser$add_argument("infile", nargs=1, help="Input file name, must be gzip file")

# The input file used to create the plots
parser$add_argument("outfile", nargs=1, help="Output file name (with no file extension)")

args <- parser$parse_args()
#args=list()
#args$a1="allele1"
#args$a2="allele0"
#args$chr="chr"
#args$pos="ps"
#args$pval="p_score"

## do the job
library(data.table)

readcmd=paste("zcat ", args$infile, sep=" ")

outqq=paste(args$outfile, ".qq.pdf", sep="")
outman=paste(args$outfile, ".man.pdf", sep="")

d=fread(readcmd, select=c(args$chr,args$pos,args$a1,args$a2,args$pval))
#d=fread(readcmd, select=c(1,3,5,6,14))

## QQ PLOT
ret=qqplot(d[,args$pval,with=FALSE][[1]])

nn=nrow(d)
upper=rep(NA, nrow(ret))
k=0;for(i in ret$order){k=k+1;upper[k]=qbeta(0.95, i, nn-i+1)}
lower=rep(NA, nrow(ret))
k=0;for(i in ret$order){k=k+1;lower[k]=qbeta(0.05, i, nn-i+1)}

#pdf(outqq)
#plot(ret$x, ret$y, pch=20, col="darkslategray", type="n", xlab="Expected quantiles", ylab="Observed quantiles")
#xx =  -log10((ret$order)/(nn+1))
#polygon(c(xx, rev(xx)), c(-log10(upper), -log10(rev(lower))), border=NA, col="gray80")
#lines(xx, -log10(upper), col="gray", lty=2, lwd=2)
#lines(xx, -log10(lower), col="gray", lty=2, lwd=2)
#lambdavalue=lambdaCalc(d[,args$pval,with=FALSE][[1]])
#text(substitute(paste(lambda, "=", lambdaval), list(lambdaval=lambdavalue)), x=1, y=max(ret$y)-1, cex=1.5)
#abline(a=0, b=1, col="firebrick", lwd=2)
#points(ret$x, ret$y, pch=20, col="dodgerblue4")
#dev.off()


pdf(outman, width=10, height=6)
retm=mhp(d[,args$chr,with=FALSE][[1]], d[,args$pos,with=FALSE][[1]], d[,args$pval,with=FALSE][[1]])
print("HELLO")
peaks=get_peaks_to_annotate(retm)

context=apply(peaks, 1, function(x){
  u=unlist(get_variant_context(as.numeric(x["chr"]), as.numeric(x["ps"]), x["a1"], x["a2"]))
  if(length(u)<3){u[3]="unknown"};
  return(u);
  })

context=as.data.frame(t(context))

colnames(context)=c("gene", "distance", "consequence")
context$distance=as.numeric(as.character(context$distance))
peaks=cbind(peaks, context)
peaks$truelabels=peaks$gene
peaks$truelabels[peaks$dist>0]=paste(peaks$gene, paste(" (", ceiling(peaks$distance/1000), "kbp)", sep=""))peaks$pch=15
peaks$col="forestgreen"
lof=c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant")
high=c("stop_lost", "start_lost", "transcript_amplification")
exonic=c("inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant")
low=c("splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant")
intronic=c("intron_variant")
intergenic=c("intergenic_variant")
peaks$pch[peaks$consequence %in% lof]=4
peaks$col[peaks$consequence %in% lof]="firebrick"
peaks$pch[peaks$consequence %in% high]=17
peaks$col[peaks$consequence %in% high]="orange"
peaks$pch[peaks$consequence %in% exonic]=25
peaks$col[peaks$consequence %in% exonic]="goldenrod"
peaks$pch[peaks$consequence %in% low]=25
peaks$col[peaks$consequence %in% low]="brown"
peaks$pch[peaks$consequence %in% intronic]=18
peaks$col[peaks$consequence %in% intronic]="brown"
peaks$pch[peaks$consequence %in% intergenic]=19
peaks$col[peaks$consequence %in% intergenic]="darkgray"

plot_manhattan(retm, peaks)
dev.off()


