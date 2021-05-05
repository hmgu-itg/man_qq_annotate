#!/usr/bin/env Rscript
# Source in packages and functions
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(data.table))


#print(getSrcDirectory(function(x) {x}))
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
print(paste("Sourcing",script.name))

suppressPackageStartupMessages(source(paste(script.basename, "manqq_functions.R", sep="/")))

# run ./run_manqq.R --help or ./run_manqq -h to display help message

# create parser object
parser <- ArgumentParser(description="A program to plot Manhattan and QQ plots")


parser$add_argument("--chr-col",
                    type="character",
                    default="chr",
                    help="The column NAME for the chromosome column, default chr",
                    metavar="[character]")

parser$add_argument("--pval-col",
                    type="character",
                    default="p_score",
                    help="The column NAME for the chromosome column, default p_score",
                    metavar="[character]")

parser$add_argument("--pos-col",
                    type="character",
                    default="ps",
                    help="The column NAME for the chromosome column, default ps",
                    metavar="[character]")

parser$add_argument("--a1",
                    type="character",
                    default="allele1",
                    help="The column NAME for the effect allele column, default allele1",
                    metavar="[character]")

parser$add_argument("--a2",
                    type="character",
                    default="allele0",
                    help="The column NAME for the non-effect column, default allele0",
                    metavar="[character]")

parser$add_argument("--build",
                    type="integer",
                    default=38,
                    help="The genome build the positions refer to",
                    metavar="[integer]")

parser$add_argument("--image",
                    type="character",
                    default="pdf",
                    help="The filetype to save plots to (png or pdf)",
                    metavar="[character]")

parser$add_argument("--af-col",
                    type="character",
                    default="af",
                    help="The column NAME for the allele frequency column, default af",
                    metavar="[character]")

parser$add_argument("--maf-filter",
                    type="double",
                    default=-0.0,
                    help="The significance threshold for MAF filter, default 0.0.",
                    metavar="[double]")

parser$add_argument("--sig",
                    type="double",
                    default=5e-8,
                    help="The significance threshold to use for peak annotation",
                    metavar="[double]")

parser$add_argument("--maxpeaks",
                    type="integer",
                    default=30,
                    help="The maximum number of peaks to annotate",
                    metavar="[integer]")

parser$add_argument("--no-qq",
                    action="store_true",
                    help="Don't plot QQ.")

parser$add_argument("--no-man",
                    action="store_true",
                    help="Don't plot Manhattan.")

parser$add_argument("--no-annot",
                    action="store_true",
                    help="Disable peak annotation even if peaks are present.")

parser$add_argument("--no-distance",
                    action="store_true",
                    help="Don't add very useful distance to gene info.")

parser$add_argument("--man-height",
                    type="integer",
                    default=6,
                    help="Force height of Manhattan in inches. Can have unpredictable consequences (some of which you may regret).",
                    metavar="[integer]")

parser$add_argument("--upper-margin",
                    type="double",
                    default=2.0,
                    help="Y limit of Manhattan plot in units of maximum data points. Even more unpredictable than the above.",
                    metavar="[double]")

parser$add_argument("--annot-cex",
                    type="double",
                    default=1.1,
                    help="Size factor for annotations.",
                    metavar="[double]")

parser$add_argument("--axes-cex",
                    type="double",
                    default=1.3,
                    help="Size factor for axes and labels.",
                    metavar="[double]")



parser$add_argument("--ylim",
                    type="integer",
                    default=-1.0,
                    help="The y-axis limit (-log10(p))",
                    metavar="[double]")


#parser$add_argument("--title",
#                    type="character",
#                    default="",
#                    help="An optional text title to add to each plot",
#                    metavar="[character]")



# The input file used to create the plots
parser$add_argument("infile", nargs=1, help="Input file name, must be gzip file")

# The input file used to create the plots
parser$add_argument("outfile", nargs=1, help="Output file name (with no file extension)")

args=parser$parse_args()

readcmd=paste("zcat ", args$infile, sep=" ")

d=fread(readcmd, select=c(args$chr,args$pos,args$a1,args$a2,args$pval,args$af))
#d=fread(readcmd, select=c(1,3,5,6,14))
setcolorder(d, c(args$chr,args$pos,args$a1,args$a2,args$pval,args$af))

setnames(d, c("chr","pos","a1","a2","p", "af"))

if (args$maf>0.0){
    d=d[af>=args$maf]
}

d=d[!(is.na(p)) & p!=0]

## QQ PLOT
if(! args$no_qq){
ret=qqplot(d[,p])

nn=nrow(d)
upper=rep(NA, nrow(ret))
k=0;for(i in ret$order){k=k+1;upper[k]=qbeta(0.95, i, nn-i+1)}
lower=rep(NA, nrow(ret))
k=0;for(i in ret$order){k=k+1;lower[k]=qbeta(0.05, i, nn-i+1)}

if(args$image=="pdf") {
    qqfile = paste(args$outfile, ".qq.pdf", sep="")
    pdf(qqfile)
} else if(args$image=="png") {
    qqfile = paste(args$outfile, ".qq.png", sep="")
    png(qqfile)
}
plot(ret$x, ret$y, pch=20, col="darkslategray", type="n", xlab="Expected quantiles", ylab="Observed quantiles")
xx =  -log10((ret$order)/(nn+1))
polygon(c(xx, rev(xx)), c(-log10(upper), -log10(rev(lower))), border=NA, col="gray80")
lines(xx, -log10(upper), col="gray", lty=2, lwd=2)
lines(xx, -log10(lower), col="gray", lty=2, lwd=2)
lambdavalue=lambdaCalc(d[,p])
# Save lambda value to a separate file
lambdafile=paste0(args$outfile, ".lambda.txt")
cat(paste(qqfile, lambdavalue, sep='\t'), file=lambdafile, sep="\t")

text(substitute(paste(lambda, "=", lambdaval), list(lambdaval=lambdavalue)), x=1, y=max(ret$y)-1, cex=1.5)
abline(a=0, b=1, col="firebrick", lwd=2)
points(ret$x, ret$y, pch=20, col="dodgerblue4")
dev.off()
}

if(! args$no_man){
## MANHATTAN PLOT
if(args$image=="pdf") {
    pdf(paste(args$outfile, ".man.pdf", sep=""), width=10, height=args$man_height)
} else if(args$image=="png") {
    png(paste(args$outfile, ".man.png", sep=""), width=10, height=args$man_height, units="in",res=300)
}

retm=mhp(d[,chr], d[,pos], d[,p],signif=args$sig)
print("Finding peaks...")
peaks=get_peaks_to_annotate(retm,d,build=args$build,signif=args$sig)
print(paste0("Number of peaks: ",nrow(peaks)))

if(nrow(peaks)==0 | args$no_annot) {
    peaks=NULL
}else{
#    print("PEAKS")
#    print(peaks)

    # if there are a lot of hits, annotate only the most significant ones
    if(nrow(peaks)>args$maxpeaks) {
      setorder(peaks,p)
#      peaks.col.only=peaks[MAX_NUM_PEAKS:nrow(peaks),]
#      peaks.col.only[,gene:=NA]
#      peaks.col.only[,distance:=NA]
#      peaks.col.only[,consequence:=NA]
#      peaks.col.only[,act:="c"]
      peaks=peaks[1:args$maxpeaks,]
#      peaks[,act:="a"]
    }
#    print("PEAKS")
#    print(peaks)

    context=apply(peaks, 1, function(x){
      u=unlist(get_variant_context(as.numeric(x["chr"]), as.numeric(x["ps"]), x["a1"], x["a2"],build=args$build))
      if(length(u)<3){u[3]="unknown"};
      return(u);
      })

    context=as.data.frame(t(context))
    colnames(context)=c("gene", "distance", "consequence")
    context$distance=as.numeric(as.character(context$distance))
    peaks=cbind(peaks, context)
    if(exists("peaks.col.only")) {
        peaks=rbind(peaks,peaks.col.only)
    }
    peaks$truelabels=as.character(peaks$gene)
    if(args$no_distance){
    	peaks$truelabels[peaks$dist>0]=as.character(peaks$gene[peaks$dist>0])
    }else{
	peaks$truelabels[peaks$dist>0]=paste(as.character(peaks$gene[peaks$dist>0]), paste(" (", ceiling(peaks$distance[peaks$dist>0]/1000), "kb)", sep=""))
    }
    print("FINAL PKS")
    print(peaks)
    peaks$pch=15
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
}

plot_manhattan(retm, peaks, signif=args$sig, MAX_NUM_PEAKS=args$maxpeaks)
dev.off()
}
