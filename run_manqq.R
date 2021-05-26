#!/usr/bin/env Rscript
# Source in packages and functions
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(manqq))


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

# The input file used to create the plots
parser$add_argument("infile", nargs=1, help="Input file name, must be gzip file")

# The input file used to create the plots
parser$add_argument("outfile", nargs=1, help="Output file name (with no file extension)")

args=parser$parse_args()

manqq::run_manqq(
    infile = args$infile,
    outfile = args$outfile,
    chr = args$chr_col,
    pos = args$pos_col,
    a1 = args$a1,
    a2 = args$a2,
    pval = args$pval_col,
    af = args$af_col,
    maf= args$maf,
    signif = args$sig,
    maxpeaks = args$maxpeaks,
    no_qq = args$no_qq,
    no_man = args$no_man,
    no_annot = args$no_annot,
    no_distance = args$no_distance,
    man_height = args$man_height,
    upper_margin = args$upper_margin,
    annot_cex = args$annot_cex,
    axes_cex = args$axes_cex,
    ylim = args$ylim,
    build = args$build,
    image = args$image
)
