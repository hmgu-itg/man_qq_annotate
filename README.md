# man_qq_annotate
Creates Manhattan and QQ plots with annotated peaks for sequencing-based GWAS outputs, by thinning the dataset to what the eye can see.

## Prerequisites
This package needs R along with the libraries `argparse`, `zoo`, `data.table`.

## Installation
`git clone` the repo, then optionally modify the shebang (1st, starting with `#!`) line of `run_manqq.R` to point it to your local R installation. If you do do this, don't forget to `chmod +x` it.

## Example :
For a GCTA output, use the following:
```
./run_manqq.R --chr-col Chr --pval-col p --pos-col bp --a1 A1 --a2 A2 --build 38 --image png --af-col Freq input.assoc.txt.gz output.prefix
```

The default column names are configured for a GEMMA output file. Input files can be gzipped or plain. Run without arguments for a list of options, run with `--help` for detailed options:

```bash
usage: ./run_manqq.R [-h] [--chr-col [character]] [--pval-col [character]]
                     [--pos-col [character]] [--a1 [character]]
                     [--a2 [character]] [--build [integer]]
                     [--image [character]] [--af-col [character]]
                     [--maf-filter [double]] [--sig [double]]
                     [--maxpeaks [integer]] [--no-qq] [--no-man] [--no-annot]
                     [--no-distance] [--man-height [integer]]
                     [--upper-margin [double]] [--annot-cex [double]]
                     [--axes-cex [double]] [--ylim [double]]
                     infile outfile

A program to plot Manhattan and QQ plots

positional arguments:
  infile                Input file name, must be gzip file
  outfile               Output file name (with no file extension)

optional arguments:
  -h, --help            show this help message and exit
  --chr-col [character]
                        The column NAME for the chromosome column, default chr
  --pval-col [character]
                        The column NAME for the chromosome column, default
                        p_score
  --pos-col [character]
                        The column NAME for the chromosome column, default ps
  --a1 [character]      The column NAME for the effect allele column, default
                        allele1
  --a2 [character]      The column NAME for the non-effect column, default
                        allele0
  --build [integer]     The genome build the positions refer to
  --image [character]   The filetype to save plots to (png or pdf)
  --af-col [character]  The column NAME for the allele frequency column,
                        default af
  --maf-filter [double]
                        The significance threshold for MAF filter, default
                        0.0.
  --sig [double]        The significance threshold to use for peak annotation
  --maxpeaks [integer]  The maximum number of peaks to annotate
  --no-qq               Don't plot QQ.
  --no-man              Don't plot Manhattan.
  --no-annot            Disable peak annotation even if peaks are present.
  --no-distance         Don't add very useful distance to gene info.
  --man-height [integer]
                        Force height of Manhattan in inches. Can have
                        unpredictable consequences (some of which you may
                        regret).
  --upper-margin [double]
                        Y limit of Manhattan plot in units of maximum data
                        points. Even more unpredictable than the above.
  --annot-cex [double]  Size factor for annotations.
  --axes-cex [double]   Size factor for axes and labels.
  --ylim [double]       The y-axis limit (-log10(p))

```
