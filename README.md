# man_qq_annotate
Creates Manhattan and QQ plots with annotated peaks for sequencing-based GWAS outputs, by thinning the dataset to what the eye can see.

## Example :
For a GCTA output, use the following:
```
./run_manqq.R --chr-col Chr --pval-col p --pos-col bp --a1 A1 --a2 A2 --build 38 --image png --af-col Freq input.assoc.txt.gz output.prefix
```

The default column names are configured for a GEMMA output file. Input files can be gzipped or plain. run without arguments for a list of options, run with `--help` for detailed options.
