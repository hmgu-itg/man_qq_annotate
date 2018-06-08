# man_qq_annotate
Creates Manhattan and QQ plots with annotated peaks.

## Example :

```
/software/R-3.4.0/bin/Rscript ~sh29/repos/man_qq/run_manqq.R \
                --chr-col chr \
                --pos-col ps \
                --pval-col p_score \
                --image png \
                ${dir}/META/META.ROR1/MANOLIS.META.ROR1.maf0.001.assoc.txt.gz \
                MANOLIS.ROR1.maf0.001.assoc
```
