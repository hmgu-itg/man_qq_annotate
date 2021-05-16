
read.assoc.file = function(filepath, chr, pos, a1, a2, pval, af) {
  data = data.table::fread(filepath, select=c(chr, pos, a1, a2, pval, af))
  data.table::setcolorder(data, c(chr, pos, a1, a2, pval, af))
  data.table::setnames(data, c("chr", "pos", "a1", "a2", "p", "af"))
  return(data)
}


maf.filter = function(data, maf = 0.0) {
  return(data[af>=maf])
}

#' @export
run_manqq = function(infile, outfile, chr, pos, a1, a2, pval, af, maf=0.0, image.type='png') {
  data = read.assoc.file(infile, chr, pos, a1, a2, pval, af)
  data = maf.filter(data, maf)

  data = data[!(is.na(p)) & p!=0] # Exclude variants with p value NA or 0

  manqq.qqplot(outfile, data[, p], image.type)
  manqq.manhattan(data, outfile, height = 6, signif = 5e-8, maxpeaks = 30, build = 38, image.type = 'png', no_distance = FALSE, no_annot = FALSE)
}


###############################################################################
# A generalised function for running a query against the ensembl rest API,this#
# just has a bunch of error handling stuff in it                              #
###############################################################################
runEnsemblQuery=function(query,allow.tries=2) {
  data=NULL
  tries=0
  retry=TRUE
  while(tries != allow.tries && retry==TRUE) {
    retry=FALSE
    data=tryCatch({
              data=jsonlite::fromJSON(query)
              },
             warning=function(war) {
                write(sprintf("[WARN] The query %s generated the following warning: %s",query,war$message),stdout())
                return(NULL)
             },
             error=function(err) {
               write(sprintf("[ERROR] The query %s generated the following warning: %s",query,err$message),stdout())
               return(NULL)
             }
    )

    if (is.null(data)==TRUE) {
      retry=TRUE
    }
    tries=tries+1
  }
  return(data)
}




## Function added by Arthur
isColor <- function(x) {
     sapply(x, function(X) {
         tryCatch(is.matrix(col2rgb(X)),
                  error = function(e) FALSE)
         })
     }

