
read.assoc.file = function(filepath, chr, pos, a1, a2, pval, af) {
  data = data.table::fread(filepath, select=c(chr, pos, a1, a2, pval, af))
  data.table::setcolorder(data, c(chr, pos, a1, a2, pval, af))
  data.table::setnames(data, c("chr", "pos", "a1", "a2", "p", "af"))
  return(data)
}

maf.filter = function(data, maf = 0.0) {
  ret=data[af>=maf]
  return(ret)
}

refine_data = function(data) {
  data[,p:=as.numeric(p)]
  data = data[!(is.na(p)) & p!=0]
  return(data)
}


#' @export
manqq_cli = function(infile,
                     outfile,
                     chr,
                     pos,
                     a1,
                     a2,
                     pval,
                     af,
                     maf=0.0,
                     signif=5e-8,
                     maxpeaks=30,
                     no_qq=FALSE,
                     no_man=FALSE,
                     no_annot=FALSE,
                     no_distance=FALSE,
                     man_height=6,
                     upper_margin=2.0,
                     annot_cex=1.1,
                     axes_cex=1.3,
                     ylim=-1,
                     build=38,
                     image='png') {
  # Prepare run config storage
  conf.file = paste0(outfile, '.run_conf')
  cat.conf.file = function(text, append = TRUE) {

    cat(paste0(text, '\n'), file = conf.file, append = append)
  }
  # Store run info in a configuration file
  now = format(Sys.time(), "%Y-%m-%d %H-%M-%S")
  cat.conf.file(paste0('Running ManQQ version: ', packageVersion('manqq')), append = FALSE)
  cat.conf.file(paste0('Timestamp: ', now))
  cat.conf.file(paste0('Input: ', infile))
  cat.conf.file(paste0('Output: ', outfile))
  cat.conf.file(paste0('chr: ', chr))
  cat.conf.file(paste0('pos: ', pos))
  cat.conf.file(paste0('a1: ', a1))
  cat.conf.file(paste0('a2: ', a2))
  cat.conf.file(paste0('pval: ', pval))
  cat.conf.file(paste0('af: ', af))
  cat.conf.file(paste0('maf: ', maf))
  cat.conf.file(paste0('signif: ', signif))
  cat.conf.file(paste0('maxpeaks: ', maxpeaks))
  cat.conf.file(paste0('no_qq: ', no_qq))
  cat.conf.file(paste0('no_man: ', no_man))
  cat.conf.file(paste0('no_annot: ', no_annot))
  cat.conf.file(paste0('no_distance: ', no_distance))
  cat.conf.file(paste0('man_height: ', man_height))
  cat.conf.file(paste0('upper_margin: ', upper_margin))
  cat.conf.file(paste0('annot_cex: ', annot_cex))
  cat.conf.file(paste0('axes_cex: ', axes_cex))
  cat.conf.file(paste0('ylim: ', ylim))
  cat.conf.file(paste0('build: ', build))
  cat.conf.file(paste0('image: ', image))

  # Load input data
  data = read.assoc.file(infile, chr, pos, a1, a2, pval, af)
  # MAF filter data
  data = maf.filter(data, maf)
  # Exclude variants with p value NA or 0
  data = refine_data(data)
  # Make QQ-Plot 
  if (!no_qq) save_qqplot(outfile, data[, p], image)
  # Make Manhattan Plot
  if (!no_man) save_manhattan(data, outfile, height = man_height, signif = signif, maxpeaks = maxpeaks, build = build, image.type = image, no_distance = no_distance, no_annot = no_annot)
}


# #' @export
# run_manqq.gcta = function(infile,
#                           outfile,
#                           maf=0.0,
#                           signif=5e-8,
#                           maxpeaks=30,
#                           no_qq=FALSE,
#                           no_man=FALSE,
#                           no_annot=FALSE,
#                           no_distance=FALSE,
#                           man_height=6,
#                           upper_margin=2.0,
#                           annot_cex=1.1,
#                           axes_cex=1.3,
#                           ylim=-1,
#                           build=38,
#                           image='png') {
#   run_manqq(infile,
#             outfile,
#             'Chr',
#             'bp',
#             'A1',
#             'A2',
#             'p',
#             'Freq',
#             maf,
#             signif,
#             maxpeaks,
#             no_qq,
#             no_man,
#             no_annot,
#             no_distance,
#             man_height,
#             upper_margin,
#             annot_cex,
#             axes_cex,
#             ylim,
#             build,
#             image)
# }

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

