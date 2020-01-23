#!/software/R-3.3.0/bin/Rscript

library(httr)
library(jsonlite)
library(data.table)

mhp = function(chr, ps, p, X_RES=2000, Y_RES=1000, signif=5e-8) {
  ##Manhat plot
  ## Expects data object to be a list containing three named columns
  ## chr, ps and p_lrt, representing
  obspval <- as.numeric(p)
  chr <- as.numeric(chr)
  pos <- as.numeric(ps)
  print(length(chr))
  print(length(pos))
  obsmax <- trunc(max(-log10(obspval)))+1
  sort.ind <- order(chr, pos)
  chr <- chr[sort.ind]
  pos <- pos[sort.ind]
  obspval <- obspval[sort.ind]

  ## Two main vectors for new coordinates, should not be more than picture resolution
  newx=rep(NA, X_RES*Y_RES)
  newy=rep(NA, X_RES*Y_RES)
  size=0
  mi=0
  ma=0
  t=list()
  posdict=list()


  xres=(3000000000/X_RES)*2
  yres=(obsmax/Y_RES)*2
  breaksy=seq(0, obsmax, by=yres)
  ## Transforming the p-values and defining plot colors
  locY = -log10(obspval)
  col=rep(NA, X_RES*Y_RES)
  col1=rgb(0,0,108,maxColorValue=255)
  col2=rgb(100,149,237,maxColorValue=255)
  col3=rgb(0,205,102,maxColorValue=255)
  coli=rgb(255, 255, 255, maxColorValue=255, alpha=0)
  posi=1
  s=1


  ## Initializing variables for the main loop
  size=0
  numpoints=0
  labpos=0


  for ( i in 1:22 ){
    curchr=which(chr==i)
    curcol=ifelse (i%%2==0, col1, col2)

    t[[i]]=pos[curchr]

    # min and max pos and size of cur chr
    mi[i]=min(t[[i]])
    ma[i]=max(t[[i]])
    size[i+1]=ma[i]-mi[i]
    numpoints[i+1]=length(t[[i]])

    ## Correcting positions: subtracting the start offset and adding length of previous chromosomes.
    ## Elements of T should now be continuous (notice that i is a sum.)
    offset=sum(size[1:i])

    ## the reassignment permanently destroys all positions so we have to save them
    posdict[[i]]=data.frame(pos=t[[i]])
    t[[i]]=(t[[i]]-mi[i])+offset+i
    posdict[[i]]$newpos=t[[i]]

    ## Label positions (for later)
    labpos[i]=offset+(max(t[[i]])-min(t[[i]]))/2+((i-1)*12000000)

    ## Create x grid for current chromosome
    topvalue=(offset+size[i+1]+i)
    breaks=seq((offset+i), topvalue, by=xres)

    ## seq does not go till the end if by is specified
    if(breaks[length(breaks)] != topvalue){breaks=c(breaks, topvalue)}

    ## compute histogram of SNPs according to grid
    h=hist(t[[i]], breaks=breaks, plot=FALSE)

    ## add in the hist coordinates for reference (multiplies exec time by 2 :/ - sadness!)
    posdict[[i]]$poscat=cut(posdict[[i]]$newpos, breaks=breaks, labels=h$mids)
    posdict[[i]]$poscat=as.numeric(as.character(posdict[[i]]$poscat))+((i-1)*12000000)
    posdict[[i]]$chr=i
    ## For each interval in this chromosome:
    ## -get all corresponding y values
    ## -compute histogram along y grid
    ## -fill non-zero intervals with single middle value

    baseoffset= sum(numpoints[1:i])

    for( j in 1:(length(h$counts)) ){
      suboffset= sum(h$counts[1:j])-h$counts[j]+baseoffset
      subset= locY[(suboffset+1):(suboffset+h$counts[j])]


      hy= hist(subset, breaks=breaksy, plot=FALSE)


      addendum= hy$mids[hy$counts>0]
      l= length(addendum)
      if(l==0){next}
      newx[posi:(posi+l-1)]=rep(h$mids[j], l)+((i-1)*12000000)
      newy[posi:(posi+l-1)]=addendum
      colvect=rep(curcol, l)
      colvect[addendum>-log10(signif)]=col3
      col[posi:(posi+l-1)]=colvect


      posi= posi+l;
    }


  }
  posdict=do.call("rbind", posdict)
  u=aggregate(posdict$pos, by=list(posdict$chr, posdict$poscat), FUN=min)
  v=aggregate(posdict$pos, by=list(posdict$chr, posdict$poscat), FUN=max)
  m=merge(u, v, by=c("Group.1", "Group.2"))
  posdict=m
  colnames(posdict)=c("chr", "coord", "min", "max")

  # Remove trailing NAs
  library(zoo)
  newx=na.trim(newx)
  newy=na.trim(newy)
  col=na.trim(col)
return(list(newcoords=data.frame(x=newx, y=newy, col=col), posdict=posdict, labpos=labpos))

}

get_peaks_to_annotate=function (manhattan_object,assoc, signif=5e-8, build=38){
  # expects an object from the mhp function
  ret=NULL
  retm=manhattan_object
  sig=unique(retm$newcoords$x[retm$newcoords$y>-log10(signif)])
#  print("GPTA : ")
#  print(retm$newcoords[retm$newcoords$y>-log10(signif),])
  if(length(sig>1)){
    for(xpos in sig){
      dict_entry=retm$posdict[retm$posdict$coord==xpos,];
#      print(paste("GPTA : for entry", xpos))
#      print(dict_entry)
      mmin=unlist(dict_entry$min)[1];
      mmax=unlist(dict_entry$max)[1];
      chr_ref=unlist(dict_entry$chr)[1];
      #assoc$chr=as.numeric(as.character(assoc$chr))
      peakdata=assoc[(assoc$chr==chr_ref) & (assoc$pos>mmin) & (assoc$pos<mmax),];
#      print("")
#      print(peakdata)
#      print("")
      peakdata=peakdata[peakdata$p==min(peakdata$p),];
#      print("GPTA: extracted minimum")
#      print(peakdata)
      peakdata=peakdata[1,];
      peakdata$plotpos=xpos
      peakdata$ploty=-log10(peakdata$p)
      ret=rbind(ret,peakdata)
    }

  ret=data.table(chr=ret$chr, ps=ret$pos, a1=ret$a1, a2=ret$a2, plotpos=ret$plotpos, ploty=ret$ploty,p=ret$p)
  ret$build=build
  return(ret)
  }else{
    return(data.table(chr=numeric(), ps=numeric(), a1=numeric(), a2=numeric(), plotpos=numeric(), ploty=numeric()))
  }
}


plot_manhattan = function(manhattan_object, annotation_object=NULL, signif=5e-8, MAX_NUM_PEAKS=30){

  yl=ifelse(!is.null(annotation_object), args$upper_margin*max(manhattan_object$newcoords$y), max(manhattan_object$newcoords$y))

  yl=ifelse(args$ylim>-1,args$ylim, yl)
  print(max(manhattan_object$newcoords$y))



  plot(manhattan_object$newcoords$x, manhattan_object$newcoords$y,
    pch=20, col=as.character(manhattan_object$newcoords$col), xlab="",
    axes=F,bty="n", ylim=c(0, yl), yaxt="n", ylab="")

  axis(2, las=1, cex=1.5,at=seq(0, 2*(max(manhattan_object$newcoords$y)%/%2+1), by=2))
 mtext(expression(paste("-log"[10], "(p)")), side=2, line=2.5, at=(max(manhattan_object$newcoords$y)%/%2+1), cex=args$axes_cex)

  if(!is.null(annotation_object)){

    # sh29: split the peaks into the ones to annotate and the ones to only colour in
#    peaks.col.only=annotation_object[act=="c"]
#    annotation_object=annotation_object[act=="a"]

    segments(x0=annotation_object$plotpos,
    y0=annotation_object$ploty,
    y1=1.2*max(manhattan_object$newcoords$y), lty=2, lwd=2, col="lightgray")
    espacement=(max(manhattan_object$newcoords$x)-min(manhattan_object$newcoords$x))
    labelslots=seq(min(manhattan_object$newcoords$x), max(manhattan_object$newcoords$x),by=espacement/MAX_NUM_PEAKS)

      labelpos=apply(annotation_object, 1, function(x) {
        if(length(labelslots)==0) {
          print("Error: too many peaks.")
          return(x["plotpos"])
        }
        slotdist=abs(labelslots-as.numeric(x["plotpos"]))
        idx=(1:length(slotdist))[slotdist==min(slotdist)]
        #print(idx)
        ret=labelslots[idx]
        #print(ret)
        #print(labelslots)
        labelslots<<-labelslots[-idx]
        #print(labelslots)
        return(ret)
        })

      segments(x0=annotation_object$plotpos,x1=labelpos,
        y0=1.2*max(manhattan_object$newcoords$y), y1=1.3*max(manhattan_object$newcoords$y),
        lty=2, lwd=2, col="lightgray")
      text(annotation_object$truelabels, x=labelpos-1e7, y=1.32*max(manhattan_object$newcoords$y),
        srt=45, cex=args$annot_cex, pos=4, font=4)
      points(x=labelpos, y=rep(1.3*max(manhattan_object$newcoords$y), length(labelpos)),
        pch=annotation_object$pch, col=annotation_object$col, font=2, cex=args$annot_cex)
  }
  abline(h=-log10(signif), lwd=2, col="lightgray", lty=3)
    #axis(2,las=1,cex=1.5)
  for (i in 1:22){
    pp=ifelse(i %% 2 == 0, 0, 1)
    mtext(i,1, line=pp,at=manhattan_object$labpos[i],cex=args$axes_cex)
  }
  mtext("Chromosome",1,at=1,cex=args$axes_cex,line=0)

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
              data=fromJSON(query)
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



getVepSnp=function(chr,pos,allele,build=38,
                   name=NULL,
                   query="vep/human/region/%i:%i-%i/%s?content-type=application/json",
                   allow.tries=2) {
  allele=toupper(allele)
  if( build == 38) {
    server="http://rest.ensembl.org"
    } else if( build == 37) {
      server="http://grch37.rest.ensembl.org"
    } else {
      print("ERROR: Invalid build supplied.")
    }

  if (is.null(name)==TRUE) {
    name=sprintf("%i:%i",chr,pos)
  }
  vep_query=sprintf(query,chr,pos,pos,allele)
  r=GET(paste(server, vep_query, sep = "/"), content_type("application/json"))
  vep_data=fromJSON(toJSON(content(r)))

if(!("error" %in% names(vep_data))) {
  return(vep_data)
}

}

get_variant_context=function(chr,pos,a1, a2,build=38) {
print(paste("getting context for ", chr, ":", pos, a1, a2))
   alleles=c(a1, a2)
alleles=toupper(alleles)
    if(build==38) {
      server="http://rest.ensembl.org"
      } else if(build==37) {
        server="http://grch37.rest.ensembl.org"
      } else {
        print("Some warning here")
      }
    ext=paste("/overlap/region/human/", chr, ":", pos, "-", pos, "?feature=gene", sep="")
    r=GET(paste(server, ext, sep = ""), content_type("application/json"))
    stop_for_status(r)
    #return(content(r))
    restr=fromJSON(toJSON(content(r), null="null"))
    ## if the snp overlaps with something else than a prot coding gene we don't care
    if(length(restr)>0 & !("protein_coding" %in% restr$biotype)) {
      restr=NULL
    }
    # if it's intergenic find closest gene
  if(length(restr)==0) {
	st=pos-1e6
	if(st<1){st=1}
      ext=paste("/overlap/region/human/", chr, ":", st, "-", pos+1e6, "?feature=gene", sep="")
      r=GET(paste(server, ext, sep = ""), content_type("application/json"))
      stop_for_status(r)
      restr=fromJSON(toJSON(content(r)))
      if(length(restr)>0 & ("protein_coding" %in% restr$biotype)) {
        restr=restr[restr$biotype=="protein_coding",]
        print("overlap no, distance")
        #print(restr)
        restr$dist1=abs(unlist(restr$start)-pos)
        restr$dist2=abs(unlist(restr$end)-pos)
        #restr$dist=ifelse(restr$dist1<restr$dist2, 1, 0)
        restr$dist=pmin(restr$dist1, restr$dist2)
        #print("====")
        #print(restr)
        #restr$dist[restr$dist==0]=restr$dist1[restr$dist==0]
        #restr$dist[restr$dist==1]=restr$dist2[restr$dist==1]
        #print(restr)
        gene=restr[restr$dist==min(restr$dist),]$external_name[[1]]
        dist=min(restr$dist)
        print(paste("minimum reached for gene", gene, "at", dist))
        #stop()
        # get consequence
        # ext=paste("/overlap/region/human/", chr, ":", pos, "-", pos, "?feature=variation", sep="")
        # r=GET(paste(server, ext, sep = ""), content_type("application/json"))
        # stop_for_status(r)
        # restsnp=fromJSON(toJSON(content(r)))
        cons=NULL
        for(i in alleles) {
          #print(getVepSnp(chr=chr,pos=pos,allele=i,build=build))
          tobind=getVepSnp(chr=chr,pos=pos,allele=i,build=build)
          #cons<<-rbind(cons, tobind, fill=T)
          #print("CONS");print(tobind)
          #print(cons)
          if(!is.null(tobind)){
            if(!is.null(cons)){warning("problem: both alleles have consequences:");print(cons); print(tobind)}
            cons=tobind$most_severe_consequence[1]}
        }
        print(cons)
        return(c(gene,dist,cons))
      }else{
        return(c("none", "0", "intergenic_variant"))
      }
        # if it's inside a gene, do stuff
    } else {
    print("overlap yes")
    ## we are sure from the above test that there is at least 1 prot coding gene
    restr=restr[restr$biotype=="protein_coding",]
    print(restr)
      restr$dist=0
      gene=paste(unlist(restr$external_name), collapse=",")
      cons=NULL
      for(i in alleles) {
          #print(paste("allele", i))
	       #print(colnames(cons))
	       topaste=getVepSnp(chr=chr,pos=pos,allele=i,build=build)
	       #print(colnames(topaste))
	       #if(!("colocated_variants" %in% colnames(topaste))){topaste$colocated_variants=""}

	        #cons=ifelse(is.null(cons), topaste, rbind(cons,topaste))
          if(!is.null(topaste)){
            if(!is.null(cons)){warning("problem: both alleles have consequences");print(cons); print(tobind)}
            cons=topaste$most_severe_consequence[1]
          }
      }
      return(c(gene,0,cons))

    }

}

lambdaCalc = function(pval,round=NULL) {
  lambda = qchisq(median(pval, na.rm=TRUE),1, lower.tail=FALSE)/0.456

  if (is.null(round)==FALSE) {
    lambda=round(lambda,round)
  }
  print(paste("lambda is", lambda))
  return(lambda)
}

## Function added by Arthur
isColor <- function(x) {
     sapply(x, function(X) {
         tryCatch(is.matrix(col2rgb(X)),
                  error = function(e) FALSE)
         })
     }

qqplot = function(data, X_GRID=800, Y_GRID=800){
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
  newx=na.trim(newx)
  newy=na.trim(newy)
  ord=na.trim(ord)
return(data.frame(x=newx, y=newy, order=ord))
}
