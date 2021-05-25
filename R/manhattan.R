

manqq.manhattan = function(data, outfile, height=6, signif = 5e-8, maxpeaks = 30, build = 38, image.type = 'png', no_distance = FALSE, no_annot = FALSE) {

  mhp.file = paste0(outfile, ".man.", image.type)
  if(image.type=="pdf") {
    pdf(mhp.file, width=10, height = height)
  } else if(image.type=="png") {
    png(mhp.file, width=10, height = height, units="in", res=300)
  }

  retm = mhp(data[,chr], data[,pos], data[,p], signif=signif)
  print("Finding peaks...")
  peaks = get_peaks_to_annotate(retm, data, build=build, signif=signif)
  print(paste0("Number of peaks: ", nrow(peaks)))

  if(nrow(peaks)==0 | no_annot) {
    peaks=NULL
  } else {
    # if there are a lot of hits, annotate only the most significant ones
    if(nrow(peaks)>maxpeaks) {
      data.table::setorder(peaks,p)
      peaks=peaks[1:maxpeaks,]
    }
    # A matrix with
    # 1st row: Overlapping gene names?
    # 2nd row: Distance
    # 3rd row: variant type (e.g. intron_variant)
    # Row-wise apply
    context = apply(peaks, 1, function(x){
      u = unlist(get_variant_context(as.numeric(x["chr"]), as.numeric(x["ps"]), x["a1"], x["a2"], build=build))
      if (length(u)<3){u[3]="unknown"}
      return(u)
    })

    context=as.data.frame(t(context))
    colnames(context)=c("gene", "distance", "consequence")
    context$distance=as.numeric(as.character(context$distance))
    peaks=cbind(peaks, context)
    if(exists("peaks.col.only")) {
      peaks=rbind(peaks,peaks.col.only)
    }
    peaks$truelabels=as.character(peaks$gene)
    if(no_distance){
      peaks$truelabels[peaks$dist>0]=as.character(peaks$gene[peaks$dist>0])
    }else{
      peaks$truelabels[peaks$dist>0]=paste(as.character(peaks$gene[peaks$dist>0]), paste(" (", ceiling(peaks$distance[peaks$dist>0]/1000), "kb)", sep=""))
    }
    print("FINAL PKS")
    print(peaks)
    peaks$pch=15
    peaks$col="forestgreen"
    lof = c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant")
    high = c("stop_lost", "start_lost", "transcript_amplification")
    exonic = c("inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant")
    low = c("splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant")
    intronic = c("intron_variant")
    intergenic = c("intergenic_variant")
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

  plot_manhattan(retm, peaks, signif=signif, MAX_NUM_PEAKS=maxpeaks)
  dev.off()
}



#' Make manhattan object
#'
#' @param chr Chromosome column of data.table
#' @param ps Position column of data.table
#' @param p P-value column of data.table
#' @param X_RES
#' @param Y_RES
#' @param signif Significance threshold (default: 5e-8)
#' @return A list with 'newcoords', 'posdict', and 'labpos' attributes
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

  newx=zoo::na.trim(newx)
  newy=zoo::na.trim(newy)
  col=zoo::na.trim(col)
  return(list(newcoords=data.frame(x=newx, y=newy, col=col), posdict=posdict, labpos=labpos))
}


#' Query Ensembl for overlapping genes in a given genomic region
#'
#' @param chr Chromosome
#' @param start Start genomic position
#' @param end End genomic position
#' @param build The assembly to query
#' @return data.frame with query results. Empty data.frame returned if no overlap
#' @examples
#' query_ensembl_gene_overlap(9, 5073770, 5073770)
#' query_ensembl_gene_overlap(9, 4973770, 5173770)
query_ensembl_gene_overlap = function(chr, start, end, build=38) {
  if(build==38) {
    server="http://rest.ensembl.org"
  } else if(build==37) {
    server="http://grch37.rest.ensembl.org"
  }
  ext = paste0("/overlap/region/human/", chr, ":", start, "-", end, "?feature=gene")
  r = httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
  httr::stop_for_status(r)
  restr = jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))

  if (!is.data.frame(restr) && is.list(restr) && length(restr)==0) {
    return(data.frame())
  }
  return(restr)
}

process_overlap_restr = function(restr) {
  prot_genes = restr[restr$biotype=="protein_coding",]
  if (is.null(prot_genes$external_name)){
    gene_names = prot_genes$id
  } else {
    gene_names = prot_genes$external_name
  }
  return (paste(gene_names, collapse = ','))
}

process_case2_restr = function(restr, pos) {
  restr = restr[restr$biotype=="protein_coding",]
  print("overlap no, distance")

  restr$dist1 = abs(unlist(restr$start)-pos)
  restr$dist2 = abs(unlist(restr$end)-pos)

  restr$dist = pmin(restr$dist1, restr$dist2)

  reordered_restr = restr[order(restr$dist), ]

  if (is.null(unlist(reordered_restr$external_name))) {
    # None of the protein coding genes have an external_name (they are all NULL).
    # Get the closest gene's Ensembl Stable ID instead.
    # (This case is probably because the region only has novel protein coding genes)
    gene = reordered_restr$gene_id[[1]]
    dist = reordered_restr$dist[[1]]
    return(list(gene = gene, dist = dist))
  }
  print(reordered_restr)
  for (i in 1:nrow(reordered_restr)) {
    gene = unlist(reordered_restr[i, 'external_name'])
    dist = reordered_restr[i, 'dist']
    if (is.null(gene)) {
      print('[INFO] Skipping gene with no external_name')
      print(reordered_restr[i,])
      next
    }
    break
  }
  print(paste("minimum reached for gene", gene, "at", dist))
  assertthat::assert_that(!is.null(gene), !is.null(dist))
  print(dist)
  return(list(gene = gene, dist = dist))
}

get_variant_context = function(chr, pos, a1, a2, build=38) {
  print(paste("getting context for ", chr, ":", pos, a1, a2))
  alleles = c(a1, a2)
  alleles = toupper(alleles)

  restr = query_ensembl_gene_overlap(chr, pos, pos, build)

  # Three possible cases:
  # 1. Variant directly overlaps region of protein coding gene(s)
  # 2. Variant is intergenic, but protein coding gene is within 2Mb window (1Mb either side)
  # 3. Variant is intergenic, and no protein coding gene present within 2Mb window (1Mb either side)
  if ("protein_coding" %in% restr$biotype) {
    # Case 1
    print("overlap yes")
    ## we are sure from the above test that there is at least 1 prot coding gene
    gene = process_overlap_restr(restr)
    dist = 0
    assertthat::assert_that(!is.null(gene))

  } else {
    # Case 2 and 3 start
    st = pos-1e6
    st = ifelse(st<1, 1, st)
    restr = query_ensembl_gene_overlap(chr, st, pos+1e6, build)
    if (nrow(restr)==0 | !("protein_coding" %in% restr$biotype)) {
      # Case 3
      # Return and finish function call as we don't need to query consequence
      return(c("none", 0, "intergenic_variant"))
    } else {
      # Case 2
      output = process_case2_restr(restr, pos)
      gene = output$gene
      dist = output$dist
    }

  }


  # get consequence
  cons = NULL
  for(i in alleles) {
    topaste=getVepSnp(chr=chr,pos=pos,allele=i,build=build)
    if(!is.null(topaste)){
      if(!is.null(cons)){
        # TODO: Find an example and ask Arthur how to handle this case.
        warning("problem: both alleles have consequences")
        print(cons)
        print(topaste)
      }
      cons = topaste$most_severe_consequence[1]
    }
  }
  assertthat::assert_that(!is.null(gene), !is.null(dist), !is.null(cons))
  return(c(gene, dist, cons))
}


#' Get peaks to annotate
#'
#' @param manhattan_object Output object of `mhp` function
#' @param assoc data.table with columns chr, pos, a1, a2, p, and af
#' @param signif Significance threshold (default: 5e-8)
#' @param build Genome assembly number (38 or 37)
#' @return A data.table with columns: plotpos, ploty, p, build
get_peaks_to_annotate = function (manhattan_object, assoc, signif = 5e-8, build = 38){
  # expects an object from the mhp function
  ret=NULL
  retm=manhattan_object
  sig=unique(retm$newcoords$x[retm$newcoords$y>-log10(signif)])
  # print("GPTA : ")
  # print(retm$newcoords[retm$newcoords$y>-log10(signif),])
  if(length(sig>1)){
    for(xpos in sig){
      dict_entry=retm$posdict[retm$posdict$coord==xpos,];
      # print(paste("GPTA : for entry", xpos))
      # print(dict_entry)
      mmin=unlist(dict_entry$min)[1];
      mmax=unlist(dict_entry$max)[1];
      chr_ref=unlist(dict_entry$chr)[1];
      # assoc$chr=as.numeric(as.character(assoc$chr))
      peakdata=assoc[(assoc$chr==chr_ref) & (assoc$pos>mmin) & (assoc$pos<mmax),];
      # print("")
      # print(peakdata)
      # print("")
      peakdata=peakdata[peakdata$p==min(peakdata$p),];
      # print("GPTA: extracted minimum")
      # print(peakdata)
      peakdata=peakdata[1,];
      peakdata$plotpos=xpos
      peakdata$ploty=-log10(peakdata$p)
      ret=rbind(ret,peakdata)
    }

    ret = data.table::data.table(chr=ret$chr, ps=ret$pos, a1=ret$a1, a2=ret$a2, plotpos=ret$plotpos, ploty=ret$ploty,p=ret$p)
    ret$build = build
    return(ret)
  } else {
    return(data.table::data.table(chr=numeric(), ps=numeric(), a1=numeric(), a2=numeric(), plotpos=numeric(), ploty=numeric()))
  }
}

plot_manhattan = function(manhattan_object,
                          annotation_object = NULL,
                          signif = 5e-8,
                          MAX_NUM_PEAKS = 30,
                          upper_margin = 2.0,
                          ylim = -1.0,
                          annot_cex = 1.1,
                          axes_cex = 1.3) {
  yl=ifelse(!is.null(annotation_object), upper_margin*max(manhattan_object$newcoords$y), max(manhattan_object$newcoords$y))
  yl=ifelse(ylim>-1, ylim, yl)
  print(max(manhattan_object$newcoords$y))

  plot(manhattan_object$newcoords$x, manhattan_object$newcoords$y,
    pch=20, col=as.character(manhattan_object$newcoords$col), xlab="",
    axes=F,bty="n", ylim=c(0, yl), yaxt="n", ylab="")

  axis(2, las=1, cex=1.5,at=seq(0, 2*(max(manhattan_object$newcoords$y)%/%2+1), by=2))
  mtext(expression(paste("-log"[10], "(p)")), side=2, line=2.5, at=(max(manhattan_object$newcoords$y)%/%2+1), cex=axes_cex)

  if(!is.null(annotation_object)){

    # sh29: split the peaks into the ones to annotate and the ones to only colour in
    # peaks.col.only=annotation_object[act=="c"]
    # annotation_object=annotation_object[act=="a"]

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
        # print(idx)
        ret=labelslots[idx]
        # print(ret)
        # print(labelslots)
        labelslots<<-labelslots[-idx]
        # print(labelslots)
        return(ret)
        })

      segments(x0=annotation_object$plotpos,x1=labelpos,
        y0=1.2*max(manhattan_object$newcoords$y), y1=1.3*max(manhattan_object$newcoords$y),
        lty=2, lwd=2, col="lightgray")
      text(annotation_object$truelabels, x=labelpos-1e7, y=1.32*max(manhattan_object$newcoords$y),
        srt=45, cex=annot_cex, pos=4, font=4)
      points(x=labelpos, y=rep(1.3*max(manhattan_object$newcoords$y), length(labelpos)),
        pch=annotation_object$pch, col=annotation_object$col, font=2, cex=annot_cex)
  }
  abline(h=-log10(signif), lwd=2, col="lightgray", lty=3)
  # axis(2,las=1,cex=1.5)
  for (i in 1:22){
    pp=ifelse(i %% 2 == 0, 0, 1)
    mtext(i,1, line=pp,at=manhattan_object$labpos[i],cex=axes_cex)
  }
  mtext("Chromosome",1,at=1,cex=axes_cex,line=0)

}
