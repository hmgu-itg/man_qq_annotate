suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggplot2))

# Get the path to my R functions
funct_path="~ag15/scripts"

# Now make sure that the function files can be found and that they can  be 
# accessed and sourced in
my_libs=c("ensembl_rest_api.R")

# Loop through all the libraries I want to source in
for (i in my_libs ) {
  # Check that the function path exists and is readable
  i_path=paste(funct_path,i,sep="/")
  
  # Make sure the files exist and we can access them
  if (file.access(i_path,4)==-1) {
    stop(paste("[FATAL] Can't find",
     i_path,
     "library (or it is not readable), ensure",
     i_path,
     "is set in your .bashrc !"),sep="")
  }
  
  # Source in the functions required to do the fast manhattan and qq plots
  # I should really drop this into a try - catch block
  source(i_path)
  #  tryCatch(source(i_path),
  #           error = function(e) e,
  #           finally = stop(sprintf("[FATAL] Problem sourcing %s!",i_path)))
}

# So this is accessable to multiple functions
label.cols=c("chr","pos","label_col","poi","peak_col","label_name","pval","xcoord","auto_peak","ycoord")

manual.label.cols=c("chr","pos","label_col","poi","peak_col","label_name")
autopeak.cols=c("chr","pos","pval","xcoord")

###############################################################################
# Calculate a lambda value from a bunch of p-values                           #
# pval = A vector of p-values                                                 #
# round = integer to round the lambda value to                                #
###############################################################################
lambdaCalc = function(pval,round=NULL) {
  print(summary(pval))
print(class(pval))
  lambda = qchisq(median(pval, na.rm=TRUE),1, lower.tail=FALSE)/0.456

  if (is.null(round)==FALSE) {
    lambda=round(lambda,round)
  }
  return(lambda)
}

## Function added by Arthur
isColor <- function(x) {
     sapply(x, function(X) {
         tryCatch(is.matrix(col2rgb(X)), 
                  error = function(e) FALSE)
         })
     }

###############################################################################
# Provides the graphical parameters for the plot. Note that this does not     #
# store any values, it simply presents the values you pass to it in the       #
# context of all the values so that they will all be available to be passed to#
# the functions                                                               #
# ... = graphical parameters, will die if they are incorrect                  #
###############################################################################
mhppar=function(...) {
  defaults=list(mhp.mar=c(5,4,2,2),
    sig.thresh=5E-08,
    sig.thresh.line=5E-08,
    mhp.sig.thresh.cols=c(rgb(0,205,102,maxColorValue=255)),
    sig.thresh.line.plot=FALSE,
    plot.pad=0.05,
    plot.vp.height=0.75,
    lines.vp.height=0.05,
    labels.vp.height=0.2,
    mhp.poi.pch=1,
    mhp.poi.cex=1.5,
    mhp.poi.col="red",
    mhp.pch=20,
    mhp.pch.cex=0.6,
    mhp.xaxis.cex=1.5,
    mhp.yaxis.cex=1.5,
    mhp.chr.cols=c(rgb(0,0,108,maxColorValue=255),rgb(100,149,237,maxColorValue=255)),
    default.label.col="grey50",
    label.cex=0.9,
    label.rot=90,
    label.pad=0.1,
    label.just="left",
    label.fontface="italic",
    label.line.lex=0.7,
    label.line.lty="dashed",
    qq.null.line.col="red",
    qq.null.line.lty="solid",
    qq.null.line.lex=1.5)
  
  passed=list(...)
  
  for (i in names(passed)) {
    if (is.null(defaults[[i]])==TRUE) {
      stop(sprintf("[FATAL] Unrecognised mp argument: %s",i))
    }
    
    if (is.null(passed[[i]])==TRUE) {
      stop(sprintf("[FATAL] no value for: %s",i))
    }
    defaults[[i]]=passed[[i]]
  }
  return(defaults)
}


###############################################################################
# A function to organise the placement of labels above the manhattan plot, it #
# optimises the labels to be as close to their peak as possible.              #
# labels - A data.table with delta_opt                                        #
# nudge.prop - The proportion of the delta_opt that the label will be moved   #
#              in each iteration                                              #
# test_it - The max number of iterations the loop will run for, used for      #
#           debugging                                                         #
###############################################################################
labelNudge = function(labels,nudge.prop=0.1,test_it=-1) {
  # I have not got around to changing the variable names in the function yet 
  variations=labels
  it=1
  # The location of the first/last SNP in the data.table
  snp_bounds=c(1,nrow(variations))
  bound_move_ok=c(-1,1)
  kill=0
  # Basically we will keep on going until nothing can move anymore
  repeat {
    moves=0
    
    loop_order=sort.int(abs(variations[,delta_opt]),index.return=TRUE,decreasing=TRUE)$ix
    
    for (i in loop_order) {
      # Calculate the direction that the label wants to move in, it it is < 0
      # then it wants to move to the left if it is > 0 then it wants to move to
      # the right
      index=sign(variations[i,delta_opt]) * -1
      direction=sign(variations[i,delta_opt])
      
      # If we are optimised, then move to the next SNP
      if (direction==0) {
        next  
      }
      
      # Calculate the distance between the current SNP and 
      # the neighbour it wants to move towards
      dist=variations[(i+bound_move_ok[index]),plot_pos+((-1*bound_move_ok[index])*name_width)]-variations[i,plot_pos]

      # If distance is NA it means that the SNP has no neighbour, so we set
      # the dist to be the delta opt
      if (length(dist)==0 || is.na(dist)==TRUE) {
        dist=variations[i,delta_opt]
      }

      # Now calculate the actual distance that the SNP will move, this will be a 
      # proportion of the delta-optimimum or the distance to the neighbour (whichever)
      # is smaller
      # Our initial move
      move_dist=abs(variations[i,delta_opt]*nudge.prop)

      # We only want to move as far as the nearest neighbour or the
      # delta opt
      target_dist=min(c(abs(dist),abs(variations[i,delta_opt])))

      move_to=target_dist-abs(variations[i,delta_opt]*nudge.prop)
      
      if (move_dist > target_dist || target_dist <= (variations[i,name_width] * 0.5) ) {
        move_dist=target_dist
      }
      
      move_dist=move_dist*direction
      
      # Now move the SNP, re-calcualte the distance from it's optimum
      variations[i,plot_pos:=plot_pos+move_dist]
      variations[i,delta_opt:=opt-plot_pos]

      # update the total moves for this repeat, when this is 0 we will break
      moves=sum(c(moves,abs(move_dist)),na.rm=TRUE)
      
    }

    # Get the indexes of the rows in order of max(abs(delta_opt))->min(abs(delta_opt))
    if (moves == 0 || it == test_it || kill==1) {
      break
    }
    it=it+1
  }
  
  return(variations)  
}


###############################################################################
# Plots any poi on the plot                                                   #
###############################################################################
plot.mhp.poi=function(label.data,mp=mhppar()) {
  # Now move back to the plot region and actually plot some points
  seekViewport("plotRegion")

  # Subset the points of interest and plot them if we still have some
  print(label.data$ycoord)
  print(sapply(label.data$ycoord,FUN=is.numeric))
  x=label.data[poi==TRUE & sapply(ycoord,FUN=is.numeric),xcoord]
  y=label.data[poi==TRUE & sapply(ycoord,FUN=is.numeric),ycoord]
#  print(sapply(label.data[,ycoord],FUN=is.numeric))
  if (length(x)>0) {
    grid.points(
          label.data[poi==TRUE & is.numeric(ycoord),xcoord],
          label.data[poi==TRUE & is.numeric(ycoord),ycoord],
          pch=mp[["mhp.poi.pch"]],
          gp=gpar(col=mp[["mhp.poi.col"]],cex=mp[["mhp.poi.cex"]])
          )
  }

}





plot.mhp.labels=function(label.data,mp=mhppar()) {
  # Navigate the the labels region
  seekViewport("labelsRegion")
  label_viewport=current.viewport()

  # Make sure we only label points that we want to label
  # That is those points where there is some data in the 
  # label_name column
#  print(label.data)
  label.data=label.data[is.null(label_name)==FALSE & 
                        is.na(label_name)==FALSE & 
                        nchar(label_name)>0,]
#  label.data=label.data[nchar(label.data[i,label_name])>0,]
#  print(label.data)
  # Initially we will scale all labels to the same size as the longest label
  # So the labels will be rotated 90 degrees
  
  # So get the longest label
  label.data[,nchar := nchar(label_name)]
  longest_label=label.data[which.max(label.data[,nchar]),label_name]
  
  # Now produce a test grob that will represent the longest label, from this we
  # can calculate the scaling required to fit the label height into the viewport
  testlabel=textGrob(longest_label,
                      rot=mp[['label.rot']],
                      gp=gpar(fontface=mp[['label.fontface']]))

  # Now get the height of the testlabel, we will scale this so it can fit 
  # into an npc of 1
  test_label_height=convertHeight(grobHeight(testlabel),"npc",valueOnly=TRUE)
  test_label_width=convertWidth(grobWidth(testlabel),"npc",valueOnly=TRUE)
   
  cex_text=mp[['label.cex']]/test_label_height

  test_label_height=test_label_height*cex_text
  test_label_width=test_label_width*cex_text

  # Now add the padding to the label width
  test_label_width=test_label_width+(test_label_width*(mp[['label.pad']]*2))
  test_label_width_native=convertWidth(unit(test_label_width,"npc"),"native",valueOnly=TRUE)

  # Now store the width in the label data as a bp unit, this will be used to 
  # calculate the total space taken up by labels
  label.data[,name_width:=convertWidth(unit(test_label_width,"npc"),"native",valueOnly=TRUE)]

   # Get the total space available in the viewport
  total_space=convertWidth(label_viewport$width,"native",valueOnly=TRUE)
  available_space=total_space-(sum(label.data[,name_width]))

  # If we do not have enough available space then we will need to remove some
  # labels from the list and then plot, the labels will be removed in the order of
  # their priority, so the lowest priority labels will be removed until we have 
  # enough labels to fit onto the plot
  if (available_space < 0) {
     no_of_snps_over=ceiling(abs(available_space)/test_label_width_native)
#      print(no_of_snps_over)  
     # Check that we still have some SNPs after pruning, if not then this is a 
     # fatal error, to be honest this should never happen, but just in case!
     if (no_of_snps_over>=nrow(label.data)) {
       stop("[FATAL] I am not sure how this has happened, but you can't fit any snp labels!")
     }

     # Make sure that the label data table is sorted to give manual labels the priority
     # and the lowest pvalues the next priority, then chr,pos
    setkey(label.data,auto_peak,pval)
#    print("****************************************************************8")
#    print(class(label.data[,pval]))
#    print(class(label.data[,mlogp]))
#    print(label.data)
#    print("****************************************************************8")
     # Now prune and re-canculate the available space, we should be able to fit
     # the remainder of the labels onto the plot
     label.data=label.data[c(1:(nrow(label.data)-(no_of_snps_over+1))),]
     available_space=total_space-(sum(label.data[,name_width]))
   }

  #  print(label_viewport$xscale)
  label.data[,xcoord:=as.numeric(xcoord)]
  setkey(label.data,xcoord)

  # Add a SNP order column
  label.data[,label_order:=seq(1,nrow(label.data),1)]

  # So we initialise the SNP positions by dividing the viewport space by the 
  # number of SNPs and plotting the SNPs evenly throughout the available space
  segment_size=total_space/nrow(label.data)
  label.data[,plot_pos:=label_viewport$xscale[1]+(label_order*segment_size)-(segment_size*0.5)]
  label.data[,opt:=(as.numeric(xcoord)+as.numeric(xcoord))/2]
  label.data[,delta_opt:=opt-plot_pos]

  label.data=labelNudge(label.data,nudge.prop=0.1)

  # Now actually plot some snps to the viewport
  for (i in 1:nrow(label.data)) {
    grid.text(label.data[i,label_name],
              x=unit(label.data[i,plot_pos],"native"),
              y=0,
              just=mp[['label.just']],
              rot=mp[['label.rot']],
              gp=gpar(cex=cex_text,
                      col=label.data[i,label_col],
                      fontface=mp[['label.fontface']]))
  }

  # Now add the connecting lines
  seekViewport("linesRegion")
  lines_viewport=current.viewport()
  
  for (i in 1:nrow(label.data)) {
    #    print(i)
    grid.lines(x=unit(c(label.data[i,plot_pos],label.data[i,opt]),"native"),
               y=c(1,0),
               gp=gpar(col=label.data[i,label_col],lty=mp[['label.line.lty']],lex=mp[['label.line.lex']]))
  }
  
  label.data[mlogp==Inf,mlogp:=0]
  seekViewport("plotRegion")
  plot_viewport=current.viewport()

  for (i in 1:nrow(label.data)) {
    grid.lines(x=unit(c(label.data[i,opt],label.data[i,opt]),"native"),
               y=unit(c(label.data[i,mlogp],plot_viewport[['yscale']][2]),"native"),
               gp=gpar(col=label.data[i,label_col],lty=mp[['label.line.lty']],lex=mp[['label.line.lex']])
               )
    
  }
#  print(label.data)
}


###############################################################################
# Process the manual labels that have been supplied by the user               #
###############################################################################
process.manual.labels=function(manual.labels,
                              gene.search.biotype=NULL) {
#  print("IN MANUAL LABELS")
  if (("data.table" %in% class(manual.labels))==FALSE) {
    stop("[FATAL] The manual labels should be a data.table")
  }
#  print(manual.labels)
  # Check that the columns of the manual labels data.table are correct
  # we should have 
  if (identical(names(manual.labels),manual.label.cols)==FALSE) {
    stop(sprintf("[FATAL] Manual labels data.table should have columns %s not: %s",
    paste(manual.label.cols,collapse=","),paste(names(manual.labels),collapse=",")))
  }
    
  # Now loop through the manual labels, if the label field is empty serach 
  # for the nearest gene using chr,pos
  if (nrow(manual.labels)>0) {
    for (i in 1:nrow(manual.labels)) {
      # If there is no data in the manual labels, then go and find the 
      # nearest gene, otherwise, use the labels
#      if (is.null(manual.labels[i,label_name])==TRUE || 
#        is.na(manual.labels[i,label_name])==TRUE || 
#        nchar(manual.labels[i,label_name])==0) {
      if (manual.labels[i,label_name]=="<AUTO>") {
        ng=getNearestGene(manual.labels[i,chr],manual.labels[i,pos],biotypes=gene.search.biotype)
          
        if (nrow(ng)>0) {
          gene_name=ng[1,external_name]
          manual.labels[i,label_name:=gene_name]
        } else {
          manual.labels[i,label_name:=""]
        }
      }
    }
    manual.labels[,pval:=0]
    manual.labels[,xcoord:=0]
    manual.labels[,chr:=as.numeric(chr)]
    manual.labels[,pos:=as.numeric(pos)]
  } else {
    stop("[FATAL] Manual labels is an empty data.table!")
  }

  # Now indicate that these labels are manual labels
  manual.labels[,auto_peak:=FALSE]

  return(manual.labels)
}


###############################################################################
# Process the auto labels that have been defined by the peak finder           #
###############################################################################
process.auto.labels=function(peaks,
                             mp=mhppar(),
                             gene.search.biotype=NULL) {
  # Create an empty autolabels data.table
  auto.labels=data.table(matrix(nrow=1,
                                ncol=length(label.cols),
                                byrow=TRUE,
                                data=NA,
                                dimnames=list(NA,label.cols) ) )
  auto.labels=auto.labels[0]
  
#  print("IN AUTO LABELS")
#  print(peaks)
  if (nrow(peaks)>0) {
    for (i in 1:nrow(peaks)) {
      ng=getNearestGene(peaks[i,chr],peaks[i,pos],biotypes=gene.search.biotype)
    
      if (nrow(ng)>0) {
        chr=peaks[i,chr]
        pos=peaks[i,pos]
        pval=peaks[i,pval]
        gene_name=ng[1,external_name]
        auto.labels=rbindlist(list(auto.labels,
        as.list(c(chr,
                  pos,
                  mp[["default.label.col"]],
                  NA,
                  NA,
                  gene_name,
                  pval,
                  peaks[i, xcoord],
                  TRUE,
                  0))))
      }
    }
  } else {
    stop("[FATAL] (process.auto.labels) No data in peaks")
  }

  return(auto.labels)
}


###############################################################################
# Creates the viewports to hold all the plottin components                    #
# xrange = a vector of 2 with the minimum,maximum range on the x axis         #
# yrange = a vector of 2 with the minimum,maximum range on the y axis         #
# chr.ticks = a vector of xaxis coordinates where chr numbers will be added   #
# mp = the graphical parameters                                               #
###############################################################################
create.mhp.viewports=function(xrange,yrange,chr.ticks,add.label.vp=FALSE,mp=mhppar()) {
  # First create the viewport for the whole region. This uses the plotViewport
  # shortcut function that takes the margins to define a plt region
  pushViewport(plotViewport(mp[["mhp.mar"]],name="wholeRegion"))

  # The default height if we are not adding labels
  plot.vp.height=1
  ext=mp[["plot.pad"]]
  x.ext.range=(xrange[2]-xrange[1])*ext
  y.ext.range=(yrange[2]-yrange[1])*ext
  xscale=c(xrange[1]-x.ext.range,xrange[2]+x.ext.range)
  yscale=c(yrange[1]-y.ext.range,yrange[2]+y.ext.range)

  # If we want to add "space" for gene labels
  if (add.label.vp==TRUE) {
 #   print("IN ADDING LABEL VIEWPORT")
    plot.vp.height=mp[["plot.vp.height"]]
    lines.vp.height=mp[["lines.vp.height"]]
    labels.vp.height=mp[["labels.vp.height"]]

    # Now add the lines track and the labels track
    pushViewport(viewport(xscale=xscale,
                          yscale=yscale,
                          height=unit(lines.vp.height,"npc"),
                          just=c(0,0),
                          x=unit(0,"npc"),
                          y=unit(plot.vp.height,"npc"),
                          name="linesRegion"))
#    grid.rect(gp=gpar(col="red"))
    upViewport()

    pushViewport(viewport(xscale=xscale,
                          yscale=yscale,
                          height=unit(labels.vp.height,"npc"),
                          just=c(0,0),
                          x=unit(0,"npc"),
                          y=unit(plot.vp.height+lines.vp.height,"npc"),
                          name="labelsRegion"))
#    grid.rect(gp=gpar(col="green"))
    upViewport()
  }

  # Add the main plot region viewport
  pushViewport(viewport(xscale=xscale,
                        yscale=yscale,
                        height=unit(plot.vp.height,"npc"),
                        just=c(0,0),
                        x=unit(0,"npc"),
                        y=unit(0,"npc"),
                        name="plotRegion"))
#  grid.rect(gp=gpar(col="blue"))

  # Add a y axis will will add the chromsomes "x" axis manually as it will be 
  # bit more involved
  # I want to make some more fine grained 
  grid.yaxis(gp=gpar(cex=mp[["mhp.yaxis.cex"]]))

  # Now add the chromosome axis, we get the dimentions for this from the whole region
  seekViewport("wholeRegion")
  vp=current.viewport()
  xcoords=vp$x
  ycoords=vp$y - unit(3,"line")
  width=vp$width

  # Seek out the root viewport, it is one above whole region, the chr axis is
  # plotted into the margin space
  upViewport()

  # Create a viewport to hold the chromsome axis
  pushViewport(viewport(xscale=xscale,
                        height=unit(3,"line"),
                        width=width,
                        just=c(0,0),
                        x=xcoords,
                        y=ycoords,
                        name="chrAxis"))

  # Add a chromosome label to the axis  
  grid.text("CHR:",
            x=0,
            y=unit(0.5,"npc"),
            just=c(0.0,0.5),
            gp=gpar(cex=mp[['mhp.xaxis.cex']]))

  # Now add the chromsome numbers
  for (i in 1:length(chr.ticks)) {
    pp= ifelse(i %% 2 == 0, 0.5, 0)
    grid.text(i,
              x=unit(chr.ticks[i],"native"),
              y=unit(pp,"npc"),
              gp=gpar(cex=mp[['mhp.xaxis.cex']]),
              just=c(0.5,0))
  }
  

}


###############################################################################
# data must have coluns p_lrt,chr,ps. if file is NULL then outout will be to  #
# the open graphics device.                                                   #
# mp = graphical parameters as returned by mhppar                             #
# manual.labels = A data.table where each row is a chr,pos to be labelled and #
#                 how they will be labelled                                   #
# auto.labels = Shall we find the peaks and label them                        #
# auto.labels.sigcut = The significance cutoff for labels to process          #
#                      automatically. i.e. don't autolabel any peaks less     #
#                      significant than this                                  #
# gene.search.biotype = vector of gene biotypes to restrict the nearest genes #
###############################################################################
mhp.annotate.peak = function(data,
                            mp=mhppar(),
                            manual.labels=NULL,
                            auto.labels=TRUE,
                            auto.labels.sigcut=5E-8,
                            gene.search.biotype=NULL,
                            padding.factor=0.004) {
#  print(mp)
  # Remove NAs and zero pvals
  data=data[complete.cases(data),]
  data=data[p_lrt!=0,]
  print(summary(data$p_lrt))

  # Make characher chromsomes numeric
  data[chr=="X",chr:="23"]
  data[chr=="Y",chr:="24"]
  data[chr=="MT",chr:="25"]
  data[,chr:=as.numeric(chr)]
  
  # The numbers of data points that will be represented in the x and y plane
  xdp=800
  ydp=600

  # The genome size in mbp, this is used to calcualte the xresolution
#  genome.size=3000000000
#              2950878797
  genome.size=sum(as.numeric(data[,.SD[,max(ps) - min(ps)],.SDcols="ps",by=chr]$V1))
#  bp.padding=12000000
  bp.padding=genome.size*padding.factor

  # The chromosomes that we will loop throug
#  chrs=c(1:23)
  chrs=sort(unique(data[,chr]))

  
  # Make sure we have some chromsomes
  if (length(chrs)==0) {
    stop("[FATAL] There are no chromsomes to plot!")
  }

  # If the user has supplied manual labels
  # These will take priority over any auto.labels
  if (is.null(manual.labels)==FALSE) {
    manual.labels=process.manual.labels(manual.labels,
                                        gene.search.biotype=gene.search.biotype)

    # Now we want to check for duplicated positions in our data.table
    lab.counts=manual.labels[,.N,by=list(chr,pos)]
    lab.counts=lab.counts[N>1,]

    if (nrow(lab.counts)>0) {
      for (i in 1:nrow(lab.counts)) {
        manual.labels[chr==lab.counts[i,chr] & pos==lab.counts[i,pos],pos:=pos+rep(c(0:lab.counts[i,N]-1),length.out=lab.counts[i,N])]
      }
    }

    # This has to be done as R is bloody stupid and if you give col2rgb
    # an NA or NULL it will give you something back, if you give it a ""
    # it will error, so not we can catch any errors!!!
    manual.labels$peak_col=as.character(manual.labels$peak_col)
    manual.labels[sapply(peak_col, FUN=is.na),peak_col:=""]
    manual.labels[sapply(peak_col, FUN=is.null),peak_col:=""]
    manual.labels[is.na(label_col) | is.null(label_col),label_col:="black"]
  }
#  print("manual.labels")
#  print(manual.labels)

  ## Expects data object to be a data.table containing three named columns
  ## chr, ps and p_lrt, representing 
  obspval <- as.numeric(data$p_lrt)
  chr <- as.numeric(data$chr)
  pos <- as.numeric(data$ps)
  
  # This produces an integer -log10 pvalue + 1
  # and is the maximum observed signal in the data set
  obsmax <- trunc(max(-log10(obspval), na.rm=T))+1

  # This returns the index of the positions if the 
  # chr,pos sort and is used to get orderd chr,pos
  # data and pvalue data
  sort.ind <- order(chr, pos) 
  chr <- chr[sort.ind]
  pos <- pos[sort.ind]
  obspval <- obspval[sort.ind]
  
  # Two main vectors for new coordinates, should not be more than picture 
  # resolution. So this vector represents the "resolution" of the plot. This
  # is not the image resolution, rather it is the number of data points that
  # will be represented on the plot. These are pre-created as it will be faster
  # updating the elements of a predefined vector than concatinating
  newx=rep(NA, xdp*ydp)
  newy=rep(NA, xdp*ydp)
  
  # Computes the Xresolution, essentially the number of bp in 1 data-point 
  # space
  # Computes the Y resolution, the number of -log10_pvalues in 1 data-point
  # I am not sure what the *2 does here
  xres=(genome.size/xdp)*2
  yres=(obsmax/ydp)*2
  ## Computes the middles of the Y grid intervals, which will serve as points.
  print(obsmax)
  print(yres)
  breaksy=seq(0, obsmax, by=yres)
  ## Transforming the p-values
  locY=-log10(obspval)
  # Initiallise a vector that will hold the colours for each of the data points
  # that will be plotted
  col=rep(NA, xdp*ydp)
  col1=mp[['mhp.chr.cols']][1]
  col2=mp[['mhp.chr.cols']][2]
  col3=mp[['mhp.sig.thresh.cols']][1]

  # Not sure what this shit is
  coli=rgb(255, 255, 255, maxColorValue=255, alpha=0)
  posi=1
  s=1
  
  # Initializing variables for the main loop, not sure what they do
  # Size store the total base pair range of the data points for each
  # chromsome
  size=0
  size[1:max(unique(chr))]=0
  numpoints=0
  numpoints[1:max(unique(chr))]=0
  labpos=0
  alv=0

  # mi and ma store the minimum and maximum bp positions of the data points
  # on each chromsome
  mi=0
  ma=0

  # t is a list which will be initialised with all the names of the chromsome
  # and the values will be all the base pair positions for each data point
  # on a chromsome
  t=list()

  # Create an empty auto.label.genes data.table, this will be filled with auto
  # peak data or manual peak data depending on the parameters passed
  auto.label.genes=data.table(matrix(nrow=1,
                                ncol=length(label.cols),
                                byrow=TRUE,
                                data=NA,
                                dimnames=list(NA,label.cols) ) )
  auto.label.genes=auto.label.genes[0]
  
  # Create an empty data.table that will store all the automatically annotated
  # peak information, the peaks are "called" during the loop below
  autopeaks=data.table(matrix(nrow=1,
                              ncol=length(autopeak.cols),
                              data=NA,
                              byrow=TRUE,
                              dimnames=list(NA,autopeak.cols)))
  autopeaks=autopeaks[0]

  # Convert the signal threshold to the -log10 scale
  sig.thresh=-log10(mp[['sig.thresh']])

  # Now loop through each chromosome
  for ( i in chrs ) {
    # Get the indexes of the data data.points that are on the current chromsome
    print(i)
    curchr=which(chr == i)

    # If there is no data for the current chromsome, do some stuff that I do
    # not understand and skip the rest of the loop
    if (length(curchr)==0) {
      # Set the range (or size in bp) for the chromsome to 0
      # not sure yet why this is offset by 1
      size[i+1]=0

      # Also set the number of data points for the current chromosome
      # to 0
      numpoints[i+1]=0
      next
    }
    
    # Now assign the colour for the points on the current chromosome
    curcol=ifelse (i%%2==0, col1, col2)

    if(is.null(manual.labels)==FALSE) {
    # make sure that any manual labels that do not have the peak colour 
    # defined are set to be coloured the default chromsome colour
    manual.labels[chr==i & !isColor(peak_col),peak_col:=curcol]
    }
    # We store all the bp positional information for the data points on
    # the current chromsome in a list entry with the name of the current 
    # chromosome
    t[[i]]=pos[curchr]

    # We also store the minimum and maximum base pair positions of the data
    # points for the current chromsome
    mi[i]=min(t[[i]])
    ma[i]=max(t[[i]])

    # Store the bp range for the data points on the current chromsome, I am
    # not 100% sure why this is offset by 1
    size[i+1]=ma[i]-mi[i]

    # Get the number of data points for the current chromosome, I am not sure
    # why this is offset by 1
    numpoints[i+1]=length(t[[i]])

    # Correcting positions: subtracting the start offset and adding length of
    # previous chromosomes. Elements of T should now be continuous (notice that
    # i is a sum.)
    offset=sum(size[1:i]) # The length of all the previous chromsomes

    # So this first normalises all the positions so the first position starts
    # at 0 and then adds the length of all the chromsomes so the positions take
    # on pseudo genomic coordinates (not sure why the current chromosome is 
    # number added)
    t[[i]]=(t[[i]]-mi[i])+offset+i
    
    # The x co-ordinates for the chromsome labels, this is halfway between
    # the current chromosome range. I am not sure why the offset is added
    # as the t[[i]] co-ordinates should already be on the genomic scale?
    # Also, not sure on the (i-1)*bp.padding, this looks like some sort
    # of offset that needs to be bigger as we move up the genome
    labpos[i]=offset+(max(t[[i]])-min(t[[i]]))/2+((i-1)*bp.padding)
    
    # Create x grid for current chromosome, topvalue is the end coordinate 
    # for the current chromsome, again not sure about the +1, The breaks are
    # the position of each element in the grid on the xaxis
    topvalue=(offset+size[i+1]+i)
    breaks=seq((offset+i), topvalue, by=xres)
    
    # seq does not go till the end if by is specified
    if (breaks[length(breaks)] != topvalue) {
      breaks=c(breaks, topvalue)
    }
    
    # Compute histogram of SNPs for the current chromsome according to the
    # grid on the xaxis
    h=hist(t[[i]], breaks=breaks, plot=FALSE)

    # Not sure what this stuff does. I think readmids, is the central 
    # position 
    realmid=h$mids-i-offset+mi[i]
    realpeakval=NULL
    realpeaky=NULL

    # For each interval in this chromosome
    #    -get all corresponding y values
    #    -compute histogram along y grid
    #    -fill non-zero intervals with single middle value
    baseoffset=sum(numpoints[1:i])
    
    # Loop through each of the bins for the current chromsome
    for( j in 1:(length(h$counts)) ){
      # I have no idea about any of this!!!!!!
      print(baseoffset)
      suboffset=sum(h$counts[1:j])-h$counts[j]+baseoffset


      subset=locY[(suboffset+1):(suboffset+h$counts[j])]

      # This is positions of all the SNPs that are in the current "xwindow"
      # or grid column
      subsetx=pos[(suboffset+1):(suboffset+h$counts[j])]
      hy=hist(subset, breaks=breaksy, plot=FALSE)

      # The y coordinates of the points drawn in the current grid
      # window. The x-coordinates are all the same for each grid window
      # so we do not need to define a vector for them
      addendum=hy$mids[hy$counts>0]

      # This is the number of points drawn in the current grid window
      l=length(addendum)

      if (l==0) {
        next
      }

      # newx & newy hold all the points for the whole plot NOT just the 
      # current grid window
      # This defines our xcoords for the current grid window, this is basically
      # repeating the xgrid coordinates for the number of points that we have
      newx[posi:(posi+l-1)]=rep(h$mids[j], l)+((i-1)*bp.padding)
      newy[posi:(posi+l-1)]=addendum

      # The colours of the points for the current grid window, we will modify 
      # the colours for points that the user wants to colour manually
      colvect=rep(curcol, l)

      # If we have some signals that are greater than our significance 
      # threshold then we will add them to the auto peaks so the top
      # signals are available for labelling if we need them
      if (length(addendum[addendum>sig.thresh])>0) {
        # This gets the x&y cordinates for all the points that are above
        # our siginicance threshold so they are available for the auto labeller
        subpos=subsetx[subset>sig.thresh]
        subpval=subset[subset>sig.thresh]
        
        # Modify the colour vector to define different coloured points that are
        # greater than our significace threshold
        colvect[addendum>sig.thresh]=col3 
        # print(autopeaks)
        if(length(subpval)>0){
        appended=cbind(c(i), # chr
                                     c(subpos[subpval==max(subpval)]), # pos
                                     c(10^(-1*(subpval[subpval==max(subpval, na.rm=T)]))), #pval
                                     c(h$mids[j]+((i-1)*bp.padding)) # xcoord
                                     ) 
        if(nrow(appended)>1){appended=appended[1,]}
      autopeaks=rbindlist(list(autopeaks, as.list(appended)) )
      } 
      }

      # If we have supplied some manual labels, then colour everything as the 
      # user wants!
      if(is.null(manual.labels)==FALSE) {
        # Get the extremities of the current grid column
        minx=min(subsetx)
        maxx=max(subsetx)

        # This gets a logical vector for for the manual labels that are in the 
        # current grid
        labels.in.grid=manual.labels[,chr]==i & manual.labels[,pos]>=minx & manual.labels[,pos]<=maxx
        
        # If we have some manual labels that fall in the current grid then get
        # a copy of them
        if(nrow(manual.labels[labels.in.grid,])>0){
          mt=manual.labels[labels.in.grid ,]

          # This tests that the actual base pair psoition of the manual labels in the 
          # are in the data set. if not then we will throw a warning
          ycoords.data=subset[subsetx%in%manual.labels[labels.in.grid, pos]]
          if (nrow(mt) != length(ycoords.data) ) {
            write("[FATAL] One of the following manual.labels is not in the data:",stderr())
            print(mt)
            stop("")
          }
          
          # if we have more than 1 manual label in the current grid then we do
          # not have enough resolution to plot them separately so we generate a 
          # warning and take the first one          
          if (nrow(mt)>1) {
            warning("[WARNING] Peaks too close to label separately! Some information will be lost!")
#            mt=mt[1,]
          }

          manual.labels[labels.in.grid, ycoord:=ycoords.data]

          # ASSERTION: There is always only one manual label per peak. If not and
          # they have different colors, it gonna look like Bazooko Circus
          manual.labels[labels.in.grid, xcoord:=h$mids[j]+((i-1)*bp.padding)]

          # This adds the "downscaled" pval data for any manual labels in
          # the current grid column
          manual.labels[labels.in.grid, pval:=10^(-1*max(addendum))]

          # Make the manual peak the required colour (if one has been specfied)
          #print(mt)
          if (mt[1,peak_col]!=curcol) {
            colvect=rep(mt[1,peak_col], l)
          }
        }
      }
      #print("Final colors transmitted to plot")
      #print(unique(colvect))      
      col[posi:(posi+l-1)]=colvect
      posi=posi+l;
    }
  }


   # if we want to auto label some peaks
  if (auto.labels==TRUE && nrow(autopeaks)>0) {
    auto.label.genes=process.auto.labels(autopeaks,
                                    mp=mp,
                                    gene.search.biotype=gene.search.biotype)
  }

  # If we have submitted some manual labels and have some with data in them
  # then we will join them to the auto labels
  if (is.null(manual.labels)==FALSE && nrow(manual.labels[nchar(label_name)>0,])>0) {
#    print(auto.label.genes)
    auto.label.genes=rbindlist(list(auto.label.genes,manual.labels))
  }

  # At the end of all that label processing then check we have some data, it 
  # is possible that we don't have anything
  if (nrow(auto.label.genes)>0) {
    # Somewhere along the line the pvalues get converted to characters
    # so this makes sure they are numeric
    auto.label.genes[,pval:=as.numeric(pval)]
    auto.label.genes[,xcoord:=as.numeric(xcoord)]
    auto.label.genes[,ycoord:=as.numeric(ycoord)]

    auto.label.genes[,mlogp:=-log10(pval)]
    setkey(auto.label.genes,mlogp)
    auto.label.genes[,order:=seq(nrow(auto.label.genes),1,-1)]
    setkey(auto.label.genes,order)
#    print("GOT HERE3")
  }

  #print("BEFORE UNIQUE")
  #print(auto.label.genes)
  # It is possible that if manual labels and auto labels are both defined that
  # we will get the same peaks replicated twice so, I will exclude peaks with
  # the same chr,pos,label_name to prevent this
  setkey(auto.label.genes,auto_peak,chr,pos)
#  print(auto.label.genes)
  auto.label.genes=unique(auto.label.genes,by=c("chr","label_name"))

#  print("AFTER UNIQUE")
#  print(auto.label.genes)

  newx=na.trim(newx)
  newy=na.trim(newy)
#  out<-file 
#  print(auto.label.genes)
  # If we have labels to add to the plot we need to make space for them
  # in the plotting space
  if ((auto.labels==TRUE || is.null(manual.labels)==FALSE) && nrow(auto.label.genes)>0) {
#    print("IN CREATE MHP LABELS")
    # Create the whole region viewport taking into account the
    # labels
    create.mhp.viewports(c(min(newx),max(newx)),
                         c(min(newy),max(newy)),
                         labpos,
                         add.label.vp=TRUE,
                         mp=mp )

    # Now add the gene labels to the plot, if we have some
    plot.mhp.labels(auto.label.genes,mp=mp)
  } else {
    # If we do not have any labels to add then create the plot
    # region without the label tracks
    create.mhp.viewports(c(min(newx),max(newx)),
                         c(min(newy),max(newy)),
                         labpos,
                         add.label.vp=FALSE,
                         mp=mp )
  }
#  print(auto.label.genes)
  # Now move back to the plot region and actually plot some points
  seekViewport("plotRegion")

  # Now plot the points in the various required colours. I want to plot the 
  # chromosome colours first followed by the other colours that might be 
  # present on the plot
  other.cols=setdiff(na.omit(unique(col)),mp[['mhp.chr.cols']])
  chr.cols=setdiff(na.omit(unique(col)),other.cols)
  #print(auto.label.genes)
  for (i in c(chr.cols,other.cols) ) {
    if (is.na(i)==FALSE) {
      plot_elements=which(col==i)
      grid.points(
        newx[plot_elements],
        newy[plot_elements],
        pch=mp[["mhp.pch"]],
        gp=gpar(col=i,fill=i,cex=mp[["mhp.pch.cex"]])
      )
    }
  }

  # If we want to add a significance threshold line
  if (mp[['sig.thresh.line.plot']]==TRUE) {
#    print("IN HERE!!!!!!!!!!!!!!!!!")
#    print(class(mp[['sig.thresh.line']]))
   grid.lines(y=unit(c(-log10(mp[['sig.thresh.line']]),-log10(mp[['sig.thresh.line']])),"native"),gp=gpar(lty="dashed",col=mp[['default.label.col']]) )
  }

  # Finally we need to add any points that the user wants to add
  # ARTHUR: if there is anything to plot
 if(nrow(auto.label.genes)>0){
  plot.mhp.poi(auto.label.genes,mp=mp)
}
}


###############################################################################
# A function to draw a qqplot                                                 #
# file - A file name to write the plot to or if NULL then it will output to   #
#        the current graphics device                                          #
# lambda - The lambda value to add to the plot or NULL for no lambda value to #
#          be added                                                           #
###############################################################################                          
qqplot = function(data,file=NULL,lambda=NULL,vp=NULL,mp=mhppar()){
  ### QQ plot
  ## Expects data to be a vector of p values
  data=data[complete.cases(data)]
  data=as.numeric(data)
  obspval <- sort(data)
  nrows=length(data)
  logobspval <- -(log10(obspval))
  exppval <- c(1:length(obspval))
  logexppval <- -(log10( (exppval-0.5)/length(exppval)))
  obsmax <- trunc(max(logobspval))+1
  expmax <- trunc(max(logexppval))+1
  
  yres=(max(logobspval)-min(logobspval))/600
  xres=(max(logexppval)-min(logexppval))/600
  ymax=max(logobspval, na.rm=T)
  xmax=max(logexppval, na.rm=T)
  index=1
  indey=1
  newx=rep(NA, 600)
  newy=rep(NA, 600)
  lowx=xmax
  lowy=ymax
  i=1
  # Here logobspval looks like the correct thing...
  print(summary(logobspval))
  while(index<nrows){
    lowx=lowx-xres;
    lowy=lowy-yres;
    before=index;
#    print(lowx)
#print(lowy)
#print(index)
#print("------")
#print (logexppval[index])
#print("------")
#print(logobspval[index])
#print("-----TRUEFALSE-")
#print (logexppval[index]>=lowx & logobspval[index]>=lowy)
#print("------")
#print("NROWS")
#    print(nrows)
#    print(index>nrows)
#print("-----------------------")
    while(logexppval[index]>=lowx & logobspval[index]>=lowy){
#      print(index)
      index=index+1
      while((is.na(logobspval[index]) | is.na(logexppval[index])) & index<=nrows){index=index+1;}
      if (index>nrows){break;}
    }
    lowx=logexppval[index]
    lowy=logobspval[index]
    
    if(before==index){next;}
    newx[i]=logexppval[before]-0.5*xres
    newy[i]=logobspval[before]-0.5*yres
    i=i+1;
  }
  #  library(zoo)
  newx=na.trim(newx)
  newy=na.trim(newy)

  #  newx=na.omit(newx)
  #  newy=na.omit(newy)

  #  out="HP-noNormQQ.png"
  out=file

#  if (is.null(file)==FALSE) {
#    png(out,height=600,width=600)
#  }
  #  print(mp)
  #  op=par(cex.lab=1.40,cex.axis=1.5)
  #  plot(c(0,expmax), c(0,expmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value", xlim=c(0,expmax), ylim=c(0,obsmax), las=1, xaxs="i", yaxs="i", bty="l", main="")
  #  points(x=newx, y=newy, pch=23, cex=.4, bg="black")
  qp=ggplot(data.frame("exp"=newx,"obs"=newy),aes(x=exp, y=obs)) + 
  geom_point() + 
  geom_abline(slope = 1,
              intercept = 0,
              colour = mp[['qq.null.line.col']]) +
  #              size = mp[['qq.null.line.lwd']],
  #              linetype = mp[['qq.null.line.lty']]) +
  xlab(expression(Expected -log[10](italic(p)))) +
  ylab(expression(Observed -log[10](italic(p))))
  if (is.null(lambda)==FALSE) {
  #    print("IN LAMBDA")
  #    print(lambda)
    qp=qp+annotate("text",
                   label=paste(" lambda == ",lambda),
  #                   label="BOLLOCKS",
                   x=0,y=Inf,vjust=1,hjust=0,parse=TRUE)
  #    text(0.1,obsmax - 1,adj = c(0,0),labels=bquote(lambda == .(lambda)),cex=1.7)
  }
  qp=qp+theme(panel.border = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"))
  print(qp)
#  if (is.null(vp)==FALSE) {
  #    print("IN HERE")
  #    print(vp)
#    print(qp,vp=vp)
#  } else {
#    print(qp)
#  }
  # If the user supplies a lambda vale add it to the plot
  #  if (is.null(lambda)==FALSE) {
  #    text(0.1,obsmax - 1,adj = c(0,0),labels=bquote(lambda == .(lambda)),cex=1.7)
  #  }
  #  par(op)
  if (is.null(file)==FALSE) {
    dev.off()
  }
}

