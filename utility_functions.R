###############################################################################
# General R utility functions to perform simple "checking" tasks that you     #
# might find common in R scripts                                              #
###############################################################################
suppressPackageStartupMessages(library(data.table))



###############################################################################
# A wrapper around system2 that will fail if the exit code != 0               #
# The command should be a vector with each command component in it            #
###############################################################################
system3 = function(command,progname="") {
  # Now run plink
  ret=system2(command,stdout=TRUE,stderr=TRUE)
  
  # if the length of this is greater than 0 it means
  # an error has been generated so the call failed, if the call
  # could not be made then R will kill the script anyway
  if (length(attr(ret,"status" )) > 0) {
    stop(paste("[FATAL] EXIT CODE",attr(ret,"status" ),progname,"failed!",sep=" "))
  }
}


###############################################################################
# Check that a binary exists and is executable                                #
###############################################################################
binOk = function(binary) {
  #
  binpath=Sys.which(binary)
  
  # Now check that the file path actually exists
  # If not die
  if (file.access(binpath,mode=1)==-1) {
    stop(paste("[FATAL] Can't find ",binary," binary or it is not executable at:",binpath))
  }
  
  return(binpath)
}




###############################################################################
# When given a regexp of column names or column numbers this searches for the #
# column names that match and will return them. If given column numbers 
searchCols = function(col_names,searches,ignore.case=FALSE,dieOnAbsent=TRUE,dieOnRepeated=TRUE) {
print(paste("Searching column ", searches, "in", col_names))
  name_int=FALSE
  # Check to see if the searches are integers
  if (grepl("^[1-9][0-9]*$",searches)==TRUE) {
    name_int=TRUE  
  }
  
  search=grepl(paste("^(",searches,")$",sep=""),col_names,ignore.case=ignore.case)
  sc=table(search)

  if (all(!search) && name_int==FALSE && dieOnAbsent==TRUE) {
    stop(paste("[FATAL] Can't find the required column in ",searches,sep=""))
  } else if ((names(sc) %in% "TRUE")==TRUE && sc["TRUE"]>1 && name_int==FALSE && dieOnRepeated==TRUE) {
    stop(paste("[FATAL] Found multiple columns in",searches,sep=""))
  }
  
  
  if (name_int==TRUE) {
    search=as.integer(searches)
    
    if ( search > length(col_names) ) {
      stop(paste("[FATAL] Column number",search," out of range ",length(col_names)))
    }
  }
  
  return(col_names[search])
}


###############################################################################
# This checks if all the required columns are present in a list of column     #
# names. If dieOnError=TRUE, then it will stop execution if some names are    #
# absent                                                                      #
###############################################################################
checkCols=function(col_names,req_cols,dieOnError=FALSE) {
  names_absent=which(req_cols %in% col_names==FALSE)

  if (length(names_absent)>0) {
    miss_col=paste(req_cols[names_absent],collapse=",")
    
    
    if (dieOnError==TRUE) {
      stop(paste("[FATAL] (checkCols) Can't find the following columns:",miss_col))
    }
    
    warning(paste("[WARN] (checkCols) Can't find the following columns:",miss_col))
  }
  return(names_absent)
}


###############################################################################
# Given a list of names and a data.table. This will make sure that the data in#
# those column names are the type requested                                   #
###############################################################################
castCols=function(data,colNames,type="integer") {
  colNames=colNames[!(colNames %in% colNames[checkCols(names(data),colNames)])]
#  print(colNames)
  if (type=="integer") {
    for (n in colNames) {
#      print(n)
#      print(data[,n,with=FALSE])
      data[,n := as.numeric(unlist(data[,n,with=FALSE])),with=FALSE]
    }
  } else {
    stop(paste("[FATAL] Unknown type:",type,"(possibly I have not implemented this yet!)"))
  }
  
  return(data)
}


###############################################################################
# This checks to see if the given file extension is on the end of the file    #
# name. The extension match is a whole match i.e. .txt. The . doesn't have to #
# be supplied in the extension. If it isn't, it will be added, if it is then  #
# the existing filname will just be returned.                                 #
###############################################################################
checkExt=function(file,ext) {
  # Add a dot onto the end of the file (assuming
  # no dot exists)
  ext=sprintf(".%s",sub("^\\.","",ext,perl=TRUE))
  
  # Create a regex for the file extension, so I can
  # test that if it is present
  ext_reg=paste('\\',ext,"$",sep="")
  
  # If there is no match for the file extension, then
  # add the file extension and return
  if (regexec(ext_reg,file)==-1) {
    return(sprintf("%s%s",file,ext))
  }
  
  # If we get here then it has matched and so we just
  # return the file name
  return(file)
}


###############################################################################
# This splits the string on the delimiter, it also gives the option to remove #
# spaces from the string before splitting                                     #
###############################################################################
splitDelimiter=function(string,delimiter=",",remove_space=TRUE) {
  if (remove_space==TRUE) {
    return( unlist(strsplit(gsub("\\s+","",string,perl=TRUE),delimiter) ) )
  } else {
    return( unlist(strsplit(string,delimiter) ) )
  }
}


###############################################################################
# When given a file and some column names, this will use fread to load only   #
# those columns from the file. It is a fatal error if any of the columns are  #
# absent. The column names can be regexps. It is also a fatal error if >1     #
# column name is identified                                                   #
###############################################################################
smartFread=function(file,colnames=NULL,map=NULL,key=NULL,ignore.case=TRUE) {
#  print(map)
#  print(colnames)
#  print("======")
  # If we have not supplied column names
  if (is.null(colnames)==TRUE) {
    stop("[FATAL] You must supply column names!")
  }
  
  # Make sure the file is readable
  if (file.access(file,4)==-1) {
    stop(paste("[FATAL] The file does not exist or you don't have read permissions:", file))
  }
  
  # If we are providing a map it has to be the same length as 
  if (is.null(map)==FALSE && length(map) != length(colnames)) {
    stop("[FATAL] The map is not the same length as the colnames")
  }
  
  # This step does a quick read of the first line of the data.table
  # so we can get the header information without loading the whole file
#  pd=data.table()
  
  # If the file is gzipped try to use zcat to read it in
  if (length(grep(".gz$",file))==1 || length(grep(".bgz$",file))==1) {
    header=names(fread(paste("zcat",file),nrows=0))
  } else {
    header=names(fread(file,nrows=0))
  }
  
  # This will hold a vector of column names we can't to extract
  colnos=c()
  
  # This will hold the actual column names we will extract, we need
  # these later as we may have supplied regexps to serach for column 
  # names so we want to no which ones we have. We also want to make sure
  # that we have not pulled the same column name twice, this may happen if
  # the regular expression is too loose
  found_columns=c()
  
  # Loop through all the column names we require
  for (cn in colnames) {
    # search for a match with the column name in the file header
    # If we find more than one match it should die, also if we find
    # nothing it should die
    col=searchCols(header,
                   cn,
                   ignore.case=ignore.case,
                   dieOnAbsent=TRUE,
                   dieOnRepeated=TRUE)

    # Find the column number for the matching column in the header
    colnos=c(colnos,which(header==col))
    
    # Store the column name that we have found for later
    found_columns=c(found_columns,col)
  }
  
  # now check that we have not got repeats of any columns, in theory if any 
  # single regexp column name supplied give duplicated columns then it will be 
  # caught by searchCols. However, if multiple regexp column names return the
  # same column then we would not know about it so that is why we are checking 
  # it here
  counts=table(found_columns)
#  print(names(counts))
  repeated_cols=names(counts)[counts > 1]
  if (length(repeated_cols>0)) {
    stop(sprintf("[FATAL] The following columns have been requested twice (regexps too loose?): %s",paste(repeated_cols,collapse=",")) )
  }
  
#  print(found_columns)
#  print(colnos)
  # Now we get the data.table selecting only the column names we require
  pd=data.table()
  if (length(grep(".gz$",file))==1 || length(grep(".bgz$",file))==1) {
    pd=fread(paste("zcat",file),select=colnos)
  } else {
    pd=fread(file,select=colnos)
  }
#  print(names(pd))
  # If we have supplied a additional set of column names we want to be on
  # our data.table then we need to work out the order they occur. This may
  # not be nessesary but I will do it just in case.
  if (is.null(map)==FALSE) {
    # Now map the map to the existing colnames so I can rename
    map_col_names=c()
    
    # Loop through all the existing column names in the data table
    for (i in names(pd)) {
      map_col_names=c( map_col_names,map[which(found_columns==i)] )
    }
#    print(map_col_names)
#    print(map_col_names)
    setnames(pd,names(pd),map_col_names)
  }
  
  # If we have supplied a key column then set it before returning
  if (is.null(key)==FALSE) {
    setkeyv(pd,key)
  }
  
  return( pd )
}


## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}



