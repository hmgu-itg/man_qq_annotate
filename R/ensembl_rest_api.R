# All the rest API calls are done via rjson
# TOOD: class/type checking within functions


# These are for "Pretty" Population Names rather than the horrible ones that are returned by ensembl
longpop=c(
  "1000GENOMES:phase_1_AFR",
  "1000GENOMES:phase_1_ALL",
  "1000GENOMES:phase_1_AMR",
  "1000GENOMES:phase_1_ASN",
  "1000GENOMES:phase_1_ASW",
  "1000GENOMES:phase_1_CEU",
  "1000GENOMES:phase_1_CHB",
  "1000GENOMES:phase_1_CHS",
  "1000GENOMES:phase_1_CLM",
  "1000GENOMES:phase_1_EUR",
  "1000GENOMES:phase_1_FIN",
  "1000GENOMES:phase_1_GBR",
  "1000GENOMES:phase_1_IBS",
  "1000GENOMES:phase_1_JPT",
  "1000GENOMES:phase_1_LWK",
  "1000GENOMES:phase_1_MXL",
  "1000GENOMES:phase_1_PUR",
  "1000GENOMES:phase_1_TSI",
  "1000GENOMES:phase_1_YRI",
  "BUSHMAN:BANTU",
  "CLINSEQ_SNP:CSAgilent",
  "COMPLETE_GENOMICS:YRI",
  "ENSEMBL:ENSEMBL_Watson",
  "ILLUMINA-UK:YRI",
  "SEQUENOM:CEPH",
  "1000GENOMES:pilot_1_CEU_low_coverage_panel",
  "1000GENOMES:pilot_1_CHB+JPT_low_coverage_panel",
  "1000GENOMES:pilot_1_YRI_low_coverage_panel",
  "CGAP-GAI:POOLED_CEPH",
  "COMPLETE_GENOMICS:PGP",
  "CORNELL:AGI_ASP_population",
  "CSHL-HAPMAP:HAPMAP-ASW",
  "CSHL-HAPMAP:HAPMAP-CHB",
  "CSHL-HAPMAP:HAPMAP-CHD",
  "CSHL-HAPMAP:HAPMAP-GIH",
  "CSHL-HAPMAP:HAPMAP-LWK",
  "CSHL-HAPMAP:HAPMAP-MEX",
  "CSHL-HAPMAP:HAPMAP-MKK",
  "CSHL-HAPMAP:HAPMAP-TSI",
  "CSHL-HAPMAP:HapMap-CEU",
  "CSHL-HAPMAP:HapMap-HCB",
  "CSHL-HAPMAP:HapMap-JPT",
  "CSHL-HAPMAP:HapMap-YRI",
  "ENSEMBL:ENSEMBL_Venter",
  "ENSEMBL:ENSEMBL_celera",
  "ESP6500:African_American",
  "ESP6500:European_American",
  "HUMANGENOME_JCVI:J. Craig Venter",
  "KYUGEN:CHMJ",
  "NCBI:NIHPDR",
  "PERLEGEN:AFD_AFR_PANEL",
  "PERLEGEN:AFD_CHN_PANEL",
  "PERLEGEN:AFD_EUR_PANEL",
  "WIAF-CSNP:WIAF-CSNP-MITOGPOP5",
  "NHLBI-ESP:ESP_Cohort_Populations",
  "BUSHMAN:BUSHMAN_POP2"
)

shortpop=c(
  "1000g:ph1:AFR",
  "1000g:ph1:ALL",
  "1000g:ph1:AMR",
  "1000g:ph1:ASN",
  "1000g:ph1:ASW",
  "1000g:ph1:CEU",
  "1000g:ph1:CHB",
  "1000g:ph1:CHS",
  "1000g:ph1:CLM",
  "1000g:ph1:EUR",
  "1000g:ph1:FIN",
  "1000g:ph1:GBR",
  "1000g:ph1:IBS",
  "1000g:ph1:JPT",
  "1000g:ph1:LWK",
  "1000g:ph1:MXL",
  "1000g:ph1:PUR",
  "1000g:ph1:TSI",
  "1000g:ph1:YRI",
  "BUSHMAN:BANTU",
  "CLINSEQ:CSAgi",
  "CG:YRI",
  "ENSEMBL:Watson",
  "ILLU-UK:YRI",
  "SEQUENOM:CEPH",
  "1000g:pi1:CEU",
  "1000g:pi1:CHB+JPT",
  "1000g:pi1:YRI",
  "CGAP-GAI:CEPH",
  "CG:PGP",
  "CORNELL:AGI+ASP",
  "HAPMAP:HAPMAP-ASW",
  "HAPMAP:HAPMAP-CHB",
  "HAPMAP:HAPMAP-CHD",
  "HAPMAP:HAPMAP-GIH",
  "HAPMAP:HAPMAP-LWK",
  "HAPMAP:HAPMAP-MEX",
  "HAPMAP:HAPMAP-MKK",
  "HAPMAP:HAPMAP-TSI",
  "HAPMAP:HapMap-CEU",
  "HAPMAP:HapMap-HCB",
  "HAPMAP:HapMap-JPT",
  "HAPMAP:HapMap-YRI",
  "ENSEMBL:Venter",
  "ENSEMBL:celera",
  "ESP6500:Afr_Am",
  "ESP6500:Eur_Am",
  "HG_JCVI:Venter",
  "KYUGEN:CHMJ",
  "NCBI:NIHPDR",
  "PERLEGEN:AFD+AFR",
  "PERLEGEN:AFD+CHN",
  "PERLEGEN:AFD+EUR",
  "WIAF-CSNP:MITOGPOP5",
  "NHLBI-ESP:ESP_POP",
  "BUSHMAN:POP2"
)
# poptranslate=data.table(pop=longpop,shortpop=shortpop)

###############################################################################
# return.maf = Return minor alleles and frequencies rather than everything    #
###############################################################################
getPopAF=function(s,
                  build="hg38",
                  query="variation/human/%s?content-type=application/json;pops=1",
                  pop=NULL
                  ) {

  if( build == "hg38") {
  server="http://rest.ensembl.org"
  } else if( build == "hg19") {
    server="http://grch37.rest.ensembl.org"
  }

  write(paste("[INFO] Processing ",s,sep=""),stdout() )
  popdata=data.table()
  #  print(grep("^(chr)?[1-9XY]?[0-9]+:\\d+", s, perl = TRUE))
  # Check the SNP is not chr:pos format as we can't handle them yet!
  if (length(grep("^(chr)?[1-9XY]?[0-9]+:\\d+",s,perl=TRUE))>0 ) {
    stop("[FATAL] At the moment this can only handle rs numbers not chr:pos")
  }

  # Build the JSON string. If the SNP does not exist then it will
  # crash and die. I need to make it exit gracefully at some point
  variation_query=sprintf(query,s)
  variation_query=paste(server,variation_query,sep="/")
  snpdata=fromJSON(file=variation_query)

  # Check that there is population data for the SNP, if not return NULL
  # and generate a warning
  if (length(snpdata$population)==0) {
    warning(paste("[WARNING] No population data for",s))
    return(NULL)
  }
#  if (s=="rs113216365") {
#    print("GOT HERE 0")
#    print(snpdata)
#    print(length(snpdata$population))
#  }
  # Get the allele string for the SNP and split it so we know
  # what the alleles are
  alleles=unlist(strsplit(unlist(snpdata$mappings)[['allele_string']],"/"))
#  print(alleles)
  # Now loop through all the data for the population
  for (i in snpdata$population) {
    # Not sure why I put this in here, I think it is just flattening everything out
    i=unlist(i)

    # Search for the population in the list of populations that has been passed
    # also search in the short names, if nether match then move on to the next one
#    if (s=="rs4520") {
#      print("GOT HERE")

#      print(pop)
#      print(i['population'])
#      print(poptranslate[pop==i['population'],shortpop])
#      print(i['population'] %in% pop)
#      print(poptranslate[pop==i['population'],shortpop] %in% pop)
#    }

    if (is.null(pop)==FALSE && i['population'] %in% pop==FALSE && (length(poptranslate[pop==i['population'],shortpop]) == 0 || (length(poptranslate[pop==i['population'],shortpop]) > 0 && poptranslate[pop==i['population'],shortpop] %in% pop==FALSE))) {
      next
    }

#    if (is.null(pop)==FALSE && i['population'] %in% pop==FALSE ) {
#      print("IN HERE")
#      next
#    }

    # Check to see if our population exists in the population translation table
    # if it doesn't then warn about it and add it so things do not fail below
    if (length(poptranslate[pop==i['population'],shortpop])==0L) {
      warning(paste("[WARNING] No short population name for",i['population'],"consider adding to the code"))
      poptranslate=rbindlist(list(poptranslate,data.table(i['population'],i['population'])) )
    }

    # Add the relevent details to a data.table that holds all of the
    # population data
    popdata=rbindlist(list(popdata,data.table(pop=i['population'],
                                              shortpop=poptranslate[pop==i['population'],shortpop],
                                              subid=i['submission_id'],
                                              snp=s,
                                              freq=as.numeric(i['frequency']),
                                              allele_count=as.numeric(i['allele_count']),
                                              allele=i['allele'],
                                              alleno="A0",
                                              count=0)))
  }

  # Make sure we have some data before we continue
  if (nrow(popdata)==0) {
    warning(paste("[WARNING] No population data") )
    return(NULL)
  }
#  if (s=="rs113216365") {
#    print("GOT HERE")
#    print(popdata)
#  }
#

  # Generate allele counts for each population/submission_id. The grouping
  # is by this as populations are not unique but pop/subid is unique.
  # I am doing this to count how many alleles are represented, this is
  # because monomorphic SNPs only have one allele in the data returned by
  # the rest API, we want to to have both represented
  # *** NOT SURE WHAT THE SUB ID IS AS IT IS ALLWAYS NA ***
  # EDIT 17/11/14: It seems to be a ss SNP name, it only seems to be populated
  # in non 1kg pops or in 1kg pre phase 1, So I will make all NA subid values
  # to be an empty string so when I subset and count later on I do not get unexpected
  # things happening
  popdata[is.na(subid)==TRUE,subid:=""]
#  if (all(is.na(popdata[,subid]))==FALSE) {
#    print(popdata)
#    stop("[FATAL] I didn't expect the sub ID to be populated")
#  }
  popdata[,count := as.double(.N),by=list(pop,snp,subid)]
#  popdata[,count := as.double(.N),by=list(pop,subid,snp)]

#  if (s=="rs181637183") {
#    print(popdata)
#    print(length(alleles))
#  }
#  print(popdata)
  # Get a list of monomorphics these will only have 1 count entry in the data.table
  # so will need to have the other entry (with the correct allele) inserted
  mono=popdata[count<length(alleles),]
#  if (s=="rs181637183") {
#    print("MONO")
    #  print(nrow(mono))
#    print(mono)
#  }
  # This is to place the "missing allele" records
  missing_alleles=data.table()
#if (s=="rs113216365") {
#  print("GOT HERE 2")
#  print(missing_alleles)
#}

  # If there are missing alleles then generate the missing allele with
  # frequeny of 0
  if (nrow(mono) > 0) {
    # Loop through all the monomorphic SNPs
    for (j in 1:nrow(mono)) {


      # Find the allele strings that are missing
      missing=alleles[!(alleles %in% mono[j,allele])]
#      if (s=="rs181637183") {
#        print(mono[j,])
#        print(missing)
#      }

      # Loop through all the missing allele strings and
      # generate a record for each one
      for (k in missing) {
#        print(data.table(mono[j,list(pop,shortpop,subid,snp,1-freq)],0,k,"A0",0))
        missing_alleles=rbindlist(
                                  list(
                                      missing_alleles,
                                      data.table(mono[j,list(pop,shortpop,subid,snp,1-freq)],0,k,"A0",0)
                                      )
                                  )
      }
    }
#    print("MISSING ALLELES")
#    print(missing_alleles)
    # Add the missing alleles to the population data
    popdata=rbindlist(list(popdata,missing_alleles))
  }
  # Do a recount just for good measure
  # EDIT I have disabled this as for some reason it is screwing up the counts
  popdata[,count := as.double(.N),by=list(pop,subid,snp)]
#  popdata[,count := as.double(.N),by=list(pop,snp)]
#  setkey(popdata,shortpop)
#  print(popdata)
#if (s=="rs113216365") {
#  print("GOT HERE 3")
#  print(popdata)
#}
  # Now we need to create an allele label column, i.e. A1 or A2
  # I will use the allele string to determine which allele should be
  # A1 and which should be A2
  for (i in 1:length(alleles)) {
    popdata$alleno[popdata$allele==alleles[i] & popdata$snp==s]=paste("A",i,sep="")
  }

#   # Do we want to subset some populations
#   if (is.null(args$pop)==FALSE) {
#     poprows=c()
#
#     # Loop through all the populations and grep the
#     # row numbers that match
#     #    print(args[['pop']])
#     for (p in args[['pop']]) {
#       matches=grep(p,popdata$pop)
#
#       # Warn if we can't find any matches
#       if (length(matches)==0) {
#         warning(paste("[WARNING] Population ",p," can't be found") )
#       }
#
#       # Store the matching row numbers and subset all at the end
#       poprows=c(poprows,matches)
#     }
#
#    # If none of the populations matched then this is a FATAL error
#    if (length(poprows)==0) {
#      stop("[FATAL] All Populations can't be found")
#    }
#
#    # Finally subset all the matching rows
#    popdata=popdata[unique(poprows),]
#  }

  setkey(popdata,pop,alleno)

#  if (return.maf==TRUE) {
#    return(popdata[,list(shortpop,subid,snp,freq,min(freq)),by=pop])
#  }

  return(popdata)
}



###############################################################################
# When given an ensembl gene id it builds a dataframe representing the        #
# structure of the gene, Transcripts,EXONs                                    #
###############################################################################
getGeneStructure=function(gene_id,
                          build="hg38",
                          query="lookup/id/%s?content-type=application/json;expand=1") {

  if( build == "hg38") {
    server="http://rest.ensembl.org"
    } else if( build == "hg19") {
      server="http://grch37.rest.ensembl.org"
    }

  # The columns in all the raw data that should be integers
  numeric_columns=c("start","end","strand")

  gene_query=sprintf(query,gene_id)
  full_query=paste(server,gene_query,sep="/")
#  gene_data=fromJSON(file=full_query)
  gene_data=runEnsemblQuery(full_query)
  # The fields that I am interested in for each gene,transcript,protein,exon
  gene_components=c("id","display_name","description","logic_name",
                    "species","assembly_name","biotype",
                    "seq_region_name","start","end","strand","db_type")

  transcript_components=c("id","display_name","logic_name","species",
                          "assembly_name","biotype","seq_region_name",
                          "start","end","strand","db_type","is_canonical")

  exon_components=c("id","assembly_name",
                    "seq_region_name","start","end","strand")

  protein_components=c("id","start","end","length")


  gene_info=c()
  for (i in gene_components) {
    if (is.null(gene_data[[i]])==TRUE) {
      gene_data[[i]]=NA
    }
    # Should really check that it exists before trying to add
    gene_info=c(gene_info,gene_data[[i]])
  }

  # Loop through all the transcript data and build a data.table from it
  transcript_data=data.table()
  protein_data=data.table()
  exon_data=data.table()
  # This loops through all the transcripts for the gene
  for (i in gene_data[['Transcript']]) {
    # This will hold the data for the transcript, the first element of the
    # vector will be the gene name as this will be used to merge to the gene
    # data at the end
    transcript_info=c(gene_info[1])

    # Loop through all the fields and get the data for each one
    for (tc in transcript_components) {
      transcript_info=c(transcript_info,i[[tc]])
    }

    # Some transcripts will be protein coding so we want to capture the
    # translation details as we can use these to determine UTR regions. The
    # first bit of info we store will be the transcript id
    protein_info=c(transcript_info[2])

    # If the transcript is protein coding then get the details for the
    # protein. Otherwise will will insert NA values
    if (is.null(i[['Translation']])==FALSE) {
      # Loop through the translation
      for (trc in protein_components) {
        protein_info=c(protein_info,i[['Translation']][[trc]])
      }
    } else {
      protein_info=c(protein_info,rep(NA,4))
    }

    # Store the protein data
    protein_data=rbindlist(list(protein_data,as.list(protein_info)))

    # Store the transcript data
    transcript_data=rbindlist(list(transcript_data,
                                   as.list(transcript_info)
                                   ))

    # Now grab the exon data for the gene
    for (e in i[['Exon']]) {
      exon_info=c(transcript_info[2])
      for (ec in exon_components) {
        exon_info=c(exon_info,e[[ec]])
      }

      # Bind to the exon data
      exon_data=rbindlist(list(exon_data,as.list(exon_info)))
    }
  } # end of looping through the transcripts

  gene_components=paste("gene",gene_components,sep="_")
  protein_components=paste("prot",protein_components,sep="_")
  transcript_components=paste("trans",transcript_components,sep="_")
  exon_components=paste("exon",exon_components,sep="_")

  setnames(protein_data,c("trans_id",protein_components))
  setnames(transcript_data,c("gene_id",transcript_components))
  setnames(exon_data,c("trans_id",exon_components))

  # merge the data on the transcript id
  transcript_data=merge(x=transcript_data,y=protein_data,by="trans_id")

  # merge the exon data on the transcript id
  transcript_data=merge(x=transcript_data,y=exon_data,by="trans_id")

  # Make the gene_data into a data.table and set the names
  gene_data=rbindlist(list(as.list(gene_info)))
  setnames(gene_data,gene_components)

  # Now merge the whole lot into the gene
  gene_data=merge(x=gene_data,y=transcript_data,by="gene_id")

  # Now I want to make sure that all the correct integer columns
  # are case to integers from characters
  cast_cols=c(paste("exon",numeric_columns,sep="_"),
              paste("gene",numeric_columns,sep="_"),
              paste("trans",numeric_columns,sep="_"),
              paste("prot",numeric_columns,sep="_") )

  return(castCols(gene_data,cast_cols,type="integer"))
}


###############################################################################
# Won't work as expected - this looks at the overlap endpoint                 #
###############################################################################
getOverlaps=function(gene_id,
                     server="h38",
                     query="overlap/id/%s?feature=gene;feature=transcript;feature=exon;content-type=application/json") {

    if( build == "hg38") {
    server="http://rest.ensembl.org"
    } else if( build == "hg19") {
      server="http://grch37.rest.ensembl.org"
    }

  gene_query=sprintf(query,gene_id)
  full_query=paste(server,gene_query,sep="/")
  #  print(full_query)
  gene_data=fromJSON(file=full_query)
  #  print(gene_data)

  # Create some data.tables to hold the
  # gene data
  # transcript data
  # exon data
  genes=data.table()
  transcripts=data.table()
  exons=data.table()

  for (n in 1:length(gene_data)) {
    #    print("================================================================")
    # Feature data
    fd=gene_data[[n]]
    #    print(fd)
    #    print(fd$feature_type)
    #    print(fd)

    if (fd$feature_type=="gene") {
      genes=rbindlist( list(genes,data.table("-",fd$id,fd$seq_region_name,fd$start,fd$end,fd$strand) ) )
    } else if (fd$feature_type=="transcript") {
      transcripts=rbindlist( list(transcripts,data.table(fd$Parent,fd$id,fd$seq_region_name,fd$start,fd$end,fd$strand) ) )
    } else if (fd$feature_type=="exon") {
      exons=rbindlist( list(exons,data.table(fd$Parent,fd$id,fd$rank,fd$seq_region_name,fd$start,fd$end,fd$strand) ) )
    }
    #    print(fd)
  }
  #  print(genes)

  setnames(genes,c("parent_id","ens_gene_id","gene_chr","gene_start","gene_end","gene_strand"))
  setkey(genes,"ens_gene_id")

  setnames(transcripts,c("ens_gene_id","ens_trans_id","trans_chr","trans_start","trans_end","trans_strand"))
  setkey(transcripts,"ens_gene_id")

  genes=merge(genes,transcripts,by="ens_gene_id")
  setnames(exons,c("ens_trans_id","ens_exon_id","rank","exon_chr","exon_start","exon_end","exon_strand"))
  setkey(exons,"ens_trans_id")

  genes=merge(genes,exons,by="ens_trans_id")
  genes=genes[ens_gene_id==gene_id,]
  #  print(exons)
  return(genes)
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
#                retry=TRUE
                return(NULL)
             },
             error=function(err) {
               write(sprintf("[ERROR] The query %s generated the following warning: %s",query,err$message),stdout())
#               retry=TRUE
               return(NULL)
             }#,
#             finally={
#               print("finally queried")
#               print(retry)
#               print(tries)
#             }
    )

    if (is.null(data)==TRUE) {
      retry=TRUE
    }
    tries=tries+1
  }
  return(data)
}


###############################################################################
# Uses the chr,pos,allele and strand information to get a VEP prediction for  #
# a variation                                                                 #
###############################################################################
getVepSnp = function(chr, pos, allele, build=38,
                   name = NULL,
                   query = "vep/human/region/%i:%i-%i/%s?content-type=application/json",
                   allow.tries = 2) {
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
  r=httr::GET(paste(server, vep_query, sep = "/"), content_type("application/json"))
  vep_data=jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))

  if(!("error" %in% names(vep_data))) {
    return(vep_data)
  }

}


###############################################################################
# A functoion to get SNPs                                                     #
###############################################################################
getSnps=function(snpid,
                 build="hg38",
                 query="variation/human/%s?content-type=application/json",
                 allow.tries=2) {

    if( build == "hg38") {
    server="http://rest.ensembl.org"
    } else if( build == "hg19") {
      server="http://grch37.rest.ensembl.org"
    }

  all_snp_data=data.frame()
  # Build the JSON string. If the SNP does not exist then it will
  # crash and die. I need to make it exit gracefully at some point
  snp_query=sprintf(query,snpid)
  full_query=paste(server,snp_query,sep="/")
#  print(full_query)
  snp_data=runEnsemblQuery(full_query,allow.tries=allow.tries)
#  print(snp_data)
  if (is.null(snp_data)==TRUE) {
    return(snp_data)
  }
#  print(names(snp_data))
  #"synonyms","evidence"
  snp_cols=c("source","name","MAF","ambiguity","var_class","ancestral_allele","most_severe_consequence")

  base_data=c()
  for (c in snp_cols) {
    sd=""
    if (is.null(snp_data[[c]])==FALSE) {
      sd=snp_data[[c]]
    }
    base_data=c(base_data,sd)
#    base_data=c(base_data,snp_data[[c]])
  }


  mapping_cols=c("location","assembly_name","seq_region_name","start","end","strand")
  # Now loop through all the mappings adding the location information
  for (m in snp_data$mappings) {
    mapping_data=c()
#    print(base_data)
#    print(m)
    for (mc in mapping_cols) {
      md=""
      if (is.null(m[[mc]])==FALSE) {
        md=m[[mc]]
      }
      mapping_data=c(mapping_data,md)
#      mapping_data=c(mapping_data,m[[mc]])
    }
#    print(c(base_data,mapping_data))
#    print(t(data.table(c(base_data,mapping_data))))
    all_snp_data=rbindlist(list(all_snp_data,as.list(c(base_data,mapping_data))))
  }


#  print(c(snp_cols,mapping_cols))
  setnames(all_snp_data,names(all_snp_data),c(snp_cols,mapping_cols))
  all_snp_data[,most_severe_consequence := gsub("\\s+","_", tolower(most_severe_consequence)) ]
  return(all_snp_data)
#  base_data=
}


###############################################################################
# This adds a priority column to a SNP list that have been fetched (or at     #
# least) is in the same format as the data returned by getSnps()              #
#                                                                             #
###############################################################################
assignSnpPriority = function(variations,so_terms) {
#  setkey(so_term,"so_term")
#  print(names(variations))
#  print(names(so_terms))

  # Re-name the so_term column name to most_severe_consequence
  setnames(so_terms,"so_term","most_severe_consequence")
  so_terms=so_terms[,list(most_severe_consequence,priority)]

  setkey(so_terms,"most_severe_consequence")
  setkey(variations,"most_severe_consequence")

  variations=merge(variations,so_terms,all.x=TRUE)

  # Anything that we don't know about gets the lowest priority
  variations[is.null(priority)==TRUE,priority:=max(so_terms[,priority])+1]
  return(variations)
}


###############################################################################
# Will return ensembl gene ids overlapping a region                           #
###############################################################################
getGenesOverlapRegion = function(chr,
                                 start,
                                 end,
                                 build="hg38",
                                 species="human",
                                 query="overlap/region/%s/%s:%i-%i?feature=gene;content-type=application/json",
                                 allow.tries=2) {

  if( build == "hg38") {
    server="http://rest.ensembl.org"
    } else if( build == "hg19") {
      server="http://grch37.rest.ensembl.org"
    }


  chr=as.character(chr)
  if (chr=="23") {chr="X"}
  if (chr=="24") {chr="Y"}
  if (chr=="25") {chr="MT"}
  all_snp_data=data.frame()
  # Build the JSON string. If the SNP does not exist then it will
  # crash and die. I need to make it exit gracefully at some point
  region_query=sprintf(query,species,chr,start,end)
  full_query=paste(server,region_query,sep="/")
  #  print(full_query)
  region_data=runEnsemblQuery(full_query,allow.tries=allow.tries)
  #  print(snp_data)
  if (is.null(region_data)==TRUE) {
    return(region_data)
  }

  # Expcols
  gene_cols=c("source",
              "logic_name",
              "version",
              "feature_type",
              "external_name",
              "description",
              "assembly_name",
              "biotype",
              "start",
              "end",
              "seq_region_name",
              "strand",
              "id")

  all_gene_data=data.table()

  # If we get here then we have some data
  for (i in region_data) {
#    print(names(i))
    # To hold the data for the current gene
    gene_data=c()

    # Loop through all the columns that we expect
    for (n in gene_cols ) {
      # If there is no data for that column then make it NA
      if (is.null(i[[n]])==TRUE) {
        i[[n]]=NA
      }

      # Add to the data vector
      gene_data=c(gene_data,i[[n]])
    }

    # Finally when we have built the data.vector add it to the total
    # list of genes we have in the region
    all_gene_data=rbindlist(list(all_gene_data,as.list(gene_data)))

  }

  if (nrow(all_gene_data)>0) {
    # Add some column names fot the gene data and return
    setnames(all_gene_data,names(all_gene_data),gene_cols)
  }

  return(all_gene_data)
}



###############################################################################
# Return the gene that is nearest to a single genomic point                   #
###############################################################################
getNearestGene = function(chr,
                          pos,
                          test.flank.bp=100000,
                          test.flank.inc=100000,
                          biotypes=NULL,
                          species="human",
                          build="hg38",
                          allow.tries=2) {

  if( build == "hg38") {
    server="http://rest.ensembl.org"
    } else if( build == "hg19") {
      server="http://grch37.rest.ensembl.org"
    }

  #  chr=as.numeric(chr)
  pos=as.numeric(pos)
  test.flank.bp=as.numeric(test.flank.bp)
  test.flank.inc=as.numeric(test.flank.inc)
  #  print(chr)
  #  print(pos)
  # The is the maximum flank test size to look for genes
  # if it is beyond this and no genes have been found then
  # we will return NULL and a warning
  max.test.size=2000000

  # Will hold the genes nearest the genomic position
  genes=data.table()

  repeat {
    startpos=max(c(1,pos-test.flank.bp))
    endpos=pos+test.flank.bp
    # print(sprintf("[INFO] Searching %i:%i-%i",chr,startpos,endpos))
    # print(startpos)
    # print(endpos)

    # Get the genes overlapping the defined region
    genes=getGenesOverlapRegion(chr,
                                startpos,
                                endpos,
                                build=build,
                                species=species,
                                allow.tries=allow.tries)
    # print(genes)
    if (nrow(genes) > 0 && is.null(biotypes)==FALSE) {
      genes=genes[biotype %in% biotypes,]
    }

    # If we have some genes, then find out which one is closest. If threre
    # is only one gene then by definition it must be the closest
    if (nrow(genes)>0) {
      # Calculate the distances between our genomic position and the start and
      # end of each gene
      genes[,enddist:=abs(pos-as.numeric(end))]
      genes[,startdist:=abs(pos-as.numeric(start))]

      # Now make sure the genomic position is not in the gene in which case
      # I will set the distance to -1
      genes[pos<=as.numeric(end) & pos>=as.numeric(start),enddist:=-1]
      genes[pos<=as.numeric(end) & pos>=as.numeric(start),startdist:=-1]

      # Calculate the lowest distance between the start and end distance
      genes[,mindist:=min(.SD[,]),.SDcols=c("startdist","enddist"),by=id]

      # Now get the gene with the overall minimum distance
      genes=genes[which.min(mindist),]
    }

    # See if we need to exit the loop, Keep on testing while we have no genes
    # and are below our max test size
    if (test.flank.bp >= max.test.size || nrow(genes)>0) {
      break
    }

    # If we get to here we will test a wider flank
    test.flank.bp=test.flank.bp+test.flank.inc
  }

  # Issue a warning if we have 0 genes
  if (nrow(genes)==0) {
    biostr=""
    if (length(biotypes)>0) {
      biostr=sprintf(" with biotypes: %s",paste(biotypes,collapse=","))
    }

    warning(sprintf("[WARNING] There are 0 genes in the region %i:%i-%i %s",chr,startpos,endpos,biostr))
#    print(genes)
  }

  if (nrow(genes)>0) {
    genes[,search_chr:=chr]
    genes[,search_pos:=pos]
  }

  return(genes)
}
