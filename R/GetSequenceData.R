#==============================================================
# The polaRRtools package
#       polaRRtools is a collection of user-friendly tools to download explore and analize 
#       High Throughput Amplicon sequencing data for Biodiversity research, using R. 
#       The aim is to lower the treshold for high quality, reproducible and accessible 
#       analysis of nucleotide data in the versatile R-environment, combined with data re-use.
#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2019-09-20)
#
#   This is version includes utilities to:
#     - download sequence data from INSDC (NCBI's SRA to be more precise)
#     - download meta data and environmental data fron INSDC (NCBI's SRA to be more precise)
#     - Quality controll meta data and environmental data
# 
#==============================================================

#library(devtools)
#library(usethis)
TermsLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/TermsLibrary.csv", 
                      header=TRUE)
TermsSyn <- sapply(as.character(TermsLib$synonyms), function(x){strsplit(x, ";")})
names(TermsSyn) <- TermsLib$name 
#usethis::use_data(MIxS_Voc, internal = FALSE)

ENA_checklistAccession<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/ENA_checklistAccession.csv", 
                   header=TRUE)

TaxIDLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/TaxIDLibrary.csv", 
                                 header=TRUE, row.names=1)

require(tidyverse)
#--------------------------------------------------------------
# 0. setting classes
#--------------------------------------------------------------
setClass("metadata.MIxS", slots=list(data="data.frame", #a dataframe with the metadata (samples are rows, variables are columns)
                                     section="character", #a character vector containing a MIxS section for each term (variable)
                                     units="character", #a character vector containing a unit for each variable
                                     type="character", #versatile (loose following of MIxS) or strict.MIxS (following all the MIxS rules)
                                     env_package="character" #the environmental package. Can be not_specified or multiple_packages.
                                     )
         )

setMethod("show",
          "metadata.MIxS",
          function(object) {
            N <- nrow(object@data)
            C <- ncol(object@data)
            cat(paste("a metadata.MIxS class object with", as.character(N), "samples and",
                      as.character(C), "variables\n"))
            }
          )

check.valid.metadata.MIxS <- function(d){
  valid <- TRUE
  #class must be metadata.MIxS 
  if(!class(d)=="metadata.MIxS"){
    valid <- FALSE
  }else{
    dcol <- ncol(d@data)
    drow <- nrow(d@data)
    #The meta must have at least 1 sample (row) and 1 variable (column)
    if(dcol<=0 & drow<=0){
      valid <- FALSE
    }
    # each variable must have an associated unit
    if(length(d@units) < dcol | 
       length(d@section) < dcol |
       length(d@units) != length(d@section)){
      valid <- FALSE
    }
    # only 2 types allowed: strict MIxS or a versatile form that not strictly follows the MIxS rules
    if(!d@type %in% c("versatile", "strict.MIxS")){
      valid <- FALSE
    }
    # the environmental package must be valid
    if(!d@env_package %in% c("air", "built_environment", "host_associated", "human_associated",
                             "human_gut", "human_oral", "human_skin", "human_vaginal", 
                             "microbial_mat_biofilm", "miscellaneous_natural_or_artificial_environment",
                             "plant_associated", "sediment", "soil", "wastewater_sludge", "water",
                             "not_specified", "multiple_packages")){
      valid <- FALSE
    }
  }
  return(valid)
}


#--------------------------------------------------------------
# 1. Tools for formatting data
#--------------------------------------------------------------
process.metadata <- function(metadata = NA, add_to = NA, strict.MIxS = FALSE, 
                             out.format="metadata.MIxS", ask.input=TRUE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description  process.metadata takes a dataframe with raw metadata downloaded from an INSDC database (see get.BioProject.metadata.INSDC and get.sample.attributes.INSDC) and performs a basis Quality Controll. (see further for details)
  #' 
  #' @usage process.metadata(metadata, add_to, strict.MIxS = FALSE)
  #' 
  #' @param metadata data.frame. The raw metadata downloaded from INSDC to be cleaned up. Rows are samples, columns variables
  #' @param add_to a metadata.MIxS object. An already present dataset with quality-comtrolled metadata. must be formatted as metadata.MIxS to ensure the correct input format of the data.
  #' @param strict.MIxS boolean. If TRUE, variable names in the output will only be MIxS terms. Any variable with no MIxS counterpart is discarded. If FALSE, than variables that have no appropriate MIxS term wil be added as extra "miscellaneous" collumns. Default is FALSE
  #' @param out.format character. The format (object class) of the output, either metadata.MIxS or data.frame. Default metadata.MIxS
  #' @param ask.input boolean. If TRUE, console input will be requested to the user when a problem occurs. Default TRUE
  #' @details Any sequencing project typically has important additional data associated with it.
  #' This goes from laboratory protocols, sequencing platform settings or environmental measurements.
  #' While (most of) these data get deposited on the INSDC databases alongside the sequences, the
  #' quality and format of this data is typically determined by the data provider rather than international
  #' standards, like MIMARKS, MIxS and Darwin Core.
  #' Therfore, the process.metadata function was develloped to sort through these metadata, and
  #' perform a basic quality controll, correcting the most common mistakes, like incorrectly 
  #' formatting the geographic coordinates, typing errors in variable names, etc. To do this, process.metadata 
  #' makes use of a build-in dictionary of (MIxS) terms and their synonyms (that is: spelling errors, 
  #' writing differences, true synonyms,...). Note that it is possible some terms are not recognized. 
  #' In that case contact the author to update the dictionary in the upcomming version.
  #' 
  #' The metadata argument is typically a dataframe generated using the get.BioProject.metadata.INSDC or 
  #' get.sample.attributes.INSDC functions. 
  #' If the argument add_to is not NA or an empty dataframe, the output will be aded to fit the add_to dataframe input.
  #' For any input, however, rows must be samples, columns variables.
  #' 
  #' process.metadata adresses the following issues for the metadata argument:
  #' 1) variable names are converted to accepted MIxS standard terms if typos or incorrect notations are found
  #' 2) Dates, if recognized, are returned in the YYYY-MM-DD format.
  #' 3) Latitudes and longitudes are returned in decimal notation. Degrees are converted.
  #' 4) length/depth measurements, if recognized, deci- and centimeters are converted to meter
  #' 
  #' @seealso get.BioProject.metadata.INSDC, get.sample.attributes.INSDC
  #' 
  #' @return a dataframe that is compatible with the MIxS standard, or that stricly follows the MIxS standard
  #' 
  #' @example metadataQC_PRJNA369175 <- process.metadata(metadata = metadata_PRJNA369175, strict.MIxS = FALSE)
  
  warningmessages<-c()
  
  if(!is.data.frame(metadata)){
    stop("The input must be a dataframe, with samples as rows/ variables as columns")
  }
  if(!is.na(add_to) && !check.valid.metadata.MIxS(add_to)){
    stop("The input for the add_to argument must be a metadata.MIxS object to ensure correct merging of the datasets")
  }
  if(!is.na(add_to) && out.format=="data.frame"){
    stop("when datasets need to be merged, the output must be formatted as metadata.MIxS. If necessary, the metadata can be extracted from this object as a dataframe by accessing the \"data\"-slot")
  }
  if(out.format=="data.frame"){
    warningmessages <- multi.warnings("units are not included in the output when the out.format argument is set to data.frame. To include units in the output set out.format to metadata.MIxS", warningmessages)
  }
        
  # remove empty columns
  metadata <- metadata[,colSums(is.na(metadata)) < nrow(metadata)] 
  
  # make an empty output file to fill along the way
  New_metadata <- data.frame(row.names=rownames(metadata))
  
  # 1. looking for the original sample name
  item_shared <- intersect(TermsSyn["original_name"][[1]], colnames(metadata))
  if(length(item_shared)>=1){ 
    likely_sampNames <- item_shared[1] #order of item_shared is in decreasing likelyness
    if(ask.input){
      cat(paste("The original sample names could be in the \"",likely_sampNames,"\" column.\n", sep="")) 
      cat(paste("The first five names in this column are ", paste(metadata[1:5,likely_sampNames], collapse="; "), " ...\nDoes this seems correct? (y/n)\n", sep="")) 
      doNext <- readline() 
      if(doNext %in% c("y", "Y", "yes", "YES", "Yes")){
        New_metadata$original_name <- as.character(metadata[,likely_sampNames])
        warningmessages <- multi.warnings(paste("assumed the \"",likely_sampNames,"\" column contained the original sample names", sep=""), warningmessages)
        tryCatch({
          New_metadata$INSDC_SampleID <- rownames(New_metadata)
          rownames(New_metadata) <- New_metadata$original_name
          },
                 error=function(x){
                   warningmessages <- multi.warnings("duplicate or missing original sample names", warningmessages)
                 })
      }else if(!doNext %in% c("n", "N", "no", "NO", "No")){
        stop("incorrect input... only yes or no allowed.")
      } 
    }else{
      New_metadata$original_name <- as.character(metadata[,likely_sampNames])
      warningmessages <- multi.warnings(paste("assumed the \"",likely_sampNames,"\" column contained the original sample names", sep=""), warningmessages)
    }
  }else{
    if(.row_names_info(Heind)>0){#rownames provided by the user might be the original names
      New_metadata$original_name <- row.names(metadata)
      warningmessages <- multi.warnings("no original sample names found, used the rownames instead", warningmessages)
    }else{
      warningmessages <- multi.warnings("no original sample names found", warningmessages)
    }
  }
  

  
  # 2. some basic info from insdc
  TermsSyn_insdc<-TermsSyn[as.character(TermsLib[TermsLib$insdc==TRUE,]$name)]
  for(item in names(TermsSyn_insdc)){
    item_shared <- intersect(TermsSyn_insdc[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    }
  }
  
  # 3. dealing with latitude-longitude, and it's many possible formats...
  TermsSyn_latlon<-TermsSyn[as.character(TermsLib[TermsLib$name=="lat_lon",]$name)]
  TermsSyn_lat<-TermsSyn[as.character(TermsLib[TermsLib$name=="DwC_decimalLatitude",]$name)]
  TermsSyn_lon<-TermsSyn[as.character(TermsLib[TermsLib$name=="DwC_decimalLongitude",]$name)]
  
  QClatlon <- dataQC.LatitudeLongitudeCheck(metadata, 
                                            latlon.colnames=list(TermsSyn_latlon[[1]],
                                                                 TermsSyn_lat[[1]],
                                                                 TermsSyn_lon[[1]]))
  warningmessages<-c(QClatlon$warningmessages, warningmessages)
  New_metadata$lat_lon <- QClatlon$values
  New_metadata$DwC_decimalLatitude <- sapply(QClatlon$values, function(x){strsplit(x, " ")[[1]][1]})
  New_metadata$DwC_decimalLongitude <- sapply(QClatlon$values, function(x){strsplit(x, " ")[[1]][2]})
  
  New_metadata[is.na(New_metadata)] <- ""


  # 4. dealing with the collection date, and putting it in the YYYY-MM-DD format
  TermsSyn_date<-TermsSyn[as.character(TermsLib[TermsLib$name=="collection_date",]$name)]
  QCDate <- dataQC.dateCheck(metadata, TermsSyn_date[[1]])
  warningmessages<-c(QCDate$warningmessages, warningmessages)
  if(length(QCDate$values)==nrow(metadata)){
    New_metadata$collection_date <- QCDate$values
  }
  
  
  # 5. the core MIxS terms
  TermsSyn_MIxS <- TermsSyn[as.character(TermsLib[TermsLib$core>0,]$name)]
  TermsSyn_MIxS <- TermsSyn_MIxS[!names(TermsSyn_MIxS) %in% c("lat_lon", "collection_date")]
  
  for(item in names(TermsSyn_MIxS)){
    item_shared <- intersect(TermsSyn_MIxS[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    } else if(length(item_shared)>1){
      New_metadata[,item] <- metadata[,item_shared[[1]]]
    }
  }
  
  # 6. the package MIxS terms (just one package allowed for strict.MIxS, otherwise multiple packages allowed)
  # 6.1 find the best package/ask user if no package was specified
  if(!"env_package" %in% colnames(New_metadata)){
    if(ask.input){
      cat("No env_package was specified.\nPlease specify what to do next:\n1) Make an educated guess based on the data\n2) Ask user for the package name\n3) stop executing\n(type 1, 2 or 3)\n") 
      doNext <- readline() 
      if(doNext==1){
        env_package <- dataQC.guess.env_package.from.data(metadata)
        warningmessages<-c(warningmessages, env_package$warningmessages)
        env_package <- env_package$values
      }else if(doNext==2){
        if(strict.MIxS){
          cat("One MIxS environmental package is required for a strict.MIxS output. Found none or more than one.\nPlease provide a single package or type none to turn off strict.MIxS\nThe choices are: air, built_environment, host_associated, human_associated,human_gut,\nhuman_oral, human_skin, human_vaginal,microbial_mat_biofilm,\nmiscellaneous_natural_or_artificial_environment,\nplant_associated, soil, sediment, wastewater_sludge, water\n") 
        } else{
          cat("Please provide a single MIxS environmental package\nThe choices are: air, built_environment, host_associated, human_associated,human_gut,\nhuman_oral, human_skin, human_vaginal,microbial_mat_biofilm,\nmiscellaneous_natural_or_artificial_environment,\nplant_associated, soil, sediment, wastewater_sludge, water\n") 
        }
        env_package <- readline() 
        if(env_package=="none"){
          strict.MIxS<-FALSE
        } else if(!env_package %in% colnames(TermsLib)){
          stop("incorrect environmental package provided. Be sure to use underscores and lowercase letters")
        } else{
          env_package <- rep(env_package, nrow(metadata))
        }
      }else if(doNext==3){
        stop("you chose to interrupt execution.")
      }else{
        stop("incorrect input. Interrupted execution.")
      }
    }else{
      env_package <- dataQC.guess.env_package.from.data(metadata)
      warningmessages<-c(warningmessages, env_package$warningmessages)
      env_package <- env_package$values
    }
  }else{# 6.2 check if package is valid or if the ENA checklist number needs to be converted to a package
    if(length(setdiff(unique(New_metadata$env_package), ENA_checklistAccession$env_package))>0 |
       length(setdiff(unique(New_metadata$env_package), ENA_checklistAccession$ena_package))>0){
      #no correct package name, check if it is one of the ENA accession numbers
      if(sum(grepl("ERC", unique(New_metadata$env_package)))==length(unique(New_metadata$env_package))){
        #convert all the ENA checklist accession numbers to MIxS packages
        ENA_checklistAccession$ena_checklist_accession
        for(i in 1:nrow(New_metadata)){
          pk <- ENA_checklistAccession[ENA_checklistAccession$ena_checklist_accession == New_metadata[i,]$env_package,]$env_package
          if(length(pk)>0 && ENA_checklistAccession[ENA_checklistAccession$env_package==pk,]$MIxS==TRUE){
            New_metadata[i,]$env_package <- pk
          } else{
            New_metadata[i,]$env_package <- ""
            warningmessages<-multi.warnings("For some samples, the given environmental package did not correspond to a MIxS package", warningmessages)
          }
        }
      }else{
        warningmessages<-multi.warnings("It seems the way the environmatla package is written is not allowed, better correct", warningmessages)
      }
    }
  }
  
  
  if(is.null(env_package) & !strict.MIxS){
    warningmessages<-multi.warnings("No env_package could be inferred", warningmessages)
  } else{
    New_metadata$env_package <- env_package
  }
  
  if(strict.MIxS){
    if(is.null(env_package) || length(unique(env_package)) != 1){
      if(ask.input){
        cat("Found none or more than one environmental packages. Only one allowed when strict.MIxS is TRUE\nContinue by setting strict.MIxS to FALSE? (y/n)\n") 
        ctn <- readline() 
        if(ctn %in% c("y", "Y", "yes", "YES", "Yes")){
          strict.MIxS<-FALSE
        }else if(ctn %in% c("n", "N", "no", "NO", "No")){
          stop("Could not output metadata strictly following the MIxS rules.")
        } else{
          stop("incorrect input... only yes or no allowed.")
        }
      }else{
        stop("Found none or more than one environmental packages. Only one allowed when strict.MIxS is TRUE")
      }
    }else{
      TermsSyn_MIxSpackage<-TermsSyn[as.character(TermsLib[TermsLib$official_MIxS==TRUE & TermsLib[,env_package]>0,]$name)]
    }
  }else{
    TermsSyn_MIxSpackage<-TermsSyn[as.character(TermsLib[TermsLib$official_MIxS==TRUE & TermsLib$core==0,]$name)]
  }
    
  for(item in names(TermsSyn_MIxSpackage)){
    item_shared <- intersect(TermsSyn_MIxSpackage[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    } else if(length(item_shared)>1){
      New_metadata[,item] <- metadata[,item_shared[[1]]]
    }
  }
  
  #check if minimal terms are present
  if(strict.MIxS){
    TermsSyn_MIxSpackage_minimal <- c(names(TermsSyn[as.character(TermsLib[TermsLib[,env_package]==2,]$name)]),
                                      names(TermsSyn[as.character(TermsLib[TermsLib[,"core"]==2,]$name)]))
    
    if(length(intersect(names(TermsSyn_MIxSpackage_minimal), names(TermsSyn_MIxSpackage))) != length(names(TermsSyn_MIxSpackage_minimal))){
      if(ask.input){
        cat("Some of the required terms for strict.MIxS are missing...\nContinue with strict.MIxS turned off? (y/n)\n") 
        ctu <- readline() 
        if(ctu %in% c("y", "Y", "yes", "YES", "Yes")){
          strict.MIxS<-FALSE
        }else if(ctu %in% c("n", "N", "no", "NO", "No")){
          stop("Could not output metadata strictly following the MIxS rules.")
        } else{
          stop("incorrect input... only yes or no allowed.")
        }
      }else{
        stop("Some of the required terms for strict.MIxS were missing.")
      }
    }
    
  }
  

  # 7. any other additionalinformation
  TermsSyn_add<-TermsSyn[as.character(TermsLib[TermsLib$insdc==FALSE & TermsLib$official_MIxS==FALSE,]$name)]
  for(item in names(TermsSyn_add)){
    item_shared <- intersect(TermsSyn_add[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    }
  }
  
  # 8. some additional quality controll
  if("investigation_type" %in% colnames(New_metadata)){
    New_metadata$investigation_type <- sapply(New_metadata$investigation_type, FUN=function(x){
      if(x %in% c("AMPLICON", "amplicon", "metabarcode")){x<-"mimarks-survey"
      }else if(x=="WGS"){x<-"metagenome"}
      return(x)})
  }else{
    if(ask.input){
      cat("No investigation_type was found...\nPlease provide an investigation_type. Common ones include mimarks-survey or metagenome. Type n to ignore.\n") 
      invtype <- readline() 
      if(! invtype %in% c("n", "N")){
        New_metadata$investigation_type <- rep(invtype, nrow(New_metadata))
      }
    }
  }
  if("target_gene" %in% colnames(New_metadata)){
    New_metadata$target_gene <- sapply(New_metadata$target_gene, function(x){
      if(grepl("16S", x)){x<-"16S ssu rRNA"
      }else if(grepl("18S", x)){x<-"16S ssu rRNA"
      }else if(grepl("ITS", x)){x<-"ITS"}})
  }
  if("specific_host" %in% colnames(New_metadata)){
    host_val <- setdiff(unique(New_metadata$specific_host), c(NA, "NA", "-", "not applicable"))
    if(identical(host_val, character(0))){
      New_metadata <- New_metadata[,!colnames(New_metadata) %in% "specific_host"]
    }
  }
  
  # QC on length/depth/size measurements
  # if multiple units found: do nothing
  # if a string found that is not one of the expected units: leave it and work with the rest
  size_related <- c("depth", "elev", "alt_elev", "filter_size", "tot_depth_water_col")
  items_shared <- intersect(size_related, colnames(New_metadata))
  QC_units<-c()
  if(length(items_shared)>0){
    alternative_units<-data.frame(name=items_shared, unit_full=rep(NA, length(items_shared)))
    possible_units <- c("nm", "um", "??m", "mm", "cm", "dm", "m", "km")
    for(name_unit in items_shared){
      vals <- as.vector(as.character(New_metadata[,colnames(New_metadata) %in% name_unit]))
      names(vals)<-row.names(New_metadata)
      #extract the units
      val_units <- unlist(lapply(vals, function(v){gsub("[0-9]|\\.|,|-", "", v)} ))
      val_units <- intersect(unique(val_units), possible_units) #danger in this step: discards any text that is not an expected unit
      vals<-unlist(lapply(vals, function(v){gsub(",", ".", v)} ))
      if(length(val_units)==1){
        for(v in 1:length(vals)){
          if(grepl(val_units[1], vals[v])){
            vals[v]<-gsub(val_units[1], "", vals[v])
          }
        }
        QC_units[name_unit]<-val_units[1]
        New_metadata[,colnames(New_metadata)==name_unit] <- c(vals[rownames(New_metadata)])
      }
    }
  }
  
  # 9. finalizing and formatting the output
  # 9.1 the units
  New_metadata_units <- c()
  for(i in 1:ncol(New_metadata)){
    New_metadata_units[colnames(New_metadata)[i]] <- as.character(TermsLib[TermsLib$name %in% colnames(New_metadata)[i],]$expected_unit)
  }
  if(length(QC_units)>0){
    for(u in 1:length(QC_units)){
      New_metadata_units[names(QC_units)[u]] <- QC_units[u]
    }
  }
  
  # 9.2 if strict.MIxS: concatenate all non-MIxS terms in the misc_param term
  if(strict.MIxS){
    MIxS_terms <- setdiff(as.character(TermsLib[TermsLib$official_MIxS==TRUE,]$name), "misc_param")
    misc_units <- New_metadata_units[which(!colnames(New_metadata) %in% MIxS_terms)]
    New_metadata_units <- New_metadata_units[which(colnames(New_metadata) %in% MIxS_terms)]
    misc_metadata <- New_metadata[,!colnames(New_metadata) %in% MIxS_terms, drop=FALSE]
    New_metadata <- New_metadata[,colnames(New_metadata) %in% MIxS_terms, drop=FALSE]
    if(ncol(misc_metadata)>0){
      for(cl in 1:ncol(misc_metadata)){
        if(!grepl("alphanumeric", misc_units[cl])){
          un <- paste("(", misc_units[cl], ")", sep="")
        }else{
          un<-""
        }
        misc_metadata[,cl]<-paste(colnames(misc_metadata)[cl], ":", misc_metadata[,cl], un, sep="")
      }
      New_metadata$misc_param<-apply(misc_metadata, 1, function(x) paste(x, collapse = ";"))
      New_metadata_units["misc_param"] <- "alphanumeric"
    }
  }
  
  # 9.3 the section
  New_metadata_section <- c()
  for(i in 1:ncol(New_metadata)){
    New_metadata_section[colnames(New_metadata)[i]] <- as.character(TermsLib[TermsLib$name %in% colnames(New_metadata)[i],]$section)
  }
  
  # 9.4 env_package
  if("env_package" %in% colnames(New_metadata)){
    env_package <- unique(New_metadata$env_package)
    if(length(env_package)>1){
      env_package <- "multiple_packages"
    } else if(is.na(env_package) | is.null(env_package)){
      env_package <- "not_specified"
    }
  } else{
    env_package <- "not_specified"
  }

  # 9.5 warning messages
  if(length(warningmessages)>0){
    for(i in 1:length(warningmessages)){
      warningmessages[i] <- paste(i, warningmessages[i], sep=". ")
    }
    warningmessages <- c("Please consider the following warning messages carefully before proceding:", warningmessages)
    warning(paste(warningmessages, collapse='\n'))
  }
  
  # 9.6 convert to the right output format (data.frame or metadata.MIxS)
  if(out.format=="metadata.MIxS"){
    if(strict.MIxS){
      New_metadata <- new("metadata.MIxS",
                          data   = New_metadata,
                          section = New_metadata_section,
                          units      = New_metadata_units,
                          env_package = env_package,
                          type = "strict.MIxS"
      )
    }else{
      New_metadata <- new("metadata.MIxS",
                          data   = New_metadata,
                          section = New_metadata_section,
                          units      = New_metadata_units,
                          env_package = env_package,
                          type = "versatile"
      )
    }
    if(!is.na(add_to)){
      New_metadata <- combine.data(New_metadata, add_to)
    }
  }
    return(New_metadata)
}

#--------------------------------------------------------------
# 2. Tools for downloading data
#--------------------------------------------------------------
get.BioProject.metadata.INSDC <- function(BioPrjct=NA, just.names=FALSE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description get.BioProject.metadata.INSDC downloads a minimal set of sample metadata of all samples ("Runs") within a BioProject (argument BioPrjct) from the INSDC databases.
  #' 
  #' @usage get.BioProject.metadata.INSDC(BioPrjct, just.names=FALSE)
  #' 
  #' @param BioPrjct  A chracter string. A single Bioproject ID.
  #' @param just.names  Boolean. If TRUE, only the INSDC sample names are returned, else all data are returned. default FALSE
  #' 
  #' @details BioProjects combine all biological nucleotide sequence data related to a single initiative, 
  #' originating from a single organization. With each sample ("Run") within a BioProject, there is additional
  #' data associated that is crucial to the correct interpretation of the nucleotide sequence data, but is
  #' not automatically downloaded along with it. 
  #' The get.BioProject.metadata.INSDC function will fetch the most basic metadata of a BioProject from the INSDC
  #' repositories to complete the nucleotide sequence dataset, using E-utils API function of NCBI.
  #' These basic metadata typically include Run number, relsease date, load date, spots, bases, av_MB and download path
  #' 
  #' Note that the data returned by get.BioProject.metadata.INSDC does not include all the metadata associated
  #' with a BioProject. Other information, like coordinates, sampling dates or environmental measurements may also
  #' be available, but require the user to register at NCBI and request a personal API-key (this is required by
  #' NCBI to acces their data since 2017)
  #' 
  #' The complete set of additional data can be downloaded using the get.sample.attributes.INSDC 
  #' function, given a user-specified API-key. see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
  #' 
  #' This function is only valid for BioProjects that are registered on NCBI's SRA database.
  #' The BioPrjct argument should be a character string of an existing BioProject number.
  #' 
  #' @seealso get.sample.attributes.INSDC
  #' 
  #' @return if get.BioProject.metadata.INSDC(just.names=FALSE) (default) a data.frame with n rows an m columns is returned, n being the number of samples in the BioProject,
  #' and m being the number of variables found. If get.BioProject.metadata.INSDC(just.names=TRUE), a character vector of length n is returned with the sample numbers ("Run numbers", SRR numbers)
  #' 
  #' @references Sayers, E. (2009) The E-utilities In-Depth: Parameters, Syntax and More, https://www.ncbi.nlm.nih.gov/books/NBK25499/
  #' 
  #' @example metadata_PRJNA369175 <- get.BioProject.metadata.INSDC(BioPrjct="PRJNA369175", just.names=FALSE)
  
  if(! is.character(BioPrjct) | c(NULL,NA) %in% BioPrjct | length(BioPrjct)>1){
    stop("incorrect BioProject provided. 
         Expected input is a single Bioproject ID as a character string.")
  }
  
  sra_url <- paste("http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=", BioPrjct, sep="")
  tmpFile <- tempfile()
  download.file(sra_url, destfile = tmpFile, "wget", quiet = TRUE)
  if(readLines(tmpFile, n=1)==""){
    stop(paste("No metadata found for ", BioPrjct, ". Could not download the project metadata.", sep="")) 
  } else{
    RawMetadata <- read.csv(tmpFile, header=TRUE)
  }
  file.remove(tmpFile)
  if(just.names){
    return(c(as.character(RawMetadata$Run)))
  } else{
    return(RawMetadata)
  }
}

get.sample.attributes.INSDC <- function(sampleID=NA, apiKey=NA, BioPrjct=NA){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description get.sample.attributes.INSDC downloads all sample attributes (that is, additional environmental or other associated data) from INSDC. This is an expanded version of the get.BioProject.metadata.INSDC function and will cover all possible metadata variables found on INSDC, but requires a user-specified API-key to acces the INSDC databases. see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ to generate an API-key
  #' 
  #' @param sampleID a list. a list of one ore more SRA sample IDs (that is "Run" numbers). This argument can be left blank if input is provided for the BioPrjct argument (see further)
  #' @param apiKey a character string. A personal API-key to the acces the NCBI databases, and required to use the Entrez Programming Utilities (E-utilities). An API-key (API stands for application programming interface) is a unique identifier used to authenticate a user. A personal API-key can easily be generated at requested at https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
  #' @param BioPrjct a character string or a vector with character strings. Providing the associated BioProject numbers helps to give more understandable error messages. Alternativey, if the sampleID argument is empty, all sample IDs from the given BioProject numbers will be given to the sampleID argument.
  #'
  #' @details Each sequence data sample ("Run") typically has additional measurements or metadata associated
  #' with it. However, these are difficult to find, and cannot be downloaded together with the nucleotide 
  #' sequence data. The get.sample.attributes.INSDC function will fetch the all the metadata of the samples given to the
  #' sampleID (or BioPrjct) argument from the INSDC. 
  #' 
  #' To do this, get.sample.attributes.INSDC will use the Entrez Programming Utilities (E-utilities) from NCBI,
  #' and will access the databases with the API-key provided by the user. 
  #' see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ for more info
  #' on getting an API-key.
  #' 
  #' This function is only valid for samples and BioProjects that are registered on NCBI's SRA database.
  #' requires the packages xml2 and reshape2
  #' 
  #' @return a dataframe with the data found
  #' 
  #' @example m <- get.sample.attributes.INSDC(sampleID=c("ERR2226417", "ERR2226463"), apiKey=user_specified_apiKey)
  #' metadata_PRJNA369175 <- get.sample.attributes.INSDC(apiKey=user_specified_apiKey, BioPrjct=c("PRJNA369175"))
  
  require(xml2) #to parse the xml file
  require(reshape2) # long to wide format
  
  if(is.na(apiKey)){stop("No apiKey provided. 
                         A personal API key for the Entrez Programming Utilities (E-utilities) should be requested at the NCBI website.")
  }
  
  # check if param sampleID is present. if not, check for BioPrjct
  if(class(sampleID)!="character" | length(c(sampleID))==0 | 
     c(NA, NULL) %in% sampleID | c("") %in% sampleID){
    if(class(BioPrjct)!="character" | length(c(BioPrjct))==0 | 
       c(NA, NULL) %in% BioPrjct | c("") %in% BioPrjct){
      stop("No sampleID or BioPrjct provided, 
           or one of the sampleIDs/BioPrjct was NA, NULL or \"\"")
    }else{
      sampleID<-c()
      for(BP in BioPrjct){
        sampleID <- c(sampleID, get.BioProject.metadata.INSDC(BP, just.names=TRUE))
      }
    }
  }
  
  # see how many requests (sample ID's) will need to be handled
  # SRA can't handle large requests (lets say >50 at a time), so the IDs will be split into chunks of no mare than 50
  sampleID_setlist <- split(sampleID, ceiling(seq_along(sampleID)/50))
  env_metadata_full <- data.frame()
  for(sampleID_set in sampleID_setlist){
    sampleID_setSring <- paste(c(unlist(sampleID_set)), collapse="+")
    request_url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=",
                         sampleID_setSring, "&retmax=99999&api_key=", apiKey, sep="")
    metaXML <- tryCatch(read_xml(request_url[1]),
                        error = function(x){return(NULL)}
    )
    if(is.null(metaXML)){
      warning("No additional data found for one or more samples: better check this out.")
      env_metadata <- NULL
      next
    } else{
      # parse the xml file
      recs <- xml_find_all(metaXML, "//EXPERIMENT_PACKAGE/SAMPLE")
      attr_sample <- trimws(xml_attr(recs, "accession"))
      attr_name <- xml_find_all(metaXML, "//EXPERIMENT_PACKAGE/RUN_SET/RUN")
      attr_name <- trimws(xml_attr(attr_name, "accession"))
      attr_vals <- xml_find_all(metaXML, "//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/*[self::VALUE]")
      attr_vals <- trimws(xml_text(attr_vals))
      attr_cols <- xml_find_all(metaXML, "//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/*[self::TAG]")
      attr_cols <- xml_text(attr_cols)
      
      attr_BioPrj <- xml_find_all(metaXML, "//STUDY/IDENTIFIERS/EXTERNAL_ID")
      attr_BioPrj <- xml_text(attr_BioPrj)
      
      n_samples <- as_list(xml_find_all(metaXML, "//SAMPLE_ATTRIBUTES"))
      n_samples <- lapply(n_samples, function(x){length(x)})
      n_samples <- unlist(n_samples)
      namesVec <- c()
      for(ni in 1:length(attr_name)){
        namesVec <- c(namesVec, rep(attr_name[ni], n_samples[ni]))
      }
      env_metadata <- data.frame(attr_name=namesVec, 
                                 attr_col=attr_cols, value=attr_vals)
      
      env_metadata <- reshape2::dcast(env_metadata, attr_name ~ attr_cols, value.var="value")
      env_metadata$BioProject <- attr_BioPrj
      env_metadata$SRA_sample <- attr_sample
    }
    if(nrow(env_metadata_full)==0){
      env_metadata_full <- env_metadata
    }else{
      env_metadata_full <- combine.data(env_metadata_full, env_metadata, 
                                              fill=NA, variables.as.cols=TRUE)
    }
  }
  row.names(env_metadata_full)<-env_metadata_full$attr_name
  return(env_metadata_full)
  }

download.sequences.INSDC <- function(BioPrj = c(), destination.path = NA, apiKey=NA, 
                                     unzip = FALSE, keep.metadata = TRUE, download.sequences = TRUE){
  
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description download.sequences.INSDC downloads high throughput nucleotide sequence datasets that are deposited on INSDC, also taking care of any possible metadata and environmental data point of entry to INSDC is the SRA database from NCBI. 
  #' 
  #' @usage download.sequences.INSDC(BioPrj, destination.path, apiKey, unzip = FALSE, keep.metadata = TRUE, download.sequences = TRUE)
  #' 
  #' @param BioPrj a list with character strings. A list of one or more BioProject numbers to be downloaded. Required argument.
  #' @param destination.path a character string. The path to the directory where all the downloaded sequence data needs to go
  #' @param apiKey a character string. Only required if download.sequences.INSDC(keep.metadata=TRUE). A personal API-key to the acces the NCBI databases, and required to use the Entrez Programming Utilities (E-utilities). An API-key (API stands for application programming interface) is a unique identifier used to authenticate a user. You can easily generate an API-key: see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
  #' @param unzip boolean. If TRUE, the all *.fastq.gz files in the destination.path will unzipped. Default FALSE
  #' @param keep.metadata boolean. If TRUE, the downloaded metadata can be saved to a file (Console), if FALSE it is discarded. Default TRUE
  #' @param download.sequences boolean. If TRUE, the sequences will be downloaded to the destination.path. If FALSE, no sequences are downloaded. Default TRUE
  #' 
  #' @details download.sequences.INSDC will write the sequence data (*.fastq.gz files) to a destiation path, 
  #' the metadata (if it should be kept) is written to the Console and should be caught in an R-varaiable that
  #' can be later written to a csv file by the user.
  #' @return the sequence data and/or the metadata are written to the destination.path
  #' @example metadataFile <- download.sequences.INSDC(BioPrj = c("PRJEB1602", "PRJNA369175"), destination.path = getwd(),
  #' apiKey=user_specified_apiKey, unzip = FALSE, keep.metadata = TRUE, download.sequences = TRUE)
  
  ### 1. some initial checks before continuing
  
  # 1.1. check destination.path 
  if(is.na(destination.path)){destination.path <- getwd()} 
  # 1.2. check BioPrj
  if(length(BioPrj) == 0){ 
    stop("No BioProject number was given (required). 
         Valid input: a list of one or more BioProject numbers")
  }
  if(keep.metadata){
    metadata_all<-data.frame()
    # 1.3. check apiKey (only needed if keep.metadata is TRUE)
    if(is.na(apiKey)){stop("No apiKey provided. 
                                         A personal API key for the Entrez Programming Utilities (E-utilities) should be requested at the NCBI website.\n")
    }
    cat("Notice!\nthe metadata will be retruned to the Console\nIf you did assign the output of this function to an R-object (using \"<-\"): better abort and restart now\n\n")
  }
  
  # 2. start going though the BioProjects
  cat("Getting the metadata ...\n")
  downloads_failed<-0
  for(BP in BioPrj){
    cat(paste("Processing BioProject ",BP,"...\n", sep=""))
    faultyBioPrj<-FALSE
    ### 2.1. Downloading the metadata from SRA to get the run numbers (=samples) of the BioProject
    ###    These run numbers are the required input to download the sequence data, as sequence data cannot be redectly downloaded via the BioProject number 
    faultyBioPrj <- tryCatch({
      RawMetadata <- get.BioProject.metadata.INSDC(BP)
      cat(paste("\t",as.character(nrow(RawMetadata))," samples (Runs)...\n", sep=""))
    }, error = function(e) {
      return(TRUE)
    })
    
    if(! is.null(faultyBioPrj) && faultyBioPrj){
      downloads_failed <- downloads_failed+1
      cat(paste("\t No data found for ", BP,". Could not download the data.\n", sep=""))
    } else{
      faultyBioPrj<-FALSE
    }
    

    
    
    ### 2.2 add the FTP path to download the sequences
    ftp_url <- paste("https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=", BP, 
                     "&result=read_run&fields=run_accession,fastq_ftp&download=txt", sep="")
    tmpFile <- tempfile()
    download.file(ftp_url, destfile = tmpFile, method="auto", quiet = TRUE)
    ftps <- read.table(tmpFile, header=TRUE)
    ftps <- ftps[ftps$run_accession[order(RawMetadata$Run)],]
    RawMetadata$ftp <- ftps$fastq_ftp
    file.remove(tmpFile)
    
    ### 2.2. if keep.metadata==TRUE, the raw metadata must be cleaned up, and if applicable added to the other metadata
    if(keep.metadata & !faultyBioPrj){
      cat("\tprocessing the metadata ...\n")
      sample_numbers<-as.character(RawMetadata$Run)
      env_data <- get.sample.attributes.INSDC(sampleID = sample_numbers, 
                                              apiKey = apiKey, BioPrjct = BP)
      rownames(RawMetadata) <- RawMetadata$Run
      RawMetadata <- RawMetadata[order(row.names(RawMetadata)),]
      rownames(env_data) <- env_data$attr_name
      env_data <- env_data[order(row.names(env_data)),]
      MetaData <- cbind(RawMetadata, env_data)
      
      if(nrow(metadata_all)>0){
        metadata_all <- process.metadata(metadata=MetaData, add_to = metadata_all, strict.MIxS = FALSE, ask.input = FALSE)
      } else{
        metadata_all <- process.metadata(metadata=MetaData, strict.MIxS = FALSE, ask.input = FALSE)
      }
    }
    
    ### 2.3. Downloading the sequence data from SRA
    #      This requires a sqlite database of what is on SRA. This database just needs to be downloaded once, but can be very big.
    #      No way of getting te data without the SRA sqlite database...
    #      Also using the sample names gathered in (1) instead of the BioProject number
    if(download.sequences & !faultyBioPrj){
      cat("\tDownloading the sequence data ...\n")
      ftps <- as.character(RawMetadata$ftp)
      ftps <- sapply(ftps, function(x){strsplit(x, ";")})
      ftps <- unname(unlist(ftps))
      for(ftp_url in ftps){
        destfile <- strsplit(ftp_url, "/")[[1]]
        destfile <- destfile[length(destfile)]
        download.file(ftp_url, destfile = file.path(destination.path, destfile),method="auto")
      }
    }else{
      downloads_failed <- downloads_failed+1
    }
  }
  
  ### 2.4. unzipping the data, if required by the user
  if(unzip == TRUE & !faultyBioPrj){
    cat("Unzipping files ...\n")
    fileList<- list.files(path=getwd(), pattern='*.fastq.gz')
    for(file in fileList){
      fileName= strsplit(file, ".gz")[[1]]
      system2(command="gunzip", args=paste(" -c ", file," > ", fileName, sep=''), stdout=TRUE)
      system2(command="rm", args=file)
    }
  }
  
  ### 3. finalizing
  if(downloads_failed < length(BioPrj)){
    outputmessage <- paste("Finished processing\nThe files are in ",destination.path,"\n", sep="")
    if(keep.metadata){
      outputmessage <- paste(outputmessage, "Please consider the warning messages for a manual quality controll of the metadata\nand environmental data\n")
    }
  } else{
    outputmessage <- paste("Finished processing\nNo files were downloaded\n", sep="")
  }
  cat(outputmessage)
  if(keep.metadata){return(metadata_all)}
}


#--------------------------------------------------------------
# 3. Tools for submitting data
#--------------------------------------------------------------
write.metadata.MIxS.as.ENA <- function(metadata, name=NULL,
                                 unique_name_prefix=NA, checklist_accession=NA,
                                 tax_name=NA, ask.input=TRUE,
                                 missing.data.as.empty.columns=TRUE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description write.metadata.MIxS.as.ENA converts metadata into the right format (a tab separated text file) to submit data to ENA (version Novembre 2019)
  #' @param metadata a metadata.MIxS class object. The object to be written as text file suited for submission to ENA.
  #' @param name a character string naming a file. Either a full path or a file name to write to the working directory. Invalid input will cause output to be printed to the console.
  #' @param checklist_accession character. The name of a MIxS environmetal package or it's ENA checklist accession number.
  #' @details ENA uses its own variant of MIxS, and has
  #' write.metadata.MIxS.as.ENA assumes the dataset has already been subjected to process the process.metadata function to standardize in input. If this is not the case, "garbage in, garbage out" is applicable. 
  #' @return tab separated .txt file
  #' @example 
  #' 

  warningmessages<-c()
  
  # 1. check input
  if(check.valid.metadata.MIxS(metadata)){
    temp_unique_name_prefix <- deparse(substitute(metadata))
    metaunits <- metadata@units
    env_package <- metadata@env_package
    metadata <- metadata@data
  }else{
    stop("metadata argument must be a metadata.MIxS class object.
         Common dataframes should be converted to metadata.MIxS with 
         the data.frame.to.metadata.MIxS function.")
  } 
  # 2. check output destination
  if(!is.character(name) | c(NULL,NA) %in% name | length(name)>1){
    stop("No correct file name or file path provided in the name argument.")
  } else{
    if(grepl("/", name)){
      file.name <- strsplit(name, "/")[[1]]
      out.dir <- paste(file.name[1:(length(file.name)-1)], collapse="/")
      file.name <- file.name[length(file.name)]
      if(!dir.exists(file.path(out.dir))){
        stop("Could not find the directory provided in the name argument.")
      }
    } else{
      file.name <- paste(name, ".txt", sep="")
      out.dir <- getwd()
    }
    full.name <- paste(out.dir, file.name, sep="/")
  }
  
  # 3. checklist_accession
  if(c(NULL,NA) %in% checklist_accession){
    checklist_accession <- env_package
    if(checklist_accession %in% c("not_specified", "multiple_packages")){
      checklist_accession <- ENA_checklistAccession[ENA_checklistAccession$env_package=="other_not_MIxS",]$ena_checklist_accession
      warningmessages <- multi.warnings("The environmental package of the input data was not specified or refered to multiple packages. The ENA default was therfore selected. This can be overridden by providing input for the checklist_accession argument.",
                                        warningmessages)
    } else{
      checklist_accession <-  as.character(ENA_checklistAccession[ENA_checklistAccession$env_package==checklist_accession,]$ena_checklist_accession)
    }
  }else if(length(checklist_accession)==1){
    if(checklist_accession %in% ENA_checklistAccession$env_package){
      checklist_accession <- as.character(ENA_checklistAccession[ENA_checklistAccession$env_package==checklist_accession,]$ena_checklist_accession)
    }else if(! checklist_accession %in% ENA_checklistAccession$ena_checklist_accession){
      stop("The input for the checklist_accession argument was invalid. 
           Must be a MIxS package (check spelling, use of underscores) or
           an ENA checklist accession number that refers to a MIxS package.")
    }
  } else{
    stop("The input for the checklist_accession argument was invalid.")
  }
  out.line01 <- paste("#checklist_accession", checklist_accession, sep='\t')
  
  # 4. unique_name_prefix
  if(c(NULL,NA) %in% unique_name_prefix | length(unique_name_prefix)>1){
    if(ask.input){
      cat("Please provide a unique name prefix (unique_name_prefix argument)\n") 
      unique_name_prefix <- readline() 
    }else{
      unique_name_prefix <- temp_unique_name_prefix
    }
  }
  out.line02 <- paste("#unique_name_prefix", unique_name_prefix, sep='\t')

  # 5. the actual data:
  ena_metadata <- data.frame(row.names=rownames(metadata)) #the values
  ena_variable <- c() #the row with the variable names (use of spaces by ENA makes is difficult to include with the values)
  ena_units <- c() #the unit row
  # 5.1 fixed term sample_alias = unique_name_prefix + a running number
  ln <- nchar(trunc(nrow(metadata)))
  lx <- c(1:nrow(metadata))
  sample_alias <- paste(unique_name_prefix, formatC(lx, width=ln, flag="0"), sep="_")
  ena_metadata$sample_alias<-sample_alias
  ena_variable <- c(ena_variable, "#sample_alias")
  ena_units <- c(ena_units, "#units")
  
  # 5.2 fixed term taxon
  if(c(NULL,NA) %in% tax_name){ #see if input was provided
    #look for a clue in the data
    c_taxa<-intersect(c("subspecf_gen_lin", "scientific_name"), colnames(metadata))
    if(length(c_taxa)>1){
      tax_name <- metadata[,colnames(metadata)==c_taxa[1]]
    }else{
      cat("Please specify what taxon was targeted in this dataset (e.g. Bacteria, Cyanobacteria, Synechococcus,...). Type n to ignore.\n")
      tax_name <- readline() 
      if(!tax_name %in% c("n", "N")){
        tax_name <- rep(tax_name, nrow(metadata))
      } else{
        tax_name <- rep("", nrow(metadata))
      }
    }
  }else{
    tax_name <- rep(tax_name, nrow(metadata))
  }
    
  # 5.3 fixed term taxID
  ena_metadata$tax_IDs <- sapply(tax_name, FUN=function(x){commonTax.to.NCBI.TaxID(x, fill.unknown="")})
  ena_variable <- c(ena_variable, "tax_id")
  ena_units <- c(ena_units, "")
  ena_metadata$tax_name <- tax_name
  ena_variable <- c(ena_variable, "scientific name")
  ena_units <- c(ena_units, "")
  
  # 5.4 fixed term sample_title	(original sample name)
  if("original_name" %in% colnames(metadata)){
    ena_metadata$sample_title<-metadata$original_name
  }else if("misc_param" %in% colnames(metadata)){
    msc <- strsplit(metadata[1,]$misc_param, ";")
    msc <- unlist(lapply(msc, function(x){strsplit(x, ":")}))
    if("original_name" %in% msc){
      sample_title <-c()
      for(stl in 1:nrow(metadata)){
        sName <- strsplit(metadata[stl,]$misc_param, ";")[[1]]
        sName <- sName[grepl("original_name:", sName)]
        sName <- strsplit(sName, ":")[[1]][2]
        if(sName=="NA" | sName=="NULL"){
          sName <- ""
        }
        sample_title <- c(sample_title, sName)
      }
      ena_metadata$sample_title <- sample_title
    }
  }else{
    ena_metadata$sample_title<-rownames(metadata)
  }
  ena_variable <- c(ena_variable, "sample_title")
  ena_units <- c(ena_units, "")
  
  # 5.5 fixed term sample_description
  if("sample_description" %in% colnames(metadata)){
    ena_metadata$sample_description<-metadata$sample_description
  }else{
    if(ask.input){
      cat("No sample description was found. Could you provide an overall sample description(y/n)?\n") 
      descr <- readline() 
      if(descr %in% c("y", "Y", "yes", "YES", "Yes")){
        cat("Please type the sample description that will be used for all samples in the dataset\n") 
        descr <- readline() 
        ena_metadata$sample_description<-rep(descr, nrow(metadata))
      }else if(descr %in% c("n", "N", "no", "NO", "No")){
        ena_metadata$sample_description<-rep("", nrow(metadata))
      }else{
        stop("incorrect input... only yes or no allowed.")
      }
    }
  }
  ena_variable <- c(ena_variable, "sample_description")
  ena_units <- c(ena_units, "")
  # 5.6 latlon DD
  if("DwC_decimalLatitude" %in% colnames(metadata) & 
     "DwC_decimalLongitude" %in% colnames(metadata)){
    ena_metadata$DwC_decimalLatitude <- metadata$DwC_decimalLatitude
    ena_metadata$DwC_decimalLongitude <- metadata$DwC_decimalLongitude
    ena_variable <- c(ena_variable, "geographic location (latitude)", "geographic location (longitude)")
    ena_units <- c(ena_units, "DD", "DD")
  }else if("lat_lon" %in% colnames(metadata)){
    metadata <- separate(metadata, "lat_lon", sep=" ", into=c("DwC_decimalLatitude", "DwC_decimalLongitude"))
    ena_metadata$DwC_decimalLatitude <- metadata$DwC_decimalLatitude
    ena_metadata$DwC_decimalLongitude <- metadata$DwC_decimalLongitude
    ena_variable <- c(ena_variable, "geographic location (latitude)", "geographic location (longitude)")
    ena_units <- c(ena_units, "DD", "DD")
  }else{
    warningmessages <- multi.warnings("The data lacks geographical coordinates.", warningmessages)
  }
  
  # 5.7 convert remaining MIxS terms to ENA variants
  redundant_cols <- which(colnames(metadata) %in% c("sample_description", "lat_lon", "DwC_decimalLatitude", 
                                                   "DwC_decimalLongitude", "original_name", "subspecf_gen_lin", 
                                                   "scientific_name"))
  metadata <- metadata[,-redundant_cols]
  metaunits <- metaunits[-redundant_cols]
  for(cl in 1:ncol(metadata)){
    clName <- colnames(metadata[cl])
    ena_name <- as.character(TermsLib[TermsLib$name==clName,]$name_variant_ENA)
    if(length(ena_name) == 0 || ena_name==""){
      ena_name <- clName
    }
    ena_metadata[,clName] <- c(as.character(metadata[,colnames(metadata) %in% clName]))
    ena_variable <- c(ena_variable, ena_name)
    ena_units <- c(ena_units, gsub("alphanumeric", "", metaunits[cl]))
  }
  
  # 6. put lines 3 to 5 together 
  out.line03 <- paste(ena_variable, collapse='\t')
  out.line04 <- paste("#template", paste(c(as.character(ena_metadata[1,][-1])), collapse="\t"), sep='\t')
  out.line05 <- paste(ena_units, collapse='\t')
  
  # 7. put line >5 together
  out.line99 <- apply(ena_metadata, 1, function(x) paste(x, collapse = "\t")) 
  out.line99 <- paste(out.line99, collapse="\n")
  
  # 8. combine everything
  output <- paste(out.line01, out.line02, out.line03, out.line04, out.line05, out.line99, sep="\n")
  
  cat(paste("The data had been written to ", out.dir, "/", file.name, "\n",sep=""))
  write.table(output, file=paste(out.dir, "/", file.name, sep=""), 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}







#--------------------------------------------------------------
# test zone
#--------------------------------------------------------------

Heind <- read.csv("/Users/msweetlove/Desktop/historic_fish/MiMARKS_Heindler_PRJEB34858.csv", row.names=1)

HeindQC <- process.metadata(metadata = Heind, strict.MIxS = FALSE, 
                            out.format="metadata.MIxS", ask.input=TRUE)

write.metadata.MIxS.as.ENA(metadata=HeindQC, name="/Users/msweetlove/Desktop/historic_fish/MiMARKS_Heindler_PRJEB34858.txt",
                           unique_name_prefix="HistoricAntarcticFishDataset_2019_", checklist_accession=NA,
                           tax_name="Bacteria", ask.input=TRUE,
                           missing.data.as.empty.columns=TRUE)






testdir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/test"


wdir2="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir"
setwd(wdir2)

#test BioProject
#BP<-"PRJNA369175"
#dataat polaRRTools::MIxSTerms.RData
#apiKey<-"5e490715dbb88b3f861565b05aab426f1408"
#test<-process.metadata(metadata=MetaData, strict.MIxS = FALSE)

#get.BioProject.metadata.INSDC(BioPrjct="PRJNA369175", just.names=TRUE)




test1<- get.sample.attributes.INSDC(sampleID=NA, apiKey="5e490715dbb88b3f861565b05aab426f1408",
                                   BioPrjct=c("PRJEB27415", "PRJNA369175"))
test1<- get.sample.attributes.INSDC(sampleID=NA, apiKey="5e490715dbb88b3f861565b05aab426f1408",
                                    BioPrjct=c("PRJNA395930"))

test2 <- process.metadata(metadata = test1, add_to = NA, strict.MIxS = FALSE, 
                             out.format="metadata.MIxS", ask.input=TRUE)
 

test3<- get.sample.attributes.INSDC(sampleID=NA, apiKey="5e490715dbb88b3f861565b05aab426f1408",
                                    BioPrjct=c("PRJEB1602"))

test4 <- process.metadata(metadata = test3, add_to = test2, 
                          strict.MIxS = FALSE)

test5 <- download.sequences.INSDC(BioPrj = c("PRJEB1602", "PRJNA369175"), 
                                  destination.path = testdir,
                                  apiKey="5e490715dbb88b3f861565b05aab426f1408",
                                              unzip = FALSE, 
                                              keep.metadata = TRUE,
                                              download.sequences = FALSE)


test5 <- download.sequences.INSDC(BioPrj = c("PRJEB1602"), 
                                  destination.path = testdir,
                                  apiKey="5e490715dbb88b3f861565b05aab426f1408",
                                  unzip = FALSE, 
                                  keep.metadata = TRUE,
                                  download.sequences = FALSE)

test5 <- download.sequences.INSDC(BioPrj = c("PRJ55555555"), 
                                  destination.path = testdir,
                                  apiKey="5e490715dbb88b3f861565b05aab426f1408",
                                  unzip = FALSE, 
                                  keep.metadata = TRUE,
                                  download.sequences = FALSE)


test<-get.BioProject.metadata.INSDC(BioPrjct="PRJNA369175", just.names=FALSE)
test<-get.BioProject.metadata.INSDC(BioPrjct="PRJNA395930", just.names=FALSE)
test<-get.BioProject.metadata.INSDC(BioPrjct="PRJEB23732", just.names=FALSE)
test<-get.BioProject.metadata.INSDC(BioPrjct=4, just.names=FALSE)

test<-get.BioProject.metadata.INSDC(BioPrjct="PRJNB23732", just.names=FALSE)#shuld be error




