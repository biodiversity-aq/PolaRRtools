### check out package down
#pkgdown::build site run loally to build site

#structure cookiecutter
#test code testy/Travis IC

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
# file encdong UTF-8
#
#   This is version includes utilities to:
#     - download sequence data from INSDC (NCBI's SRA to be more precise)
#     - download meta data and environmental data fron INSDC (NCBI's SRA to be more precise)
#     - Quality controll meta data and environmental data
# 
#==============================================================


#library(devtools)
#library(usethis)
TermsLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/TermsLibrary_v2.csv", 
                      header=TRUE)
TermsSyn <- sapply(as.character(TermsLib$synonyms), function(x){strsplit(x, ";")})
names(TermsSyn) <- TermsLib$name 
#usethis::use_data(MIxS_Voc, internal = FALSE)


ENA_checklistAccession<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/ENA_checklistAccession.csv", 
                   header=TRUE)

ENA_allowed_terms<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/ENA_allowed_terms.csv", 
                                 header=TRUE, stringsAsFactors = FALSE)
ENA_geoloc<-ENA_allowed_terms$geo_loc_name
ENA_geoloc<-ENA_geoloc[!ENA_geoloc==""]
ENA_instrument<-ENA_allowed_terms$instrument_model
ENA_instrument<-ENA_instrument[!ENA_instrument==""]
ENA_select<-ENA_allowed_terms$library_selection
ENA_select<-ENA_select[!ENA_select==""]
ENA_strat<-ENA_allowed_terms$library_strategy
ENA_strat<-ENA_strat[!ENA_strat==""]


MarsLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/MarsLibrary.csv", 
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

setClass("metadata.DwC", slots=list(name="character", #the name of the resource
                                    core="data.frame", #a dataframe with the DwC core (events/occurences are rows, variables are columns)
                                    emof="data.frame", #a dataframe with the DwC EMOF extention (=the environmental data)
                                    type="character", #type of the core: occurence or event
                                    eml_url="character" #the url to the EML file from the IPT
)
)

setMethod("show",
          "metadata.DwC",
          function(object) {
            N <- nrow(object@core)
            if(nrow(object@emof)>0){
              E <- "and an EMOF extension"
            } else{E <- ""}
            cat(paste("a metadata.DwC class object with ", as.character(N), " ",
                      object@type, "s", as.character(E), ".\n", sep=""))
          }
)

check.valid.metadata.DwC <- function(d){
  valid <- TRUE
  #class must be metadata.MIxS 
  if(!class(d)=="metadata.DwC"){
    valid <- FALSE
  }else{
    dcol <- ncol(d@core)
    drow <- nrow(d@core)
    #The core must have at least 1 sample (row) and 1 variable (column)
    if(dcol<=0 & drow<=0){
      valid <- FALSE
    }
    #type must be occurence or event
    if(d@type="occurence"){
      if(!"occurenceID" %in% colnames(d@core)){
        valid <- FALSE
      }
    }else if(d@type="event"){
      if(!"eventID" %in% colnames(d@core)){
        valid <- FALSE
      }
    }else{
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
  
  # 0. pre-process input
  # 0.1. check input data
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
        
  # 0.2 formatting
  # remove empty columns
  metadata <- metadata[,colSums(is.na(metadata)) < nrow(metadata)] 
  # clean up columnames
  colnames(metadata) <- gsub("[\\.]+", "_", colnames(metadata)) # replace double dots with underscore
  colnames(metadata) <- gsub("_$", "", colnames(metadata))  # remove trailing underscore
  colnames(metadata) <- tolower(colnames(metadata)) # all to lowercase
  
  # 0.3 check if data is in one-header table or if there are additional MiMARKS header lines
  # for additional MIxS headers: "environmental package", "units template" => only units of importance
  pre_def_units <- FALSE
  if(grepl("units", tolower(row.names(metadata)[1]))){
    units <- data.frame(var_name=c(colnames(metadata)), unit= unlist(metadata[1,]))
    metadata <- data.frame(metadata[-1,], stringsAsFactors = FALSE)
    pre_def_units <- TRUE
  } else if(grepl("units", tolower(row.names(metadata)[2]))){
    units <- data.frame(var_name=c(colnames(metadata)), unit= unlist(metadata[2,]))
    metadata <- data.frame(metadata[-c(1,2),], stringsAsFactors = FALSE)
    pre_def_units <- TRUE
  }

  
  # 0.4 prepare output data:
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
    if(.row_names_info(metadata)>0){#rownames provided by the user might be the original names
      if(ask.input){
        cat("No original sample names found...\n\tuse the rownames instead? (y/n)\n") 
        ctu <- readline() 
        if(ctu %in% c("y", "Y", "yes", "YES", "Yes")){
          New_metadata$original_name <- row.names(metadata)
          warningmessages <- multi.warnings("no original sample names found, used the rownames instead", warningmessages)
        }else if(ctu %in% c("n", "N", "no", "NO", "No")){
          cat("you chose no.\n\t Can you give the columnname with the original sample names instead? Type n to ignore\n") 
          ctu2 <- readline() 
          if(ctu2 %in% colnames(metadata)){
            New_metadata$original_name <- metadata[,ctu3]
          }else if(!ctu2 %in% c("n", "N", "no", "NO", "No")){
            stop("incorrect input... only yes or no allowed.")
          }else{
            stop(paste("could not find the column", ctu2, "in the colnames of the dataset provided..."))
          }
        } else{
          stop("incorrect input... only yes or no allowed.")
        }
      }
     }else{
      warningmessages <- multi.warnings("no original sample names found", warningmessages)
    }
  }
  
  # 2. some basic info from insdc
  TermsSyn_insdc<-TermsSyn[as.character(TermsLib[TermsLib$name_origin=="INSDC",]$name)]
  for(item in names(TermsSyn_insdc)){
    item_shared <- intersect(TermsSyn_insdc[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    }
  }
  
  # 3. dealing with latitude-longitude, and it's many possible formats...
  TermsSyn_latlon<-TermsSyn[as.character(TermsLib[TermsLib$name=="lat_lon",]$name)]
  TermsSyn_lat<-TermsSyn[as.character(TermsLib[TermsLib$name=="decimalLatitude",]$name)]
  TermsSyn_lon<-TermsSyn[as.character(TermsLib[TermsLib$name=="decimalLongitude",]$name)]
  
  QClatlon <- dataQC.LatitudeLongitudeCheck(metadata, 
                                            latlon.colnames=list(TermsSyn_latlon[[1]],
                                                                 TermsSyn_lat[[1]],
                                                                 TermsSyn_lon[[1]]))
  warningmessages<-c(QClatlon$warningmessages, warningmessages)
  New_metadata$lat_lon <- QClatlon$values
  New_metadata$decimalLatitude <- sapply(New_metadata$lat_lon, function(x){strsplit(x, " ")[[1]][1]})
  New_metadata$decimalLongitude <- sapply(New_metadata$lat_lon, function(x){strsplit(x, " ")[[1]][2]})
  
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
  if(!"env_package" %in% colnames(metadata)){
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
    New_metadata$env_package <- metadata$env_package
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
        env_package <- dataQC.guess.env_package.from.data(metadata)
        warningmessages<-c(warningmessages, env_package$warningmessages)
        env_package <- env_package$values
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
        cat("Found none or more than one environmental packages. Only one allowed when strict.MIxS is TRUE\n\tContinue by setting strict.MIxS to FALSE? (y/n)\n") 
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
        cat("Some of the required terms for strict.MIxS are missing...\n\tContinue with strict.MIxS turned off? (y/n)\n") 
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
  # 7.1 already registered terms
  TermsSyn_add <- TermsSyn[as.character(TermsLib[TermsLib$name_origin %in% c("DwC", "miscellaneous"),]$name)]
  TermsSyn_add <- TermsSyn_add[!names(TermsSyn_add) %in% c("decimalLatitude", "decimalLongitude")]
  for(item in names(TermsSyn_add)){
    item_shared <- intersect(TermsSyn_add[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    }
  }
  # 7.2 novel terms
  unknown_terms <- setdiff(colnames(metadata), unlist(TermsSyn))
  if(length(unknown_terms) > 0 & strict.MIxS & ask.input){
    cat("non-MIxS variables were encountere.\n\tContinue with strict.MIxS turned off? (y/n)\n") 
    ctu <- readline() 
    if(ctu %in% c("y", "Y", "yes", "YES", "Yes")){
      strict.MIxS<-FALSE
    }else if(ctu %in% c("n", "N", "no", "NO", "No")){
      stop("Could not output metadata strictly following the MIxS rules.")
    } else{
      stop("incorrect input... only yes or no allowed.")
    }
  }
  if(!strict.MIxS){
    if(length(unknown_terms) > 5){
      # if there are too much novel terms, don't go over all of them
      if(ask.input){
        t<-paste(unknown_terms, collapse=", ")
        cat(paste("The following unknown variables were encountered:\n",t," \n\tAdd all to the QC'd data? (y/n)\n", sep=" ")) 
        ctu <- readline() 
        if(ctu %in% c("y", "Y", "yes", "YES", "Yes")){
          for(t in unknown_terms){
            New_metadata[,t] <- metadata[,t]
          }
        }
      }else{
        warningmessages<-multi.warnings("Some unknown variables present in the data", warningmessages)
        for(t in unknown_terms){
          New_metadata[,t] <- metadata[,t]
        }
      }
    }else{
      for(t in unknown_terms){
        if(ask.input){
          cat(paste("The following unknown variable was encountered:",t," \n\tAdd to the QC'd data? (y/n)\n", sep=" ")) 
          ctu <- readline() 
          if(ctu %in% c("y", "Y", "yes", "YES", "Yes")){
            New_metadata[,t] <- metadata[,t]
          }
        }else{
          New_metadata[,t] <- metadata[,t]
          warningmessages<-multi.warnings("Some unknown variables present in the data", warningmessages)
        }
      }
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
      cat("No investigation_type was found...\n\tPlease provide an investigation_type. Common ones include mimarks-survey or metagenome. Type n to ignore.\n") 
      invtype <- readline() 
      if(! invtype %in% c("n", "N")){
        New_metadata$investigation_type <- rep(invtype, nrow(New_metadata))
      }
    }
  }
  if("target_gene" %in% colnames(New_metadata)){
    New_metadata$target_gene <- sapply(New_metadata$target_gene, function(x){
      if(grepl("16S", x)){x<-"16S ssu rRNA"
      }else if(grepl("18S", toupper(x))){x<-"18S ssu rRNA"
      }else if(grepl("ITS", toupper(x))){x<-"ITS"
      }else if(grepl("COI", toupper(x))){x<-"COI"
      }
      return(x)})
  }
  if("specific_host" %in% colnames(New_metadata)){
    host_val <- setdiff(unique(New_metadata$specific_host), c(NA, "NA", "-", "not applicable"))
    if(identical(host_val, character(0))){
      New_metadata <- New_metadata[,!colnames(New_metadata) %in% "specific_host"]
    }
  }
  #if("eventID" %in% colnames(New_metadata) & ask.input){
  #  if(length(unique(New_metadata$eventID))==nrow(New_metadata)){
  #    cat("There are the same number of events as there are samples...\n\tkeep $eventID? (y/n)\n") 
  #    ctu <- readline() 
  #    if(ctu %in% c("n", "N", "no", "NO", "No")){
  #      New_metadata <- New_metadata[,!colnames(New_metadata) %in% "eventID"]
  #    }
  #  }
  #}
  
  # QC on length/depth/size measurements
  # if multiple units found: do nothing
  # if a string found that is not one of the expected units: leave it and work with the rest
  size_related <- c("depth", "elev", "alt_elev", "filter_size", "tot_depth_water_col")
  items_shared <- intersect(size_related, colnames(New_metadata))
  QC_units<-c()
  if(length(items_shared)>0){
    alternative_units<-data.frame(name=items_shared, unit_full=rep(NA, length(items_shared)))
    possible_units <- c("nm", "um", "μm", "mm", "cm", "dm", "m", "km")
    for(name_unit in items_shared){
      vals <- as.vector(as.character(New_metadata[,colnames(New_metadata) %in% name_unit]))
      names(vals)<-row.names(New_metadata)
      #extract the units
      if(pre_def_units){
        val_units <- as.character(units[units$var_name==name_unit,]$unit)
        val_units <- rep(val_units,nrow(New_metadata))
      }else{
        val_units <- unlist(lapply(vals, function(v){gsub("[0-9]|\\.|,|-", "", v)} ))
      }
      val_units <- gsub(" ", "", tolower(val_units))
      val_units <- gsub("nano", "n", tolower(val_units))
      val_units <- gsub("micro", "u", tolower(val_units))
      val_units <- gsub("mili", "m", tolower(val_units))
      val_units <- gsub("centi", "c", tolower(val_units))
      val_units <- gsub("deci", "d", tolower(val_units))
      val_units <- gsub("kilo", "k", tolower(val_units))
      val_units <- gsub("meter", "m", tolower(val_units))
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
  if(!pre_def_units){
    # there were no pre defined units in a additional header line
    for(i in 1:ncol(New_metadata)){
      New_metadata_units[colnames(New_metadata)[i]] <- as.character(TermsLib[TermsLib$name %in% colnames(New_metadata)[i],]$expected_unit)
    }
  }else{
    # there were pre defined units in a additional header line
    # user defined units obviously have priority over assumed standard units in the TermsLib file
    for(i in 1:ncol(New_metadata)){
      if(colnames(New_metadata)[i] %in% units$var_name){
        New_metadata_units[colnames(New_metadata)[i]] <- as.character(units[units$var_name==colnames(New_metadata)[i],]$unit)
      } else{
        New_metadata_units[colnames(New_metadata)[i]] <- as.character(TermsLib[TermsLib$name %in% colnames(New_metadata)[i],]$expected_unit)
      }
    }
  }
  # QC_units
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
    if(colnames(New_metadata)[i] %in%TermsLib$name){
      New_metadata_section[colnames(New_metadata)[i]] <- as.character(TermsLib[TermsLib$name %in% colnames(New_metadata)[i],]$MIxS_section)
    } else{
      New_metadata_section[colnames(New_metadata)[i]] <- "miscellaneous"
    }
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

semantics <- function(dataset = NA, vocab=c("MIxS", "MiMARKS", "DarwinCore")){
  # use the library to shore up the semantics of a dataset using a volcabulary (DwC or MIxS)
  # make sure all the terms are the correct terms
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

# idea: prep.metadata.ENA: output csv for editing, later convert cvs to tsv
# idea: check if sequence files correspond to the names in the metadata file

prep.metadata.ENA <- function(metadata, dest.dir=NULL, file.name=NULL,
                              sample.unique_name_prefix=NA, checklist_accession=NA,
                              tax_name=NA, ask.input=TRUE,
                              insert.size=NA, library.layout=NA,
                              library.strategy=NA, library.selection=NA,
                              seq.file.extension=".fastq.gz"
                              ){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description converts metadata into the right format (a tab separated text file) to submit data to ENA (version Januari 2020)
  #' @param metadata a metadata.MIxS class object. The object to be written as text file suited for submission to ENA.
  #' @param file.name a character string. A name (without a file type extension) to use for the output files.
  #' @param dest.dir a character string. The file path to the directory where the output files must be written. If left blank files are written to the working directory.
  #' @param sample.unique_name_prefix a character string. The unique name prefix to append to the sample names (to tie them all together). Required for ENA submissions
  #' @param checklist_accession a character string. The name of a MIxS environmetal package or it's ENA checklist accession number.
  #' @param tax_name a character string. The scientific name of a taxon targeted in the data (applicable to all the samples). Same as "subspecf_gen_lin" or "scientific_name" in the metadata.MIxS input.
  #' @param ask.input boolean. Whether or not to ask for user input to make decisions or solve problems that arise during the reformatting (e.g. missing data,...)

  #' @param insert_size a character string. The size of the reads (in number of basepairs), if applicable to all samples.
  #' @param library_layout a character string. The layout of the library, either PAIRED or SINGLE, if applicable to all samples.
  #' @param library_strategy a character string. The library strategy (e.g. AMPLICON, WGS,...), if applicable to all samples.
  #' @param library_selection a character string. The method used to select for, enrich or or screen the material being sequenced (e.g. PCR)?, if applicable to all samples.
  #' @param seq.file.extension a character string. The extension for the sequence files. Default is .fastq.gz
  #' 
  #' @details This function will reformat metadata for submission to ENA. Specifically made for ecological environmental studies (e.g. amplicon sequencing, shotgun metagenomics,...), with additional QCs build in for Antarctic and Southern Ocean data.
  #' prep.metadata.ENA assumes the dataset has already been subjected to process the process.metadata function to standardize and QC in input.
  #' @return a .tsv file with the sample metadata, and a _runInfo.tsv file with the technical data
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
  if(!is.character(file.name) | c(NULL,NA) %in% file.name | length(file.name)>1 | grepl("/", file.name)){
    stop("No correct file name provided in the file.name argument.")
  }else{
    if(grepl(".txt$", file.name)){
      file.name <- gsub(".txt$", "", file.name)
    }else if(grepl(".csv$", file.name)){
      file.name <- gsub(".csv$", "", file.name)
    }
  }
  
  if(!is.null(dest.dir)){
    if(!dir.exists(file.path(dest.dir))){
      stop("Could not find the directory provided in the name argument.")
    }
    dest.dir<-gsub("/$","", dest.dir)
  }else{
    dest.dir <- getwd()
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
  ena_variable <- c(ena_variable, "sample_alias")
  ena_units <- c(ena_units, "#units")
  
  # 5.2 fixed term taxon: common_name
  # = The taxonomic classification (based on the ENA taxIDs) of the sample
  # for molecular data: assume environmental metagnemome unless enough evidence geven for organism/taxon
  if(c(NULL,NA) %in% tax_name){ #see if input was provided
    #look for a clue in the data
    c_taxa<-intersect(c("subspecf_gen_lin", "scientific_name"), colnames(metadata))
    if(length(c_taxa)>=1){
      tax_name <- metadata[,colnames(metadata)==c_taxa[1]]
    }else{
      if(ask.input){
        cat("No Target taxon found:\n\tPlease type a taxon name (e.g. Bacteria, Cyanobacteria, Synechococcus,...). Or type n to ignore.\n")
        tax_name <- readline() 
        if(!tax_name %in% c("n", "N")){
          tax_name <- rep(tax_name, nrow(metadata))
        } else{
          tax_name <- rep("not collected", nrow(metadata))
        }
      }else{
        tax_name <- rep("not collected", nrow(metadata))
      }
    }
  }else{
    tax_name <- rep(tax_name, nrow(metadata))
  }
  # common tax name
  ena_metadata$common_name <- tax_name
  ena_variable <- c(ena_variable, "common_name")
  ena_units <- c(ena_units, "")
  
  # 5.3 fixed term taxID and scientific tax name
  # taxID
  # investigation_type has priority because if sample is environmental, 
  # a unknown variety and number of organisms will be present, thus one of the environmental 
  # metagenome taxIDs should be taken
  if("investigation_type" %in% colnames(metadata)){
    ena_metadata$tax_name <- sapply(metadata$investigation_type, function(x){
      x <- gsub("mimarks-survey", "metagenome", x)
      return(x)})
    if(!all(unique(ena_metadata$tax_name)=="metagenome")){
      for(itx in 1:length(ena_metadata$tax_name)){
        if(!ena_metadata[itx,]$tax_name=="metagenome"){
          ena_metadata[itx,]$tax_name <- tax_name[itx]
        }
      }
    }
  }else{
    ena_metadata$tax_name <- as.character(tax_name) #just to create the column
  }
  # first fill in the scientific_name, already present as ena_metadata$tax_name
  ena_variable <- c(ena_variable, "scientific_name")
  ena_units <- c(ena_units, "")
  
  # then convert scientific_name to its tax_id
  ena_metadata$tax_IDs <- ena_metadata$tax_name
  ena_metadata$tax_IDs <- sapply(ena_metadata$tax_IDs, FUN=function(x){commonTax.to.NCBI.TaxID(x)})
  ena_variable <- c(ena_variable, "tax_id")
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
      cat("No sample description was found:\n\tPlease type a sample description. Or type n to ignore.\n")
      descr <- readline() 
      if(!descr %in% c("n", "N")){
         ena_metadata$sample_description<-rep(descr, nrow(metadata))
      }else{
        ena_metadata$sample_description<-rep("", nrow(metadata))
      }
    }
  }
  ena_variable <- c(ena_variable, "sample_description")
  ena_units <- c(ena_units, "")
  
  # 5.6 geographic location (to get it on the right place)
  if("geo_loc_name" %in% colnames(metadata)){
    # for geographic location, only the names in the ENA_geoloc vector are allowed, oherwise ENA trows an error at submission 
    # check if valid:
    if(any(!unique(metadata$geo_loc_name) %in% ENA_geoloc)){
      # not valid location names
      geoloc_data <- data.frame(geoloc_orig = as.character(unique(metadata$geo_loc_name)), 
                                geoloc_ena = as.character(unique(metadata$geo_loc_name)),
                                stringsAsFactors = FALSE)
      geoloc_data <- geoloc_data[!geoloc_data$geoloc_orig=="",] #empty values will be handled at the end
      if(any(grepl(":", geoloc_data$geoloc_orig)) | 
         any(grepl("|", geoloc_data$geoloc_orig)) |
         any(grepl(",", geoloc_data$geoloc_orig)) |
         any(grepl(";", geoloc_data$geoloc_orig))){
        # probably a hierarchical order of nested locations 
        geoloc_data$geoloc_orig <- sapply(geoloc_data$geoloc_orig, FUN=function(x){
          gsub("|",":",x, fixed=TRUE)})
        geoloc_data$geoloc_orig <- sapply(geoloc_data$geoloc_orig, FUN=function(x){
          gsub(",",":",x)})
        geoloc_data$geoloc_orig <- sapply(geoloc_data$geoloc_orig, FUN=function(x){
          gsub(";",":",x)})
        for(l in 1:nrow(geoloc_data)){
          loc <- strsplit(geoloc_data[l,1], ":")[[1]]
          loc <- trimws(loc, which="both")
          i_loc <- intersect(tolower(loc), tolower(ENA_geoloc))
          if(length(i_loc)<1){
            bad_input <- TRUE
          }else if(length(i_loc)>1){
            if("antarctica" %in% tolower(i_loc) & "southern ocean" %in% tolower(i_loc)){
              i_loc <- "Southern Ocean"
              bad_input <- FALSE
            } else{
              bad_input <- TRUE
            }
          }else{
            i_loc <- ENA_geoloc[which(grepl(i_loc, tolower(ENA_geoloc)))]
            bad_input <- TRUE
          }
          if(bad_input & ask.input){
            if(ask.input){
              while(bad_input){
                cat(paste("Illegal geographic location name found for \"",geoloc_data[l,1],"\".\n\tPlease type the correct name.\n\tType n to ignore.\n\tType h to see the allowed terms first.", sep=""))
                loc_input <- readline()
                if(!loc_input %in% c("n", "N", "h", "H") & loc_input %in% ENA_geoloc){
                  i_loc<- loc_input
                  bad_input <- FALSE
                }else if(loc_input %in% c("n", "N")){
                  i_loc <- ""
                  bad_input <- FALSE
                }else if(loc_input %in% c("h", "H")){
                  print(ENA_geoloc)
                }
              }
            }else{
              i_loc <- ""
            }
          }
          geoloc_data[l,2] <- as.character(i_loc)
        }
      }
      ena_metadata$geo_loc_name <- unname(unlist(sapply(as.character(metadata$geo_loc_name), FUN=function(x){
        if(x %in% geoloc_data$geoloc_orig){
          x<-geoloc_data[geoloc_data$geoloc_orig==x,]$geoloc_ena
        }else{
          x<-""
        }
        return(x)
        })))
      
    }else{
      ena_metadata$geo_loc_name <- metadata$geo_loc_name
    }
    ena_variable <- c(ena_variable, get.ENAName("geo_loc_name"))
    ena_units <- c(ena_units, "")
  } else{
    ena_metadata$geo_loc_name <- rep("not collected", nrow(ena_metadata))
    ena_variable <- c(ena_variable, get.ENAName("geo_loc_name"))
    ena_units <- c(ena_units, "")
    }
  
  # 5.7 latlon DD
  if("decimalLatitude" %in% colnames(metadata) & 
     "decimalLongitude" %in% colnames(metadata)){
    ena_metadata$decimalLatitude <- metadata$decimalLatitude
    ena_metadata$decimalLongitude <- metadata$decimalLongitude
  }else if("lat_lon" %in% colnames(metadata)){
    metadata <- separate(metadata, "lat_lon", sep=" ", into=c("decimalLatitude", "decimalLongitude"))
    ena_metadata$decimalLatitude <- metadata$decimalLatitude
    ena_metadata$decimalLongitude <- metadata$decimalLongitude
  }else{
    warningmessages <- multi.warnings("The data lacks geographical coordinates.", warningmessages)
    ena_metadata$decimalLatitude <- rep("not collected", nrow(ena_metadata))
    ena_metadata$decimalLongitude <- rep("not collected", nrow(ena_metadata))
  }
  ena_variable <- c(ena_variable, get.ENAName("decimalLatitude"), get.ENAName("decimalLongitude"))
  ena_units <- c(ena_units, "DD", "DD")
   
  # 5.8 convert remaining MIxS terms to ENA variants
  redundant_cols <- which(colnames(metadata) %in% c("sample_description", "lat_lon", "decimalLatitude", 
                                                   "decimalLongitude", "geo_loc_name", 
                                                   "original_name", "subspecf_gen_lin", "scientific_name"))
  metadata_reduced <- metadata[,-redundant_cols]
  metaunits <- metaunits[-redundant_cols]
  for(cl in 1:ncol(metadata_reduced)){
    clName <- colnames(metadata_reduced[cl])
    ena_name <- get.ENAName(clName)
    if(length(ena_name) == 0 || ena_name==""){
      ena_name <- clName
    }
    ena_metadata[,clName] <- c(as.character(metadata_reduced[,colnames(metadata_reduced) %in% clName]))
    ena_variable <- c(ena_variable, ena_name)
    ena_units <- c(ena_units, gsub("alphanumeric", "", metaunits[cl]))
  }
  
  # 6. convert missing data for required fields to ENA accepted term "not collected"
  ENArequired_terms <- c("lat_lon", "decimalLatitude", "decimalLongitude", "geo_loc_name",
                      "collection_date", "seq_meth", "investigation_type", "alt",
                      "elev", "depth")
  for(required_val in ENArequired_terms){
    if(required_val %in% colnames(ena_metadata)){
      ena_metadata[,required_val]<- sapply(ena_metadata[,required_val], function(x){
        x<-as.character(x)
        if(x==""){
          x<-"not collected"}
        return(x)})
    }
  }
  if("host_taxid" %in% colnames(ena_metadata)){
    ena_metadata$host_taxid<- sapply(ena_metadata$host_taxid, function(x){
      x<-as.character(x)
      if(x==""){
        x<-"NCBI:txid12908"}
      return(x)})
  }
  
  
  # 7. the runInfo data (= technical details of the samples)
  # 7.1. sample_alias and creating the runInfo data frame
  ena_runInfo <- data.frame(sample_alias = ena_metadata$sample_alias, 
                            instrument_model = "", library_name = "", library_source = "", 
                            library_selection = "", library_strategy = "", library_layout = "",
                            design_description = "", library_construction_protocol = "", 
                            insert_size = "", forward_file_name = "", forward_file_md5 = "", 
                            reverse_file_name = "", reverse_file_md5 = "", 
                            stringsAsFactors = FALSE)
    
  # 7.2 instrument_model
  if("seq_meth" %in% colnames(metadata)){
    ena_runInfo$instrument_model <- metadata$seq_meth
  }
  
  if(any(! unique(ena_runInfo$instrument_model) %in% ENA_instrument) |
     !"seq_meth" %in% colnames(metadata) &
     ask.input){
    bad_input <- TRUE
    while(bad_input){
      cat("No valid sequencing instrument model found:\n\tPlease type the correct instrument model.\n\tType n to ignore.\n\tType h to see the allowed instrument models first.\n")
      descr <- readline() 
      if(!descr %in% c("n", "N", "h", "H") & descr %in% ENA_instrument){
        ena_runInfo$instrument_model <- rep(descr, nrow(ena_runInfo))
        bad_input <- FALSE
      } else if(descr %in% c("h", "H")){
        print(ENA_instrument)
      } else if(descr %in% c("n", "N")){
        ena_runInfo$instrument_model <- rep("", nrow(ena_runInfo))
        bad_input <- FALSE
      }
    }
  }

  
  # 7.3. library_name	
  if("original_name" %in% colnames(metadata)){
    ena_runInfo$library_name <- metadata$original_name
  }
  
  # 7.4. library_source 
  if("library_source" %in% colnames(metadata)){
    ena_runInfo$library_source <- toupper(metadata$library_source)
  }else if(ask.input){
    bad_input<-TRUE
    while(bad_input){
      cat("Could the data be cathegorized as:\n\ta) metagenomic\n\tb) metatranscriptomic\n\tc) genomic\n\td)viral RNA\n\te)other\n\tType one of the letters (a-e) or n to ignore.\n")
      descr <- readline() 
      if(tolower(descr) %in% c("a", "b", "c", "d", "e")){
        if(tolower(descr)=="a"){
          descr <- "METAGENOMIC"
        }else if(tolower(descr)=="b"){
          descr <- "METATRANSCRIPTOMIC"
        }else if(tolower(descr)=="c"){
          descr <- "GENOMIC"
        }else if(tolower(descr)=="d"){
          descr <- "VIRAL RNA"
        }else if(tolower(descr)=="e"){
          descr <- "OTHER"
        }
        ena_runInfo$library_source <- rep(descr, nrow(ena_runInfo))
        bad_input<-FALSE
      }else if(descr %in% c("n", "N")){
        bad_input<-FALSE
      }
    }

  }else if("investigation_type" %in% colnames(metadata)){
    # metagenomic
    ena_runInfo$library_source <- gsub("mimarks-survey", "METAGENOMIC", metadata$investigation_type)
    ena_runInfo$library_source <- gsub("metagenome", "METAGENOMIC", ena_runInfo$library_source)
  }
  
  if(any(! unique(ena_runInfo$library_source)  %in% c("GENOMIC", "METAGENOMIC",
                                                      "OTHER", "TRANSCRIPTOMIC",
                                                      "METATRANSCRIPTOMIC", "VRIAL RNA",
                                                      "SYNTHETIC", ""))){
    stop("Invalid input for \"library_source\"\n\tOnly following terms allowed: GENOMIC, METAGENOMIC,\n\tOTHER, TRANSCRIPTOMIC, METATRANSCRIPTOMIC, VRIAL RNA, or SYNTHETIC")
  }

  # 7.5. library_selection
  if(!is.na(library.selection) && library.selection %in% ENA_select){
    ena_runInfo$library_selection <- rep(library.selection, nrow(ena_runInfo))
  }else if(ask.input){
    bad_input <- TRUE
    while(bad_input){
      cat("What method was used to select for, enrich or or screen the material being sequenced (e.g. PCR)?\n\tType the appropriate method.\n\tType n to ignore.\n\tType h to see all the allowed methods.\n")
      descr <- readline() 
      if(!descr %in% c("n", "N", "h", "H") & descr %in% ENA_select){
        ena_runInfo$library_selection <- rep(descr, nrow(ena_runInfo))
        bad_input <- FALSE
      } else if(descr %in% c("h", "H")){
        print(ENA_select)
      } else if(descr %in% c("n", "N")){
        ena_runInfo$library_selection <- rep("unspecified", nrow(ena_runInfo))
        bad_input <- FALSE
      }
    }
  }
  
  # 7.6. library_strategy
  if(!is.na(library.strategy) && library.strategy %in% ENA_strat){
    ena_runInfo$library_strategy <- rep(library.strategy, nrow(ena_runInfo))
  }else if("run_type" %in% colnames(metadata)){
    ena_runInfo$library_strategy <- metadata$run_type
  }else if("library_strategy" %in% colnames(metadata)){
    ena_runInfo$library_strategy <- metadata$library_strategy
  }
  if(any(! unique(ena_runInfo$library_strategy) %in% ENA_strat) &
     ask.input){
    bad_input <- TRUE
    while(bad_input){
      cat("No valid library strategy found:\n\tPlease type the correct strategy (e.g. AMPLICON, WGS,...).\n\tType n to ignore.\n\tType h to see the allowed strategies.\n")
      descr <- readline() 
      if(!descr %in% c("n", "N", "h", "H") & descr %in% ENA_strat){
        ena_runInfo$library_strategy <- rep(descr, nrow(ena_runInfo))
        bad_input <- FALSE
      } else if(descr %in% c("h", "H")){
        print(ENA_strat)
      } else if(descr %in% c("n", "N")){
        ena_runInfo$library_strategy <- rep("", nrow(ena_runInfo))
        bad_input <- FALSE
      }
    }
  }
  
  # 7.7. library_layout SINGLE or PAIRED
  if(!is.na(library.layout) && library.layout %in% c("PAIRED", "SINGLE")){
    ena_runInfo$library_layout <- rep(library.layout, nrow(ena_runInfo))
  }else if("library_layout" %in% colnames(metadata)){
    ena_runInfo$library_layout <- metadata$library_layout
  }else if(ask.input){
    bad_input<-TRUE
    while(bad_input){
      cat("What was the library layout?.\n\ta) paired-end\n\tb) single-end\n\tType the appropriate letter (a-b).\n")
      descr <- readline() 
      if(descr %in% c("a", "A", "b", "B")){
        descr <- gsub("a", "PAIRED", tolower(descr))
        descr <- gsub("b", "SINGLE", tolower(descr))
        ena_runInfo$library_layout <- rep(descr, nrow(ena_runInfo))
        bad_input<-FALSE
      }else{
        cat("wrong input\n")
      }
    }
  }
  
  # 7.8. design_description
  if("experimental_factor" %in% colnames(metadata)){
    ena_runInfo$design_description <- metadata$experimental_factor
  }
  
  # 7.9. library_construction_protocol
  if("lib_const_meth" %in% colnames(metadata)){
    ena_runInfo$library_construction_protocol <- metadata$lib_const_meth
  }
  
  # 7.10. insert_size
  ## required parameter
  if(!is.na(insert.size)){
    ena_runInfo$insert.size <- rep(insert.size, nrow(ena_runInfo))
  }else if("insert_size" %in% colnames(metadata)){
    ena_runInfo$insert_size <- metadata$insert_size
  }else if(ask.input){
    cat("What was the insert size of the reads?.\n\tType the appropriate number, or type n to ignore.\n")
    descr <- readline()
    if(grepl("[0-9]",descr)){
      ena_runInfo$insert_size <- rep(descr, nrow(ena_runInfo))
    }
  }
  
  # 7.11. forward_file_name / reverse_file_name
  if("original_name" %in% colnames(metadata)){
    # assuming the original names corresponds to the sequence data file names
    for(i in 1:nrow(ena_runInfo)){
      if(ena_runInfo[i,]$library_layout == "PAIRED"){
        # assuming _1 for forward files and _2 for reverse files
        ena_runInfo[i,colnames(ena_runInfo)=="forward_file_name"] <- paste(metadata[i,]$original_name, "_1", seq.file.extension, sep="")
        ena_runInfo[i,colnames(ena_runInfo)=="reverse_file_name"] <- paste(metadata[i,]$original_name, "_2", seq.file.extension, sep="")
      }else if(ena_runInfo[i,]$library_layout == "SINGLE"){
        ena_runInfo[i,colnames(ena_runInfo)=="forward_file_name"] <- paste(metadata[i,]$original_name, seq.file.extension, sep="")
        ena_runInfo[i,colnames(ena_runInfo)=="reverse_file_name"] <- paste(metadata[i,]$original_name, seq.file.extension, sep="")
      }
    }
  }
  

  # 8. put everything together for the sample metadata
  # 8.1 lines 3 to 5
  out.line03 <- paste(ena_variable, collapse='\t')
  out.line04 <- paste("#template", paste(rep("", ncol(ena_metadata)-1), collapse="\t"), sep='\t') #paste(c(as.character(ena_metadata[1,][-1]))
  out.line05 <- paste(ena_units, collapse='\t')
  
  # 8.2 line >5
  out.line99 <- apply(ena_metadata, 1, function(x) paste(x, collapse = "\t")) 
  out.line99 <- paste(out.line99, collapse="\n")
  
  # 8.3 combine everything
  SampleMetadata_output <- paste(out.line01, out.line02, out.line03, out.line04, out.line05, out.line99, sep="\n")
  
  # 8.4 the runinfodata
  out.runInfo1 <- paste(colnames(ena_runInfo), collapse = "\t")
  out.runInfo2 <- apply(ena_runInfo, 1, function(x) paste(x, collapse = "\t")) 
  out.runInfo2 <- paste(out.runInfo2, collapse="\n")
  runInfo_output <- paste(out.runInfo1, out.runInfo2, sep="\n")
  
  # 9. finalize
  metadataFile <- paste(dest.dir, "/", file.name, ".tsv", sep="")
  runInfoFile <- paste(dest.dir, "/", file.name, "_runInfo.tsv", sep="")
  
  write.table(SampleMetadata_output, file=metadataFile, 
              col.names = FALSE, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  cat(paste("The sample metadata has been written to ", metadataFile, "\n",sep=""))
  
  write.table(runInfo_output, file=runInfoFile, 
              col.names = FALSE, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  cat(paste("The technical runInfo data has been written to ", runInfoFile, "\n",sep=""))
  
  


  
}



#--------------------------------------------------------------
# 4. tools for getting help
#--------------------------------------------------------------

term.definition <- function(term){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description term.defnition retrieves the definition of a variable term of MIxS or DarwinCore used in the POLAAAR portal
  #' @param term a character string. The term of which you would like to see the definition.
  #' @details Standerdizing microbial sequence data, metadata and environmental data can be quite difficult given the plethora of standard terms already in existance. This function returns a definition of any term that is used on the POLAAAR portal at biodiversity.aq. 
  #' @return the definition of a term, depending on the origin of the term, which can be 1) MIxS, 2) DarwinCore or 3) the POLAAAR team

  if(term %in% TermsLib$name){
    def_out <- as.character(TermsLib[TermsLib$name==term,]$definition)
  }else{
    for(ls in TermsSyn){
      if(grepl(term, ls)){
        ### still needs to be finished!!
        print(ls)
      }
    }
  }
  
  
  
  cat(def_out)
}


#--------------------------------------------------------------
# test zone
#--------------------------------------------------------------

Heind <- read.csv("/Users/msweetlove/Desktop/historic_fish/MiMARKS_Heindler_PRJEB34858.csv", row.names=1)
####!!!!!!!!pictures on zenodo
HeindQC <- process.metadata(metadata = Heind, strict.MIxS = FALSE, 
                            out.format="metadata.MIxS", ask.input=TRUE)


prep.metadata.ENA(metadata=HeindQC, dest.dir="/Users/msweetlove/Desktop/historic_fish",
                  file.name="MiMARKS_Heindler_PRJEB34858",
                  sample.unique_name_prefix="HistoricAntarcticFishDataset_2019_", 
                  checklist_accession=NA, tax_name=NA, ask.input=TRUE,
                  insert.size = NA, library.layout = "PAIRED",
                  library.strategy = "AMPLICON", library.selection = "PCR")


testdir <- "/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/test"


cesar <- read.csv("/Users/msweetlove/Desktop/cesar.csv", row.names=1)
head(cesar)
cesarQC <- process.metadata(metadata = cesar, strict.MIxS = FALSE, 
                            out.format="metadata.MIxS", ask.input=T)

head(cesarQC@data)
write.metadata.MIxS.to.mars(cesarQC, name.prefix="cesar", out.dir=testdir,
                                        add.missing.data.columns=TRUE, ask.input=FALSE)

write.metadata.MIxS.as.ENA(metadata=HeindQC, name="/Users/msweetlove/Desktop/historic_fish/MiMARKS_Heindler_PRJEB34858.tsv",
                           unique_name_prefix="HistoricAntarcticFishDataset_2019_", checklist_accession=NA,
                           ask.input=TRUE, missing.data.as.empty.columns=TRUE)





testdir="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir/test"


wdir2="/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/R_wDir"
setwd(wdir2)

#test BioProject
#BP<-"PRJNA369175"
#dataat polaRRTools::MIxSTerms.RData
#apiKey<-"5e490715dbb88b3f861565b05aab426f1408"
#test<-process.metadata(metadata=MetaData, strict.MIxS = FALSE)

#get.BioProject.metadata.INSDC(BioPrjct="PRJNA369175", just.names=TRUE)




cesar2<- get.sample.attributes.INSDC(sampleID=NA, apiKey="5e490715dbb88b3f861565b05aab426f1408",
                                   BioPrjct=c("PRJNA541486"))
test2 <- process.metadata(metadata = test1, add_to = NA, strict.MIxS = FALSE, 
                          out.format="metadata.MIxS", ask.input=TRUE)


test1<- get.sample.attributes.INSDC(sampleID=NA, apiKey="5e490715dbb88b3f861565b05aab426f1408",
                                    BioPrjct=c("PRJNA395930"))

test2 <- process.metadata(metadata = test, add_to = NA, strict.MIxS = FALSE, 
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


########
## METHANOBASE


#Methanobase_Metadata_MIxS

# change file names
dir_methano <- "/Users/msweetlove/Desktop/methanobase_seq"
files <- list.files(path=dir_methano, pattern="*.fastq.gz")

methano <- read.csv("/Users/msweetlove/Desktop/Methanobase_Metadata_MIxS.csv", 
                    row.names=1, header=TRUE)

methanoQC <- process.metadata(metadata = methano, strict.MIxS = FALSE, 
                            out.format="metadata.MIxS", ask.input=TRUE)

prep.metadata.ENA(metadata=HeindQC, dest.dir="/Users/msweetlove/Desktop/historic_fish",
                  file.name="MiMARKS_Heindler_PRJEB34858",
                  sample.unique_name_prefix="HistoricAntarcticFishDataset_2019_", 
                  checklist_accession=NA, tax_name=NA, ask.input=TRUE,
                  insert.size = NA, library.layout = "PAIRED",
                  library.strategy = "AMPLICON", library.selection = "PCR")



for(f in files){
  f_pth <- paste(dcoi, f, sep="/")
  f2 <- gsub("-1.fastq.gz", "_1.fastq.gz", f, fixed=TRUE)
  f2 <- gsub("-2.fastq.gz", "_2.fastq.gz", f, fixed=TRUE)
  
  fn <- seqlist[seqlist$oldName==f2,]$NewName
  if(length(fn)>0){
    fn2 <- paste(dcoi, "/", fn, ".fastq.gz", sep="")
    file.rename(f_pth, fn2)
  }
  
}



###########



