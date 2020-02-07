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
#==============================================================

#--------------------------------------------------------------
# 1. Loading Datsets
#--------------------------------------------------------------

#library(devtools)
#library(usethis)

# some explanation
#' @field name character. Has the central name (as a unique ID) of the term/concept.
#' @field name_variant_ENA character. The name variant used by ENA
#' @field corresponding_DwC character. Has the name variant used in DarwinCore (note: all DarwinCore terms are in the list seperately as well)
#' @field name_origin character. The origin vocabulary or standard of the term (as in the column "name"). Either MIxS (GSC), DwC (TDWG), INSDC (SRA, ENA,...) or miscellaneous
#' @field MIxS_section character. The section the term (as in the column "name") belongs to as defined by in the MIxS standard.
#' @field DwC_section character. The section the term (as in the column "name") belongs to as defined by in the DarwinCore (DwC) standard.
#' @field expected_unit character. The expected unit.
#' @field synonyms character. A list of all the synonyms the term can have, including the one in the "name" column, as well as common typo's. Synonyms in the list are separated by ";" (without any spaces), and are written in lower case (case sensitive).
#' @field definition character. definition of a term. (no commas allowed, used ";" instead)
#' @field polaaar_emof numeric (0-1). 1 = a term that can be writen to a DarwinCore extende MeasurementOrFact (eMoF) file, 0 = a term that can not be written to eMoF
#' @field polaaar_template numeric (0-1-2). 2 = a term that is listed in the user template for data submission to POLAAAR that is required, 1 = a non-required term that in the user template, 0 = a term not listed in the user template
#' @field polaaarDB_environment numeric (0-1). 1 = a term that can be stored in the Environment table of POLAAAR, 0 = term that can not be stored in the Environment table
#' @field MIxS_core numeric. (0-1-2). 2 = a MIxS core term that is required, 1 = a non-required MIxS core term, 0 =  not a MIxS core term
#' @field DwC_Event numeric. (0-1-2). 2 = a DwC EventCore term that is required, 1 = a non-required DwC EventCore term, 0 =  not a DwC EventCore term
#' @field DwC_Occurrence numeric. (0-1-2). 2 = a DwC Occurrence core (or extension) term that is required, 1 = a non-required DwC Occurrence term, 0 =  not a DwC Occurrence term
#' @field DwC_eMoF numeric. (0-1). 1 = a DwC eMoF term, 0 = not a DwC eMoF term
#' @field air numeric. (0-1-2). 2 = a required term from the MIxS air package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field built_environment numeric. (0-1-2). 2 = a required term from the MIxS built_environment package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field host_associated numeric. (0-1-2). 2 = a required term from the MIxS host_associated package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field human_associated numeric. (0-1-2). 2 = a required term from the MIxS human_associated package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field human_gut numeric. (0-1-2). 2 = a required term from the MIxS human_gut package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field human_oral numeric. (0-1-2). 2 = a required term from the MIxS human_oral package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field human_skin numeric. (0-1-2). 2 = a required term from the MIxS human_skin package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field human_vaginal numeric. (0-1-2). 2 = a required term from the MIxS human_vaginal package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field microbial_mat_biofilm numeric. (0-1-2). 2 = a required term from the MIxS microbial_mat_biofilm package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field miscellaneous_natural_or_artificial_environment numeric. (0-1-2). 2 = a required term from the MIxS miscellaneous_natural_or_artificial_environment package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field plant_associated numeric. (0-1-2). 2 = a required term from the MIxS plant_associated package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field sediment numeric. (0-1-2). 2 = a required term from the MIxS sediment package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field soil numeric. (0-1-2). 2 = a required term from the MIxS soil package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field wastewater_sludge numeric. (0-1-2). 2 = a required term from the MIxS wastewater_sludge package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
#' @field water numeric. (0-1-2). 2 = a required term from the MIxS water package, 1 = a non-required MIxS term from the that package, 0 = term does not occur in the that package
TermsLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/TermsLibrary_v3.csv", 
                   header=TRUE)
TermsSyn <- sapply(as.character(TermsLib[TermsLib$name_origin %in% c("MIxS", "INSDC", "miscellaneous"),]$synonyms), function(x){strsplit(x, ";")})
names(TermsSyn) <- TermsLib[TermsLib$name_origin %in% c("MIxS", "INSDC", "miscellaneous"),]$name 

TermsSyn_DwC <- sapply(as.character(TermsLib[TermsLib$name_origin %in% c("DwC"),]$synonyms), function(x){strsplit(x, ";")})
names(TermsSyn_DwC) <- TermsLib[TermsLib$name_origin %in% c("DwC"),]$name 
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
# 2. metadata.MIxS
#--------------------------------------------------------------
# 2.1 metadata.MIxS Class object
#--------------------------------------------------------------
setClass("metadata.MIxS", slots=list(data="data.frame", #a dataframe with the metadata (samples are rows, variables are columns)
                                     section="character", #a character vector containing a MIxS section for each term (variable)
                                     units="character", #a character vector containing a unit for each variable
                                     type="character", #versatile (loose following of MIxS) or strict.MIxS (following all the MIxS rules)
                                     env_package="character", #the environmental package. Can be not_specified or multiple_packages.
                                     QC="logical" #Has the data been Quality comtrolled (TRUE/FALSE)
)
)

#--------------------------------------------------------------
# 2.2 metadata.MIxS Methods
#--------------------------------------------------------------
setMethod("show",
          "metadata.MIxS",
          function(object) {
            N <- nrow(object@data)
            C <- ncol(object@data)
            cat(paste("a metadata.MIxS class object with", as.character(N), "samples and",
                      as.character(C), "variables\n"))
            if(object@QC){
              cat("\tThe quality of this object has been checked.\n")
            }else{
              cat("\tThe quality of this object has NOT been checked.\n")
            }
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
# 3. metadata.DwC
#--------------------------------------------------------------
# 3.1 metadata.DwC Class Object
#--------------------------------------------------------------
setClass("DwC.event", slots=list(core="data.frame", #a dataframe with the DwC event core (events are rows, variables are columns)
                                occurrence="data.frame", #a dataframe with the occurrence extension
                                emof="data.frame", #a dataframe with the DwC eMoF extention (=the environmental data)
                                EML.url="character", #the url to the EML file from the IPT
                                QC="logical" #Has the data been Quality comtrolled (TRUE/FALSE)
)
)

setClass("DwC.occurrence", slots=list(core="data.frame", #a dataframe with the DwC event core (events are rows, variables are columns)
                                      emof="data.frame", #a dataframe with the DwC eMoF extention (=the environmental data)
                                      EML.url="character", #the url to the EML file from the IPT
                                      QC="logical" #Has the data been Quality comtrolled (TRUE/FALSE)
)
)

#--------------------------------------------------------------
# 3.2 metadata.DwC Methods
#--------------------------------------------------------------
setMethod("show",
          "DwC.event",
          function(object) {
            N <- nrow(object@core)
            if(nrow(object@emof)>0 & nrow(object@occurrence)>0){
              E <- "and an occurrence and eMoF extensions"
            } else if(nrow(object@emof)>0 & nrow(object@occurrence)==0){
              E <- "and an eMoF extension"
            }else if(nrow(object@emof)==0 & nrow(object@occurrence)>0){
              E <- "and an occurrence extension"
            }else{ E <- ""}
            cat(paste("a  DwC.event class  object with a core of ", as.character(N), 
                      " events", as.character(E), ".\n", sep=""))
            if(EML.url!=""){
              cat(paste("\tThe metadata can be found at", EML.url, 
                        "\n", sep=""))
            }
            if(object@QC){
              cat("\tThe quality of this object has been checked.\n")
            }else{
              cat("\tThe quality of this object has NOT been checked.\n")
            }
          }
)
setMethod("show",
          "DwC.occurrence",
          function(object) {
            N <- nrow(object@core)
            if(nrow(object@emof)>0){
              E <- "and an eMoF extension"
            }else{ E <- ""}
            cat(paste("a  DwC.occurrence  class object with a core of ", as.character(N), 
                      " occurrences", as.character(E), ".\n", sep=""))
            if(EML.url!=""){
              cat(paste("\tThe metadata can be found at", EML.url, 
                        "\n", sep=""))
            }
            if(object@QC){
              cat("\tThe quality of this object has been checked.\n")
            }else{
              cat("\tThe quality of this object has NOT been checked.\n")
            }
          }
)



check.valid.metadata.DwC <- function(d){
  valid <- TRUE
  #class must be DwC.event or DwC.occurrence 
  if(!class(d) %in% c("DwC.event", "DwC.occurrence")){
    valid <- FALSE
  }else{
    dcol <- ncol(d@core)
    drow <- nrow(d@core)
    #The core must have at least 1 sample (row) and 1 variable (column)
    if(dcol<=0 & drow<=0){
      valid <- FALSE
    }
    #type must be occurence or event
    if(class(d) == "DwC.occurrence" &
       !"occurrenceID" %in% colnames(d@core) |
       !"basisOfRecord" %in% colnames(d@core)){
      valid <- FALSE
    }else if(class(d) == "DwC.event" &
             !"eventID" %in% colnames(d@core)){
      valid <- FALSE
    }
  }
  return(valid)
}
