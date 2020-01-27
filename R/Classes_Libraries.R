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
# 2. metadata.MIxS
#--------------------------------------------------------------
# 2.1 metadata.MIxS Class object
#--------------------------------------------------------------
setClass("metadata.MIxS", slots=list(data="data.frame", #a dataframe with the metadata (samples are rows, variables are columns)
                                     section="character", #a character vector containing a MIxS section for each term (variable)
                                     units="character", #a character vector containing a unit for each variable
                                     type="character", #versatile (loose following of MIxS) or strict.MIxS (following all the MIxS rules)
                                     env_package="character" #the environmental package. Can be not_specified or multiple_packages.
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
setClass("metadata.DwC", slots=list(name="character", #the name of the resource
                                    core="data.frame", #a dataframe with the DwC core (events/occurences are rows, variables are columns)
                                    emof="data.frame", #a dataframe with the DwC EMOF extention (=the environmental data)
                                    type="character", #type of the core: occurence or event
                                    EML.url="character" #the url to the EML file from the IPT
)
)

#--------------------------------------------------------------
# 3.2 metadata.DwC Methods
#--------------------------------------------------------------
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
