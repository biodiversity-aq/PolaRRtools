#==============================================================
# The polaRRtools package
#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2019-09-20)
# file encdong UTF-8
#
#==============================================================
# Script with tools to convert between data formats
#
#==============================================================

MIxS.to.DwC_eMOF <- function(metadata.object)

eMoF.to.wideTable <- function(dataset){
  ## converts an extended Measurement of Fact table to a regular wide table
  
  require(tidyr)
  # check input
  if(!all(c("measurementType", "measurementValue", "measurementUnit") %in% colnames(dataset))){
    stop("Invalid DwC eMoF file.")
  }

  if("eventID" %in% colnames(dataset)){
    dataLong<-dataset[,colnames(dataset) %in% c("eventID", "measurementType", "measurementValue")]
    dataWide <- data.frame(tidyr::spread(dataLong, eventID, measurementValue))
  }else if("occurrenceID" %in% colnames(dataLong)){
    dataLong<-dataset[,colnames(dataset) %in% c("occurrenceID", "measurementType", "measurementValue")]
    dataWide <- data.frame(tidyr::spread(dataLong, occurrenceID, measurementValue))
  }else{
    stop("Invalid DwC eMoF file.")
  }
  
  rownames(dataWide) <- dataWide$measurementType
  dataWide <- data.frame(t(dataWide[,!colnames(dataWide) %in% "measurementType"]))
  
  dataUnits<-dataset[,colnames(dataset) %in% c("measurementType", "measurementUnit")]
  dataUnits <- unique(dataUnits)
  Units_list<-c()
  Units_list[as.character(dataUnits$measurementType)]<-as.character(dataUnits$measurementUnit)
  
  Method_list<-c()
  if("measurementMethod" %in% colnames(dataset)){
    dataMethods<-dataset[,colnames(dataset) %in% c("measurementType", "measurementMethod")]
    dataMethods <- unique(dataMethods)
    Method_list[as.character(Method_list$measurementType)]<-as.character(Method_list$measurementMethod)
  }

  # not yet included:
  #measurementAccuracy
  #measurementDeterminedDate
  #measurementDeterminedBy
  #measurementRemarks
  #measurementID
  #measurementTypeID
  #measurementValueID
  #measurementUnitID
  
  return(list(data=dataWide, units=Units_list, method=Method_list))
  
}


preprocess.polaaar.DwC.eMoF <- function(dataset){
  dataset <- eMoF.to.wideTable(dataset)
  return(dataset)

}






preprocess.polaaar.MIxS <- function(metadata.object){
  ## note: this is a preprocessing function, no QC is executed
  
  metaunits <- metadata.object@units
  metasection <- metadata.object@section
  dataset <- metadata.object@data
  
  # a. samplingProtocol
  if("samp_collect_device" %in% colnames(dataset)){
    dataset$samplingProtocol <- dataset$samp_collect_device
  }
  
  newTerms <- setdiff(colnames(dataset), names(metaunits))
  for(trm in newTerms){
    metaunits[trm] <- as.character(TermsLib[TermsLib$name=="samplingProtocol",]$expected_unit)
    metasection[trm] <- as.character(TermsLib[TermsLib$name=="samplingProtocol",]$MIxS_section)
  }
  
  New_metadata <- new("metadata.MIxS",
                      data   = dataset,
                      section = metasection,
                      units      = metaunits,
                      env_package = metadata.object@env_package,
                      type = metadata.object@type,
                      QC = metadata.object@QC
  )
  
  return(New_metadata)
}

preprocess.polaaar.DwC.core <- function(dataset){
  
  ## note: this is a preprocessing function, no QC is executed
  
  # a. create a collection_date and time field
  if("eventDate" %in% colnames(dataset)){
    dataset$collection_date <- dataset$eventDate
  }
  if("eventTime" %in% colnames(dataset)){
    dataset$collection_time <- dataset$eventTime
  }
  
  # b. geo_loc_name
  geoTerms <- intersect(colnames(dataset), c("continent", "country", "state_province", 
                                             "waterBody", "islandGroup", "island", "locality"))
  if(length(geoTerms)!=0){
    for(gt in geoTerms){
      dataset[,gt]<-paste(gt, "=", dataset[,gt], sep="")
      dataset[,gt]<-gsub(".+=$", "",dataset[,gt], fixed=FALSE)
    }
    dataset$geo_loc_name <- apply(dataset[,geoTerms], 1, function(x) paste(x[!is.na(x) & x != ""], collapse = ";", sep=""));
    dataset$geo_loc_name<-gsub(";;", ";",dataset$geo_loc_name, fixed=TRUE)
  }
  
  # c. dynamicProperties to individualcolumns
  # doing this cell by cell, because the content can be anything...
  # expected format: {"colNameWithUnits"=value}, {...}
  if("dynamicProperties" %in% colnames(dataset)){
    for(i in 1:nrow(dataset)){
      dc<-dataset[i,]$dynamicProperties
      if(grepl("\\{", dc)){
        # format brackets to one type
        dc<-gsub("\\[", "\\{", dc)
        dc<-gsub("\\]", "\\}", dc)
        dc<-gsub("\\(", "\\{", dc)
        dc<-gsub("\\)", "\\}", dc)
        
        #remove quotes
        dc<-gsub("\"", "", dc)
        dc<-gsub("\'", "", dc)
        
        #try to split fieds
        dc<-gsub("\\} \\{", "---", dc)
        dc<-gsub("\\}, \\{", "---", dc)
        dc<-gsub("\\},\\{", "---", dc)
        dc<-gsub("\\}\\{", "---", dc)
        dc<-gsub("\\{", "", dc)
        dc<-gsub("\\}", "", dc)
        dc<-strsplit(dc, "---")[[1]]
        
        for(dx in dc){
          name<-strsplit(dx, "=")[[1]][1]
          name<-gsub(" ", "", name)
          value<-strsplit(dx, "=")[[1]][2]
          value<-gsub(" ", "", value)
          if(!name %in% colnames(dataset)){
            dataset[,name]<-""
            dataset[i,name]<-value
          }else{
            dataset[i,name]<-value
          }
        }
        
      }
      
    }
    
  }
  
  # get associatedSequences
  
  ## associatedSequences to emof??
  
  return(dataset)
}
