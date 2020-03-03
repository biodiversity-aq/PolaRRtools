###library() the polartools package for the QC functions

MIxS_packages<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/MIxSPackageLibrary.csv", 
                  header=TRUE)

get.polaaar.EML.data <- function(EML.url){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description get.polaaar.EML.data extracts the information needed for the POLAAAR database from a given onine EML document
  #' @param EML.url the URL to the EML document
  #' @return a list with the data needed for the POLAAAR database. This includes: abstract, start_date, end_date, bounding_box, is_public and associated_references

  require(EML)
  require(emld)

  eml_file <- emld::as_emld(EML.url, from = "xml")
  
  dataset_name <- unname(unlist(eml_get(eml_file,"title", from="xml"))["title"])
  
  # abstract
  if(any(grepl("abstract", as.character(eml_file)))){
    abstract <- eml_get(eml_file,"abstract", from="xml")$para
    if(is.null(abstract)){
      abstract <- ""
    }
  }else{
    abstract <- ""
  }
  
  # start_date | end_date
  if(any(grepl("formationPeriod", as.character(eml_file)))){
    dates <- eml_get(eml_file,"formationPeriod", from="xml")
    dates <- unlist(dates)[1]
    start_date  <- strsplit(dates, " ")[[1]][1]
    end_date  <- strsplit(dates, " ")[[1]][2]
  } else if(any(grepl("temporalCoverage", as.character(eml_file)))){
    dates <- eml_get(eml_file,"temporalCoverage", from="xml")
    if("rangeOfDates" %in% names(dates)){
      start_date  <- gsub('\n| |\t', '', dates$rangeOfDates$beginDate)
      end_date  <- gsub('\n| |\t', '', dates$rangeOfDates$endDate)
    } else if("singleDateTime" %in% names(dates)){
      start_date  <- gsub('\n| |\t', '', dates$singleDateTime)
      end_date <- start_date
    }else{
      start_date  <- ""
      end_date <- ""
    }
  }
  
  # bounding_box
  if(any(grepl("boundingCoordinates", as.character(eml_file)))){
    bounding_box <- unlist(EML::eml_get(eml_file, "boundingCoordinates", from = "xml"))
    bounding_box <- paste("SRID=4326;POLYGON ((",
                          bounding_box["northBoundingCoordinate"], ", ",
                          bounding_box["eastBoundingCoordinate"], ", ",
                          bounding_box["southBoundingCoordinate"], ", ",
                          bounding_box["westBoundingCoordinate"], "))", sep="")
  }else{
    bounding_box<-""
  }
    
  # is_public
  if(any(grepl("pubDate", as.character(eml_file)))){
    is_public <- TRUE
  }else{
    is_public <- FALSE
  }

  # associated_references
  if(any(grepl("bibliography", as.character(eml_file)))){
    associated_references <- unlist(EML::eml_get(eml_file, "bibliography", from = "xml"))["citation"]
    associated_references <- unname(associated_references)
  }else{
    associated_references <- ""
  }

  output<-list(name = dataset_name,
               abstract = abstract,
               start_date = start_date,
               end_date = end_date,
               bounding_box = bounding_box,
               is_public = is_public,
               associated_references = associated_references
  )
  
  return(output)
}

metadata.to.polaaar <- function(metadata.object, EML.url=NA, user_ID="", 
                                dest.dir=getwd(), create.dir=FALSE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description metadata.to.polaaar writes a metadata.MIxS or metadata.DwC object to the polaaar database
  #' @param metadata.object a metadata.MIxS class object or metadata.DwC. The object to be formatted to fit the polaaar database schema
  #' @param EML.url a named vector. URL to the EML (IPT) associated with the dataset, written as a named vector, e.g. c(project1="http://url1", project2="http://url1"). The names in the vector must correspond to the different projects listed under "project_name" in the metadata object. If only one link is provided and is unnamed, than this will automatically be the parent project to all the samples in the metadata object.
  #' @param user_ID ID of the user that send the data
  #' @param dest.dir a character vector. The complete path to the directory where the CSV files should be written to. Default is the working diectory. If create.dir is TRUE, the given directory will be created if it does not exist.
  #' @param create.dir boolean. If TRUE, the directory path in the variable dest.dir will be created if it does not already exist. Default FALSE
  #' @details 
  #' @return data is written to polaaar

  # important! This function assumes data has been QC'd by process.metadata()
  
  require(RCurl)
  baseName <- deparse(substitute(metadata.object))

  # 1. check input
  #****************************************************************************
  if(check.valid.metadata.MIxS(metadata.object)){
    metadata.object <- preprocess.polaaar.MIxS(metadata.object)
    metaformat <- "MIxS"
    metaunits <- metadata.object@units
    metasection <- metadata.object@section
    metapackage <- metadata.object@env_package
    coredata <- metadata.object@data
    env_metadata <- coredata
  }else if(check.valid.metadata.DwC(metadata.object)){
    metaformat <- "DwC"
    metapackage <- class(metadata.object)
    coredata <- preprocess.polaaar.DwC.core(metadata.object@core)
    emof <- eMoF.to.wideTable(metadata.object@emof)
    metaunits <- emof$units
    env_metadata <- emof$data
    if(!all(rownames(coredata) %in% rownames(env_metadata))){
      newRowNames <- setdiff(rownames(coredata), rownames(env_metadata))
      newRows <- matrix("", nrow=length(newRowNames), ncol=ncol(env_metadata))
      rownames(newRows) <- newRowNames
      colnames(newRows) <- colnames(env_metadata)
      env_metadata <- rbind(env_metadata, newRows)
    }
    env_metadata<-env_metadata[rownames(coredata),]#make order the same as coredata
    if(metapackage=="DwC.occurrence"){
      occurdata <- coredata
    }else{
      occurdata <- metadata.object@occurrence
    }
    
  }else{
    stop("metadata.object argument must be a metadata.MIxS, DwC.event or DwC.occurrence class object.\n\tCommon dataframes should be converted to metadata.MIxS or metadata.DwC with \n\tthe data.frame.to.metadata.MIxS function.")
  } 
  
  # the data must be associated with an EML reccord, otherwise it should not be added to the database
  if(is.na(EML.url)){
    if(metaformat=="DwC" & !is.na(metadata.object@EML.url)){
      EML.url <- metadata.object@EML.url
    }else{
      stop("No EML.url provided.\n\tAny dataset must be associated with a project, of which the project metadata has been\n\tsubmitted to the IPT.")
    }
  }
  # multiple EML URLs must be named
  if(length(EML.url)>1 & !is.null(names(EML.url))){
    stop("If multiple EML URLs provided, they must be named with names corresponding to\n\tproject names in the \"project_name\" column of the metadata object")
  }

  # check if all EML.urls lead to an existing web page
  for(emlurl in EML.url){
    if(!RCurl::url.exists(emlurl)){
      StopMessage <- paste("The following EML page does not seem to exist:\n\t", emlurl, sep="")
      stop(StopMessage)
    }
  }
  
  num_recs <- nrow(coredata)

  # check output destination
  if(!dir.exists(file.path(dest.dir))){
    if(create.dir){
      dir.create(dest.dir)
    }else{
      stop("the output directory does not exist.\n\tCheck for typo's.")
    }
  }
  
  # 2. filling the Polaaar database (PDB) tables
  #****************************************************************************
  #------------------------------------------
  #### create Event table ####
  #------------------------------------------
  # Event $ parent_event [ref: ParentEvent $ id]
  # Event $ metadata [ref: Metadata $ id]
  # Event $ occurrence [ref: occurrence $ id]
  # Event $ environment [ref: Environment $ id]
  
  EventTerms <- c("id", "footprintWKT", "eventRemarks", "sample_name", "collection_date", 
                  "collection_time", "parent_event", "parent_sample", 
                  "samplingProtocol", "occurrence", "metadata", "environment", "metadata_exists",
                  "occurrence_exists", "environment_exists")
  Event <- data.frame(matrix(nrow=num_recs, ncol=length(EventTerms), data=""))
  colnames(Event) <- EventTerms
  rownames(Event) <- rownames(coredata)
  
  
  # Event $ id (should be generated by the PDB)
  Event$id <- 1:num_recs
  
  # Event $ footprintWKT = coordinates in well-known text format
  if("footprintWKT" %in% colnames(coredata)){
    Event$footprintWKT<-coredata$footprintWKT
  }else if("decimalLatitude" %in% colnames(coredata) &
           "decimalLongitude" %in% colnames(coredata) &
           "eventID" %in% colnames(coredata)){
    Event$footprintWKT <- dataQC.generate.footprintWKT(coredata)
    Event$footprintWKT <- paste("SRID=4326;", Event$footprintWKT)
    Event$footprintWKT <- gsub("SRID=4326;$", "", Event$footprintWKT, fixed=FALSE)
    Event$footprintWKT <- gsub("SRID=4326; $", "", Event$footprintWKT, fixed=FALSE)
  }else if("lat_lon" %in% colnames(coredata)){
    Event$footprintWKT <- paste("SRID=4326; POINT (", 
                                     gsub(" ", ", ", as.character(coredata$lat_lon)), ")", 
                                     sep="")
  }
  Event$footprintWKT <- gsub("POINT (, )", "", Event$footprintWKT)
  Event$footprintWKT <- gsub("POINT (NA, NA)", "", Event$footprintWKT)
  
  # Event $ eventRemarks / collection_date / collection_time / samplingProtocol
  for(varTerm in c("eventRemarks", "collection_date", "collection_time", "samplingProtocol")){
    if(varTerm %in% colnames(coredata)){
      Event[,varTerm] <- coredata[,varTerm]
    }
  }
  
  # Event $ sample_name
  metadataNames <- dataQC.findNames(dataset = coredata, ask.input=FALSE)
  if(!any(is.na((metadataNames$Names)$original_names))){
    Event$sample_name <- (metadataNames$Names)$original_name
    warningmessages <- multi.warnings(metadataNames$warningmessages, warningmessages)
  }else{
    Event$sample_name <- colnames(coredata) 
  }
  coredata$polaaar_name <- Event$sample_name
  env_metadata$polaaar_name <- Event$sample_name
  
  #------------------------------------------
  #### create ParentEvent table ####
  #------------------------------------------
  # starts with same number of rows as dataset, but will later be reduced to unique eventIDs
  # if there are no eventIDs, the project(s) will assume the role of parent event to all the samples
  # ParentEvent $ event_type [ref: EventType $ id]
  # ParentEvent $ project [ref: ProjectMetadata $ id]
  # ParentEvent $ id [ref: Event $ parent_event]
  # hierarchy for event_type: sample [=row in coredata] > event [groups samples] > parentEvent [groups events] > project [groups samples, events or parentEvents]
  ParentEventTerms <- c("id", "parent_event_name", "event_type", "description", "parent_event",
                        "event_creator", "created_on", "updated_on", "project")
  
  if("eventID" %in% colnames(coredata)){
    event_col <- "eventID"
  }else{
    event_col <- "polaaar_name"}
  if("parentEventID" %in% colnames(coredata)){
    parentEvent_col <- "parentEventID"
  }else{
    parentEvent_col <- NA}
  if("project_name" %in% colnames(coredata)){
    project_col <- "project_name"
  }else if(!is.null(names(EML.url)) & length(EML.url)==1){
    project_col <- "project_name"
    coredata$project_name <- rep(names(EML.url), nrow(ParentEvent))
  } else if("bioproject" %in% colnames(coredata)){
    project_col <- "bioproject"
  }
  
  # get event data from the coredata
  eventTable <- dataQC.eventStructure(coredata, eventID.col = event_col, 
                                      parentEventID.col = parentEvent_col, 
                                      project.col = project_col,
                                      complete.hierarchy=TRUE)

  ParentEvent <- data.frame(matrix(nrow=nrow(eventTable), ncol=length(ParentEventTerms), data=""))
  colnames(ParentEvent) <- ParentEventTerms
  
  ParentEvent$parent_event_name <- eventTable$eventID
  ParentEvent$event_type <- eventTable$type
  ParentEvent$parent_event <- eventTable$parenEventID
  ParentEvent$project <- eventTable$project
  
  # link to Event table (first name, change to id later)
  Event$parent_event <- as.character(eventTable[eventTable$original==TRUE,]$eventID)
  Event$parent_sample <- as.character(eventTable[eventTable$original==TRUE,]$eventID)
  
  #remove duplicated events
  ParentEvent <- unique(ParentEvent)
  # ParentEvent $ id
  ParentEvent$id <- 1:nrow(ParentEvent)
  # Event $ parent_event [ref: ParentEvent $ id]
  # replace name with ID in Event table
  Event$parent_event <- unlist(sapply(Event$parent_event, FUN = function(x){gsub(x,ParentEvent[ParentEvent$parent_event_name==x,]$id,x)}))
  
  # ParentEvent $ parent_event [ref: ParentEvent $ id]
  ParentEvent$parent_event <- unlist(sapply(as.character(ParentEvent$parent_event), 
                                            FUN = function(x){
                                              if(x!=""){
                                                gsub(x,ParentEvent[ParentEvent$parent_event_name==x,]$id,x)
                                              }else{x<-""}
                                            }))


  
  #eventTable
  
  parent_event_name <- eventTable$eventID
  
  # ParentEvent $ description
  if("sample_description" %in% colnames(coredata)){
    ParentEvent[ParentEvent$parent_event_name %in% eventTable[original==TRUE,]$eventID,]$description <- coredata$sample_description
  }
  
  # ParentEvent $ created_on / updated_on
  ParentEvent$created_on <- rep(Sys.Date(), nrow(ParentEvent))
  ParentEvent$updated_on <- ParentEvent$created_on
  
  # ParentEvent $ event_creator
  # take user ID of the user that send the data, or main author of the paper if the data was harvested that way
  ParentEvent$event_creator <- rep(user_ID, nrow(ParentEvent))
  
  #------------------------------------------
  #### create ProjectMetadata table ####
  #------------------------------------------
  # information on the project level
  # preferably 1 dataset == 1 project, if so it can be perfectly linked to the EML
  # project_qaqc is a boolean to indicate if a dataset has been QCd.
  # ProjectMetadata $ id [ref: ParentEvent $ project]
  # ProjectMetadata $ associatedReferences [ref: Reference $ id]
  # project_creator = user id of person that emailed the spreadsheet
  # boundingbox notation:
  # SRID=4326;POLYGON ((N, E, S, W))
  ProjectMetadataTerms <- c("id", "project_name", "start_date", "end_date", "EML_URL", "abstract", 
                            "bounding_box", "is_public", "associated_references", 
                            "associated_media", "created_on", "updated_on", 
                            "project_creator", "project_qaqc")
  ProjectMetadata <- data.frame(matrix(nrow=length(EML.url), ncol=length(ProjectMetadataTerms), data=""), stringsAsFactors = FALSE)
  colnames(ProjectMetadata) <- ProjectMetadataTerms
  
  # some checks first
  if(length(unique(ParentEvent$project))>0 &
     !is.null(names(EML.url)) &&
     !all(names(EML.url) %in% unique(ParentEvent$project))){  
    #EML URL names must correspond to project names (but project names can have no corresponding EML name)
    #otherwise there will be records created in the ProjectMetadata table that correspond to nothing
    stop("There are EML.url names that do not correspond to project names in the metadata object.")
  }
  
  # ProjectMetadata $ id
  ProjectMetadata$id <- 1:length(EML.url) 
  
  # ProjectMetadata $ EML_URL
  ProjectMetadata$EML_URL <- EML.url
  
  # ProjectMetadata $ project_creator
  ProjectMetadata$project_creator <- rep(user_ID, nrow(ProjectMetadata))
  
  # ProjectMetadata $ project_qaqc
  # assume this to be true if it is being added to the database
  ProjectMetadata$project_qaqc <- rep(TRUE, nrow(ProjectMetadata))
  
  # ProjectMetadata $ project_name
  if(!is.null(names(EML.url))){
    # the EML.url has names
    # all EML names correspond to project names (has been checked)
    # EML names get priority, project names that do not correspond to an EML name will be droped (only in the ParentEvent$project column)
    ProjectMetadata$project_name <- names(EML.url) 
  }else{
    # no EML names, can only happen if there was just one URL provided (was checked in the beginning)
    # in that case all samples will be linked to this EML, regardless what project name has been provided
    ProjectMetadata$project_name <- "unnamed_project"
  }
  
  #### getting data from EML webpage
   for(emlurl in EML.url){
     eml_row <- as.numeric(row.names(ProjectMetadata[match(emlurl,ProjectMetadata$EML_URL),]))
     eml_data <- get.polaaar.EML.data(emlurl)
     #ProjectMetadata $ abstract
     ProjectMetadata[eml_row,"abstract"] <- eml_data$abstract
     #ProjectMetadata $ start_date
     ProjectMetadata[eml_row,"start_date"] <- eml_data$start_date
     #ProjectMetadata $ end_date
     ProjectMetadata[eml_row,"end_date"] <- eml_data$end_date
     #ProjectMetadata $ bounding_box
     ProjectMetadata[eml_row,"bounding_box"] <- eml_data$bounding_box
     #ProjectMetadata $ is_public
     ProjectMetadata[eml_row,"is_public"] <- eml_data$is_public
     #ProjectMetadata $ associated_references
     ProjectMetadata[eml_row,"associated_references"] <- eml_data$associated_references
     
     #if no name for the project was provided, change it to the dataset name found in the EML
     if(ProjectMetadata[eml_row,"project_name"] == "unnamed_project"){
       ProjectMetadata[eml_row,"project_name"] <- eml_data$name
     }
   }
  
  if("associatedReferences" %in% colnames(coredata)){
    warning("There are additional associated references in the data table that need to be checked by hand")
  }
  
  # ProjectMetadata $ created_on / updated_on
  ProjectMetadata$created_on <- rep(Sys.Date(), nrow(ProjectMetadata))
  ProjectMetadata$updated_on <- ProjectMetadata$created_on 
  
  # ProjectMetadata $ associated_media
  ### not in EML, most likely in a DwC metadata object, still need to figure out how to put this in...
  
  ### change references (if present) to IDs of Reference table
  #if(any(ProjectMetadata$associated_references != "")){
    #------------------------------------------
    #### create Reference table ####
    #------------------------------------------
    # Reference $ id [ref: ProjectMetadata $ associated_references]
    # short author list = "et al." -version
    
    # part to be resolved


  #}

  #------------------------------------------
  #### create event_type table ####
  #------------------------------------------
  # event_type $ id [ref: ParentEvent $ event_type]
  event_type <- data.frame(matrix(nrow=length(unique(ParentEvent$event_type)), ncol=2, data=""), stringsAsFactors = FALSE)
  colnames(event_type) <- c("id", "name")

  event_type$name <- unique(ParentEvent$event_type)
  event_type$id <- 1:nrow(event_type)

  # change ParentEvent $ event_type to ID
  ParentEvent$event_type <- unlist(sapply(ParentEvent$event_type, 
                                       FUN = function(x){
                                         gsub(x,event_type[event_type$name==x,]$id,x)
                                       }))


  #------------------------------------------
  #### create Metadata table ####
  #------------------------------------------
  # extra information about the samples
  # Metadata $ id [ref: Event $ metadata]
  # Metadata $ sequence [ref: Sequences $ id]
  # Metadata $ env_biome [ref: Biome $ id]

  MetadataTerms <- c("id", "metadata_tag", "md_created_on", "metadata_creator", "license", 
                     "continent", "country", "state_province", "waterBody", "islandGroup", 
                     "island", "location", "geo_loc_name", "additional_info", "env_biome", 
                     "env_package", "env_feature", "env_material", "institutionID", 
                     "nucl_acid_amp", "nucl_acid_ext", "ref_biomaterial", "rel_to_oxygen", 
                     "rightsHolder", "samp_collect_device", "samp_store_dur", "samp_store_loc",
                     "samp_store_temp", "samp_vol_we_dna_ext", "samplingProtocol", 
                     "source_mat_id", "submitted_to_insdc", "investigation_type",
                     "isol_growth_condt", "lib_size", "sequence")
  Metadata <- data.frame(matrix(nrow=num_recs, ncol=length(MetadataTerms), data=""), stringsAsFactors = FALSE)
  colnames(Metadata) <- MetadataTerms
  
  # Metadata $ md_created_on
  Metadata$md_created_on <- rep(Sys.Date(), nrow(Metadata))
  
  # Metadata $ metadata_creator
  Metadata$metadata_creator <- rep(user_ID, nrow(Metadata))
  
  # Metadata $ license
  if("license" %in% colnames(coredata)){
    Metadata$license <- coredata$license
  }else{ # use CC BY 4.0 as default (Creative Commons BY attibution, version 4.0)
    Metadata$license <- rep("CC BY 4.0", num_recs)
  }
  
  ### geographical info
  ## remark: make this a hierarchical table, like Biome?
  # Metadata $ continent, country, state_province, waterBody, islandGroup, island, location
  geoTerms <- intersect(colnames(coredata), c("continent", "country", "state_province", 
                                              "waterBody", "islandGroup", "island", "locality"))
  if(length(geoTerms)!=0){
    for(gt in geoTerms){
      if(gt != "locality"){
        Metadata[,gt] <- coredata[,gt]
      }else{
        Metadata$location <- coredata$locality
      }
    }
  }
  # Metadata $ geo_loc_name
  if("geo_loc_name" %in% colnames(coredata)){
    Metadata$geo_loc_name <- coredata$geo_loc_name
  } 
  
  # Metadata $ env_package
  #------------------------------------------
  #### Package table ####
  #------------------------------------------
  # this table is pre-defined, and is not created. Its IDs are imported.
  # Package $ id [ref: Metadata $ env_package]
  if("env_package" %in% colnames(env_metadata)){
    Metadata$env_package <- unname(sapply(env_metadata$env_package, FUN=function(x){x<-MIxS_packages[MIxS_packages$name==x,]$id}))
  } else{
    if(metapackage != "multiple_packages"){
      Metadata$env_package <- rep(MIxS_packages[MIxS_packages$name==metapackage,]$id)
    }else{
      Metadata$env_package <- rep(16, num_recs)
    }
  }
  
  # Metadata $ all the other terms
  for(varTerm in c("env_feature", "env_material", "institutionID", 
                   "nucl_acid_amp", "nucl_acid_ext", "ref_biomaterial", "rel_to_oxygen", 
                   "rightsHolder", "samp_collect_device", "samp_store_dur", "samp_store_loc",
                   "samp_store_temp", "samp_vol_we_dna_ext", "samplingProtocol", 
                   "source_mat_id", "submitted_to_insdc", "investigation_type",
                   "isol_growth_condt", "lib_size")){
    if(varTerm %in% colnames(coredata)){
        Metadata[,varTerm] <- as.character(coredata[,varTerm])
    } else if(varTerm %in% colnames(env_metadata)){
      Metadata[,varTerm] <- as.character(env_metadata[,varTerm])
    }
  }
  
  # Metadata $ id
  Metadata$id <- 1:nrow(Metadata)
  # Event $ metadata   !! asuming an event is equivalent to a sequence run (=sample)
  Event$metadata <- Metadata$id
  # Event $ metadata_exists
  Event$metadata_exists<- rep("TRUE", nrow(Event))
  
  # Metadata $ additional_info
  NoteTerms <- intersect(colnames(coredata), c("fieldNotes", "eventRemarks"))
  if(length(NoteTerms)!=0){
    Metadata$additional_info <- paste(coredata[,colnames(coredata) %in% NoteTerms], collapse=" | ")
  }
  
  #------------------------------------------
  #### create Sequences table ####
  #------------------------------------------
  # Sequences $ id [ref: Metadata $ sequence]
  SequencesTerms <- c("id", "sequence_name", "MID", "subspecf_gen_lin", "target_gene",
                       "target_subfragment", "type", "primerName_forward", 
                       "primerName_reverse", "primer_forward", "primer_reverse", "run_type",
                       "seqData_url", "seqData_accessionNumber", "seqData_projectNumber",
                       "seqData_runNumber", "seqData_sampleNumber", "seqData_numberOfBases",
                       "seqData_numberOfSequences")
  
  Sequences <- data.frame(matrix(nrow=num_recs, ncol=length(SequencesTerms), data=""), stringsAsFactors = FALSE)
  colnames(Sequences) <- SequencesTerms
  
  # Sequences $ sequence_name
  Sequences$sequence_name <- coredata$polaaar_name
  
  # Sequences $ MID
  if("mid" %in% colnames(env_metadata)){
    if("additional_mid" %in% colnames(env_metadata)){
      Sequences$MID <- paste(env_metadata$mid, env_metadata$additional_mid, collapse="; ")
    }else{
      Sequences$MID <- env_metadata$mid
    }
  }
  
  # Sequences $ type
  Sequences$type <- find.dataset(env_metadata, TermsList = c("library_source"))

  # Sequences $ primerName_forward
  Sequences$primerName_forward <- find.dataset(env_metadata, TermsList = c("forward_primer"))
  # Sequences $ primerName_reverse
  Sequences$primerName_reverse <- find.dataset(env_metadata, TermsList = c("reverse_primer"))

  
  # Sequences $ primer_forward
  if("primer_forward_sequence" %in% colnames(env_metadata)){
    Sequences$primer_forward <- env_metadata$primer_forward_sequence
    # Sequences $ primer_reverse # assume if there is no forward, there won't be a reverse
    if("primer_reverse_sequence" %in% colnames(env_metadata)){
      Sequences$primer_reverse <- env_metadata$primer_reverse_sequence
    }
  }else if("pcr_primers" %in% colnames(env_metadata)){
    primer_split <- unlist(strsplit(as.character(env_metadata$pcr_primers), " "))
    #
    #
    # still need to be worked out
    #
    #
  }

  # Sequences $ url
  Sequences$url <- find.dataset(env_metadata, TermsList = c("seqData_url"))
  
  # Sequences $ seqData_accessionNumber | seqData_sampleNumber
  accTerms <- intersect(colnames(env_metadata), c("genbank_accession_numbers", "INSDC_SampleID", "sra_run_number"))
    if(length(accTerms)!=0){
      Sequences$seqData_accessionNumber <- paste(env_metadata[,colnames(env_metadata) %in% accTerms], collapse="; ")
      Sequences$seqData_sampleNumber <- Sequences$seqData_accessionNumber
    }
  # Sequences $ seqData_projectNumber
  Sequences$seqData_projectNumber <- find.dataset(env_metadata, TermsList = c("bioproject", "project_name"))
  
  # Sequences $ seqData_runNumber
  Sequences$seqData_runNumber <- find.dataset(env_metadata, TermsList = c("sra_run_number"))

  # Sequences $ seqData_numberOfBases
  Sequences$seqData_numberOfBases <- find.dataset(env_metadata, TermsList = c("number_of_bases_predicted"))
  

  # Sequences $ seqData_numberOfSequences 
  Sequences$seqData_numberOfSequences <- find.dataset(env_metadata, TermsList = c("lib_reads_seqd"))
  
  # remaining terms
  for(varTerm in c("subspecf_gen_lin", "target_gene",  "target_subfragment", "run_type")){
    Sequences[,varTerm] <- find.dataset(env_metadata, TermsList = c(varTerm))
  }
  SeqNumbers <- Sequences$seqData_accessionNumber # set aside for later in the Occurrence table?

  # remove samples with no sequences
  Sequences <- Sequences[Sequences$seqData_url !="" | 
                           Sequences$seqData_accessionNumber !="" |
                           Sequences$seqData_projectNumber !="" |
                           Sequences$seqData_runNumber !="" |
                           Sequences$seqData_sampleNumber !="" ,]
  if(nrow(Sequences) > 0){
    # Sequences $ id
    Sequences$id <- 1:nrow(Sequences)
    # Metadata $ sequence 
    coredata$polaaar_SeqID <- coredata$polaaar_name
    for(i in 1:nrow(coredata)){
      if(coredata[i,]$polaaar_SeqID %in% Sequences$sequence_name){
        coredata[i,]$polaaar_SeqID <- Sequences[Sequences$sequence_name == coredata[i,]$polaaar_SeqID,]$id
      }else{
        coredata[i,]$polaaar_SeqID <- ""
      }
    }
    Metadata$sequence <- coredata$polaaar_SeqID
  }
  
  # Metadata $ env_biome
  if("env_biome" %in% colnames(env_metadata)){
    Metadata$env_biome <- as.character(env_metadata$env_biome)
    # assume words of one terms are separated by a " ", ".", and "_" => change this to space
    Metadata$env_biome <-gsub(".", " ", Metadata$env_biome, fixed=TRUE)
    Metadata$env_biome <-gsub("_", " ", Metadata$env_biome, fixed=TRUE)
    Metadata$env_biome<-gsub("\\s+", " ", Metadata$env_biome, fixed=FALSE)
    # change possible hierarchy separators to ":"
    Metadata$env_biome <-gsub(",", ":", Metadata$env_biome, fixed=TRUE)
    Metadata$env_biome<-gsub(">", ":", Metadata$env_biome, fixed=TRUE)
    Metadata$env_biome<-gsub(";", ":", Metadata$env_biome, fixed=TRUE)
    Metadata$env_biome<-gsub("|", ":", Metadata$env_biome, fixed=TRUE)
    Metadata$env_biome<-gsub(":+", ":", Metadata$env_biome, fixed=FALSE)
  }
  
  #------------------------------------------
  #### create Occurrence table ####
  #------------------------------------------
  # Occurrence $ id [ref: Event $ Occurrence]
  # Occurrence $ taxon [ref: Taxa $ id]
  # Occurrence $ associated_sequences [ref: Sequence $ id]
  OccurrenceTerms <- c("id", "occurrenceID", "taxon", "occurrence_notes", "occurrence_status",
                       "occurrence_class", "catalog_number", "date_identified", "other_catalog_numbers",
                       "recorded_by", "associated_sequences" )
  Occurrence <- data.frame(matrix(nrow=0, ncol=length(OccurrenceTerms), data=NA))
  colnames(Occurrence) <- OccurrenceTerms
  #------------------------------------------
  #### create Taxa table ####
  #------------------------------------------
  # Taxa $ id [ref: occurrence $ taxon]
  TaxaTerms <- c("id", "name", "TaxonRank", "taxonID", "parent_taxa")
  Taxa <- data.frame(matrix(nrow=0, ncol=length(TaxaTerms), data=NA))
  colnames(Taxa) <- TaxaTerms

  if(metaformat == "MIxS"){
    #NOTE: works for sequencing data; what about taxonomic data of hosts?
    taxaList <- dataQC.TaxonListFromData(coredata)
    Occurrence[1:nrow(coredata),] <- ""
    
    #fill in the Occurrence table
    Occurrence$id <- 1:nrow(coredata)
    Occurrence$occurrenceID <- as.character(coredata$polaaar_name)
    Occurrence$occurrence_class <- "MachineObservation"
    Occurrence$catalogNumber <- find.dataset(coredata, TermsList=c("original_name", "INSDC_SampleID", "genbank_accession_numbers"))
    Occurrence$otherCatalogNumber <- find.dataset(coredata, TermsList=c("INSDC_SampleID", "genbank_accession_numbers", "sra_run_number", "biosample"))
    if(exists("SeqNumbers")){
      Occurrence$associated_sequences <- SeqNumbers
    }
    
    ## Event $ occurrence
    Event$occurrence <- Occurrence$id
    ## Event $ occurrence_exists
    Event$occurrence_exists <- TRUE
    
  }else if(metaformat == "DwC" & nrow(occurdata)>0){
    taxaList <- dataQC.TaxonListFromData(occurdata)
    Occurrence[1:length(taxaList),] <- ""
    
    #fill in the Occurrence table
    Occurrence$id<- 1:nrow(occurdata)
    Occurrence$occurrence_notes <- find.dataset(occurdata, TermsList=c("occurrenceID"))
    Occurrence$occurrence_notes <- find.dataset(occurdata, TermsList=c("occurrenceRemarks"))
    Occurrence$occurrence_status <- find.dataset(occurdata, TermsList=c("occurrenceStatus"))
    Occurrence$occurrence_class <- find.dataset(occurdata, TermsList=c("basisOfRecord"))
    Occurrence$other_catalog_numbers <- find.dataset(occurdata, TermsList=c("otherCatalogNumber"))
    Occurrence$catalog_number <- find.dataset(occurdata, TermsList=c("catalogNumber"))
    Occurrence$date_identified <- find.dataset(occurdata, TermsList=c("dateIdentified"))
    Occurrence$recorded_by <- find.dataset(occurdata, TermsList=c("identifiedBy"))
    Occurrence$associated_sequences <- find.dataset(occurdata, TermsList=c("associatedSequences"))
    
    #make a table to link to occurrences to the events
    if("eventID" %in% colnames(occurdata)){
      Occurrence_to_Event <- data.frame(polID=Occurrence$id, occID=Occurrence$occurrenceID, 
                                        evID=occurdata$eventID , stringsAsFactors = FALSE)
      Occurrence_to_Event <- aggregate(Occurrence_to_Event$polID, by=list(Occurrence_to_Event$evID), paste, collapse=",")
      colnames(Occurrence_to_Event)<-c("evID", "polID")

      ## Event $ occurrence
      Event$occurrence <- unname(unlist(sapply(as.character(Event$sample_name), 
                                               FUN = function(x){
                                                 gsub(x,Occurrence_to_Event[Occurrence_to_Event$evID==x,]$polID,x)
                                               })))
      ## Event $ occurrence_exists
      Event$occurrence_exists <- unname(unlist(sapply(as.character(Event$occurrence), 
                                               FUN = function(x){
                                                 if(x!=""){
                                                   x <- TRUE
                                                 }else{x<-FALSE}
                                               })))
    }
  }
  
  # add the taxa to the Occurrence table
  taxaList <- dataQC.taxaNames(taxaList)
  Occurrence$taxon<-taxaList$scientificName
  # remove taxa that are "" => these are sampels with only environmental data
  Occurrence <- Occurrence[Occurrence$taxon!="", ]
  
  # make a small table with the unique taxa to collect all the info from WORMS or GBIF
  cat("Completing the taxonomic information (this might take a while)...\n")
  taxaRanks <- c("kingdom", "phylum", "class", "order", "family", "genus", "specificEpithet")
  taxaKey <- dataQC.completeTaxaNamesFromRegistery(taxaList$scientificName)
  taxaKey$taxaTabID<-""
  
  # fill in the Taxa table
  for(tk in 1:nrow(taxaKey)){
    for(tkcol in 1:7){
      if(taxaKey[tk,taxaRanks[tkcol]] != ""){
        if(!taxaKey[tk,taxaRanks[tkcol]] %in% Taxa$name){ #new name, so add to the table
          currentRow <- nrow(Taxa)+1
          Taxa[currentRow,]$id <- currentRow
          Taxa[currentRow,]$name <- taxaKey[tk,taxaRanks[tkcol]]
          Taxa[currentRow,]$TaxonRank <- taxaRanks[tkcol]
          Taxa[currentRow,]$taxonID <- ""
          if(tkcol==1){#kingdom level has no parent
            Taxa[currentRow,]$parent_taxa <- ""
          }else{
            tkmin<-tkcol
            repeat{
              if(length(Taxa[Taxa$name==taxaKey[tk,taxaRanks[tkmin-1]],]$id) > 0){
                Taxa[currentRow,]$parent_taxa  <- Taxa[Taxa$name==taxaKey[tk,taxaRanks[tkmin-1]],]$id
                break
              }else if(tkmin==1){
                break
              }else{
                tkmin <- tkmin-1
              }
            }
          }
        }
        taxaKey[tk,]$taxaTabID <- Taxa[Taxa$name==taxaKey[tk,taxaRanks[tkcol]],]$id #keep updating ID untill highest level
      }
    }
  }

  Occurrence$taxon <- unname(unlist(sapply(as.character(Occurrence$taxon), 
                                           FUN = function(x){
                                             gsub(x,taxaKey[taxaKey$scientificName==x,]$taxaTabID,x)
                                           })))
  
    
  #------------------------------------------
  #### create Environment table ####
  #### create unitsPol table ####
  #### create sampling_method table ####
  #### create variable table ####
  #------------------------------------------
  # this is the last set of tables to be checked, and will take most of the remaining variables that have not been recognized by the system
  # again: important to know we assume the data has gone through the process.metadata() function
  # Environment $ id [ref: Event $ environment]
  # Environment $ env_methods [ref: sampling_methods $ id]
  # Environment $ env_units [ref: units $ id]
  # Environment $ sequences [ref: Sequences $ id]
  # Environment $ env_variable [ref: Variable $ id]

  ## Environment
  EnvironmentTerms <- c("id", "env_sample_name", "created_at", "Latitude", "Longitude",
                        "link_climate_info", "env_variable", "env_method", "env_units", 
                        "sequences", "env_numeric_value", "env_text_value")
  Environment <- data.frame(matrix(nrow=0, ncol=length(EnvironmentTerms), data=""), stringsAsFactors = FALSE)
  colnames(Environment) <- EnvironmentTerms
  # units
  unitsTerms <- c("id", "name", "html_tag")
  unitsPol <- data.frame(matrix(nrow=0, ncol=length(unitsTerms), data=""), stringsAsFactors = FALSE)
  colnames(unitsPol) <- unitsTerms
  # Variable
  #var_type either TXT (text) or NUM (numeric)
  VariableTerms <- c("id", "name", "var_units", "method", "var_type")
  Variable <- data.frame(matrix(nrow=0, ncol=length(VariableTerms), data=""), stringsAsFactors = FALSE)
  colnames(Variable) <- VariableTerms
  # sampling_method 
  sampling_methodTerms <- c("id", "shortname", "description")
  sampling_method <- data.frame(matrix(nrow=0, ncol=length(sampling_methodTerms), data=""), stringsAsFactors = FALSE)
  colnames(sampling_method) <- sampling_methodTerms

  
  
  ### environmental variables
  ### one line per measurement
  ### refer back to Event $ environment
  # first remove any columns that are not environmental variables:
  env_metadata <- env_metadata[,!colnames(env_metadata) %in% c(as.character(TermsLib[TermsLib$polaaarDB_environment==0,]$name))]
  
  for(envVar in setdiff(colnames(env_metadata), c("polaaar_name", "polaaar_seqID"))){
    # data for the variable table
    env_Variable <- as.character(env_metadata[,colnames(env_metadata)==envVar])
    names(env_Variable) <- env_metadata$polaaar_name
    env_units <- metaunits[envVar]
    env_Variable <- env_Variable[!env_Variable %in% c(NA, "NA", "", "not collected", "not_collected")]
    
    #determine var_type
    if(any(grepl("[[:alpha:]]+", env_Variable))){
      var_type <- "TXT"
    }else{
      var_type <- "NUM"
    }
    
    Environment_sub <- data.frame(matrix(nrow=length(env_Variable), ncol=length(EnvironmentTerms), data=""), stringsAsFactors = FALSE)
    colnames(Environment_sub) <- EnvironmentTerms
    
    # Environment  $id
    Environment_sub$id <- (nrow(Environment)+1):(nrow(Environment)+length(env_Variable))
    # Environment $ env_sample_name
    Environment_sub$env_sample_name <- c(names(env_Variable))
    # Environment $ created_at
    Environment_sub$created_at <- rep(Sys.Date(), length(env_Variable))
    
    # latitude-longitude
    if("decimalLatitude" %in% colnames(coredata) &
       "decimalLongitude" %in% colnames(coredata)){
      Environment_sub$Latitude <- coredata[coredata$polaaar_name %in% names(env_Variable),]$decimalLatitude
      Environment_sub$Longitude <- coredata[coredata$polaaar_name %in% names(env_Variable),]$decimalLongitude
    }else if("lat_lon" %in% colnames(coredata)){
      lat<-unname(sapply(coredata[coredata$polaaar_name %in% names(env_Variable),]$lat_lon, FUN=function(x){x<-strsplit(x, " ")[[1]][1]}))
      lon<-unname(sapply(coredata[coredata$polaaar_name %in% names(env_Variable),]$lat_lon, FUN=function(x){x<-strsplit(x, " ")[[1]][2]}))
      Environment_sub$Latitude <- as.numeric(lat)
      Environment_sub$Longitude <- as.numeric(lon)
    }
    
    # Environment $ link_climate_info 
    if("climate_environment" %in% colnames(coredata)){
      Environment_sub$link_climate_info <- coredata[coredata$polaaar_name %in% names(env_Variable),]$climate_environment
    }
    
    if(var_type=="TXT"){
      Environment_sub$env_numeric_value <- env_Variable
    }else{
      Environment_sub$env_text_value <- env_Variable
    }
    
    Variable <- rbind(Variable,
                      data.frame(id=(nrow(Variable)+1),
                                 name= envVar,
                                 var_units="",
                                 method="",
                                 var_type=var_type
                      ))
    #Environment $ env_variable
    Environment_sub$env_variable <- rep((nrow(Variable)+1), length(env_Variable))
    
    if(!env_units %in% unitsPol$name && !is.na(env_units)){
      unitsPol <- rbind(unitsPol,
                        data.frame(id=(nrow(unitsPol)+1),
                                   name= env_units,
                                   html_tag=""
                        ))
      rownames(unitsPol)[nrow(unitsPol)]<- names(env_units)
    }
    #Environment $ env_units
    Environment_sub$env_units <- rep(unitsPol[unitsPol$name==env_units,]$id, length(env_Variable))
    Variable$var_units <- unitsPol[unitsPol$name==env_units,]$id
    
    #Environment $ env_method
    # ususlly methods associated with environmental measurements are not given
    
    #Environment $ sequences
    Environment_sub$sequences <- Environment_sub$env_sample_name
    Environment_sub$sequences <- unlist(sapply(Environment_sub$sequences, FUN=function(x){
      x <- coredata[coredata$polaaar_name==x,]$polaaar_SeqID
    }))
    Environment<-rbind(Environment, Environment_sub)
  }
  
  ## Event $ environment, environment_exists
  env_col <- aggregate(id ~ env_sample_name, Environment, paste0, collapse = ",")
  rownames(env_col)<- env_col$env_sample_name
  env_col<-combine.data.frame(env_col, Event, fill=NA, merge.cols=FALSE)
  env_col <- env_col[,c("sample_name", "id")]
  Event$environment <- env_col[order(match(env_col$sample_name,Event$sample_name)),]$id
  Event$environment_exists <- Event$environment
  Event$environment <- unlist(sapply(Event$environment, function(x){if(is.na(x)){x<-""};return(x)}))
  Event$environment_exists <- unlist(sapply(Event$environment_exists, function(x){if(is.na(x)){x<-FALSE}else{x<-TRUE};return(x)}))

  
  #------------------------------------------
  #### write out the tables ####
  #------------------------------------------
  write.csv(Event, paste(dest.dir, "/", baseName, "_Event.csv", sep=""), na="", row.names = FALSE)
  write.csv(ParentEvent, paste(dest.dir, "/", baseName, "_ParentEvent.csv", sep=""), na="", row.names = FALSE)
  write.csv(event_type, paste(dest.dir, "/", baseName, "_event_type.csv", sep=""), na="", row.names = FALSE)
  write.csv(Metadata, paste(dest.dir, "/", baseName, "_Metadata.csv", sep=""), na="", row.names = FALSE)
  write.csv(Sequences, paste(dest.dir, "/", baseName, "_Sequences.csv", sep=""), na="", row.names = FALSE)

  write.csv(Environment, paste(dest.dir, "/", baseName, "_Environment.csv", sep=""), na="", row.names = FALSE)
  write.csv(Variable, paste(dest.dir, "/", baseName, "_Variable.csv", sep=""), na="", row.names = FALSE)
  write.csv(unitsPol, paste(dest.dir, "/", baseName, "_units.csv", sep=""), na="", row.names = FALSE)
  
  ## conditional tables: these are not always made (depends on the input data)
  if(exists("occurrence")){
    write.csv(Sequences, paste(dest.dir, "/", baseName, "occurrence.csv", sep=""), na="", row.names = FALSE)
  }
  if(exists("Taxa")){
    write.csv(Sequences, paste(dest.dir, "/", baseName, "Taxa.csv", sep=""), na="", row.names = FALSE)
  }
  if(exists("sampling_method")){
    write.csv(Sequences, paste(dest.dir, "/", baseName, "sampling_method.csv", sep=""), na="", row.names = FALSE)
  }
  if(exists("Reference")){
    write.csv(Reference, paste(dest.dir, "/", baseName, "_Reference.csv", sep=""), na="", row.names = FALSE)
  }
  
  cat(paste("The data has been be written to CSV files with \"", baseName, "_\" as base name", sep=""))

}

  




