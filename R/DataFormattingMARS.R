MarsLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/MarsLibrary.csv", 
                  header=TRUE)

write.metadata.MIxS.as.mars <- function(metadata.object, name.prefix=NULL, dest.dir=getwd(),
                                        add.missing.data.columns=TRUE, ask.input=TRUE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description write.metadata.MIxS.as.mars takes a metadata.MIxS object and outputs it as a MiMARKS and SequenceSet table to upload to mARS
  #' @param metadata.object a metadata.MIxS class object. The object to be submitted to mARS
  #' @param name a character string naming a file. Either a full path or a file name to write to the working directory. Invalid input will cause output to be printed to the console.
  #' @param checklist_accession character. The name of a MIxS environmetal package or it's ENA checklist accession number.
  #' @details ENA uses its own variant of MIxS, and has
  #' write.metadata.MIxS.as.ENA assumes the dataset has already been subjected to process the process.metadata function to standardize in input. If this is not the case, "garbage in, garbage out" is applicable. 
  #' @return tab separated .txt file
  #' @example 
  #' 
  
  warningmessages<-c()
  
  # 1. check input
  if(check.valid.metadata.MIxS(metadata.object)){
    metaunits <- metadata.object@units
    metasection <- metadata.object@section
    metapackage <- metadata.object@env_package
    metadata <- metadata.object@data
  }else{
    stop("metadata.object argument must be a metadata.MIxS class object.
         Common dataframes should be converted to metadata.MIxS with 
         the data.frame.to.metadata.MIxS function.")
  } 
  # 2. check output destination
  if(!is.character(name.prefix) | c(NULL,NA) %in% name.prefix | length(name.prefix)>1){
    stop("Need input for the name.prefix argument.")
  } 
  if(!dir.exists(file.path(dest.dir))){
    stop("Could not find the directory provided in the dest.dir argument.")
  }
  mimarks.name <- paste(dest.dir, "/MiMARKS_", name.prefix, ".csv", sep="")
  seqset.name <- paste(dest.dir, "/SeqSet_", name.prefix, ".csv", sep="")
  
  # 3. create MiMARKS tabble
  # 3.1 first row
  Row01 <- paste("section", "Structured Comment Name", "units", 
                 paste(rownames(metadata), collapse=","), sep=",")
  # 3.2 general terms
  mimarksTerms <- as.character(MarsLib[MarsLib$mars_mimarks==TRUE,]$name) #the basic terms that need to be in every dataset
  #associated units
  mimarksUnits <- as.character(TermsLib[match(setdiff(mimarksTerms, names(metaunits)), TermsLib$name),]$expected_unit)
  if(length(mimarksUnits)!=0){
    names(mimarksUnits) <- setdiff(mimarksTerms, names(metaunits))
  }
  mimarksUnits <- c(mimarksUnits, metaunits)
  #assocated MIxS sections
  mimarksSection <- as.character(TermsLib[match(setdiff(mimarksTerms, names(metasection)), TermsLib$name),]$section)
  if(length(mimarksSection)!=0){
    names(mimarksSection) <- setdiff(mimarksTerms, names(metasection))
  }
  mimarksSection <- c(mimarksSection, metasection)
  Row02<-character()
  #some additional steps for the primer field (to split forward and reverse primer in the seqset file)
  investigation_type<-character() 
  fw_primerName <- character() 
  fw_primerSeq <- character() 
  rv_primerName <- character()
  rv_primerSeq <- character()
  mimarks_data <- metadata # save data for the seqset round
  for(tx in 1:length(mimarksTerms)){
    Row0x_section <- as.character(mimarksSection[names(mimarksSection)==mimarksTerms[tx]])
    Row0x_section <- paste("MiMARKS_", Row0x_section,sep="")
    Row0x_name <- mimarksTerms[tx]
    Row0x_units <- as.character(mimarksUnits[names(mimarksUnits)==mimarksTerms[tx]])
    if(mimarksTerms[tx] %in% colnames(mimarks_data)){
      Row0x_data <- mimarks_data[,colnames(mimarks_data)==mimarksTerms[tx]]
      #remove the used data to see what columns will remain at the end
      mimarksUnits <- mimarksUnits[!names(mimarksUnits)==mimarksTerms[tx]] 
      mimarksSection <- mimarksSection[!names(mimarksSection)==mimarksTerms[tx]] 
      mimarks_data <- mimarks_data[,!colnames(mimarks_data)==mimarksTerms[tx],drop=FALSE] 
      if(mimarksTerms[tx]=="primer"){
        Row0x_data
        fw_primerSeq <- character() 
        rv_primerSeq <- character()
      }
    } else if(ask.input){
      if(mimarksTerms[tx] == "primer" && "mimarks-survey" %in% investigation_type){ ###special case
        cat(paste("no primer information found:\n\tPlease provide the forward primer name, or hit enter to leave blank.\n", sep=""))
        fw_primerName <- readline() 
        cat(paste("no primer information found:\n\tPlease provide the forward primer sequence, or hit enter to leave blank.\n", sep=""))
        fw_primerSeq <- readline() 
        cat(paste("no primer information found:\n\tPlease provide the reverse primer name, or hit enter to leave blank.\n", sep=""))
        rv_primerName <- readline() 
        cat(paste("no primer information found:\n\tPlease provide the reverse primer sequence, or hit enter to leave blank.\n", sep=""))
        rv_primerSeq <- readline() 
        Row0x_data <- rep(paste("(", fw_primerName, "):",fw_primerSeq, " (", rv_primerName, "):",rv_primerSeq, sep=""), nrow(mimarks_data))
      } else{
        cat(paste("no data found for ",mimarksTerms[tx],":\n\tPlease type the info to fill the cells, or hit enter to leave blank.\n", sep=""))
        user_input <- readline() 
        Row0x_data <- rep(user_input, nrow(mimarks_data))
      }
    }else{
      Row0x_data <- rep("", nrow(mimarks_data))
    }
    Row0x_data[is.na(Row0x_data)]<-""
    Row0x_data <- gsub("^NA$", "", Row0x_data)
    if(Row0x_name =="investigation_type"){
      investigation_type<-unique(Row0x_data)
    }
    Row0x_data <- paste(Row0x_data, collapse=",")
    Row0x <- paste(Row0x_section, Row0x_name, Row0x_units, Row0x_data, sep=",")
    Row02 <-  paste(Row02, Row0x, sep="\n")
    
  }
  
  # 3.3 remaining dataset-specific terms
  if(ncol(mimarks_data)>0){
    Row03 <- character()
    for(cx in 1:ncol(mimarks_data)){
      #assuming the data has gone through the QC of process.metadata()
      Row0x_section <- mimarksSection[colnames(mimarks_data[cx])]
      if(!Row0x_section %in% c("package", "miscellaneous")){
        Row0x_section <- paste("MiMARKS_", Row0x_section,sep="")
      }
      Row0x_name <- colnames(mimarks_data[cx])
      Row0x_units <- mimarksUnits[colnames(mimarks_data[cx])]
      Row0x_data <- mimarks_data[,cx]
      Row0x_data[is.na(Row0x_data)]<-""
      Row0x_data <- gsub("^NA$", "", Row0x_data)
      Row0x_data <- paste(Row0x_data, collapse=",")
      Row0x <- paste(Row0x_section, Row0x_name, Row0x_units, Row0x_data, sep=",")
      Row03 <- paste(Row03, Row0x, sep="\n")
    }
  }
  mimarks.data <- paste(Row01, Row02, Row03, sep="")
  
  # 4. create SeqSet tabble
  # 4.1 first row
  Row01 <- paste("unique_sequence_set_id", 
                 paste(rownames(metadata), collapse=","), sep=",")
  # 4.2 other rows
  seqsetTerms <- as.character(MarsLib[MarsLib$mars_seqset==TRUE,]$name)
  Row02<-character()
  for(tx in 1:length(seqsetTerms)){
    Row0x_name <- seqsetTerms[tx]
    if(seqsetTerms[tx] %in% colnames(metadata)){
      Row0x_data <- metadata[,colnames(metadata)==seqsetTerms[tx]]
    } else{
      mimarksEquiv <- as.character(MarsLib[MarsLib$name==seqsetTerms[tx],]$mimarks_equivalent)[1]
      if(mimarksEquiv == "investigation_type"){ ###special case
        Row0x_data <- gsub("mimarks-survey", "marker gene", Row0x_data)
      }
      if(nchar(mimarksEquiv)==0){ #case it had no mimarks equivalent or wasn't in the colnames of metadata
        if(ask.input){
          cat(paste("no data found for ",seqsetTerms[tx],":\n\tPlease type the info to fill the cells, or hit enter to leave blank.\n", sep=""))
          user_input <- readline() 
          if(user_input == ""){
            Row0x_data <- rep("", nrow(metadata))
          } else{
            Row0x_data <- rep(user_input, nrow(metadata))
          }
        }else{
          Row0x_data <- rep("", nrow(metadata))
        }
      }else{
        Row0x_data <- metadata[,colnames(metadata)==mimarksEquiv]
      }
    }
    Row0x_data[is.na(Row0x_data)]<-""
    Row0x_data <- gsub("^NA$", "", Row0x_data)
    Row0x_data <- paste(Row0x_data, collapse=",")
    Row0x <- paste(Row0x_name, Row0x_data, sep=",")
    Row02 <-  paste(Row02, Row0x, sep="\n")
  }
  seqset.data <- paste(Row01, Row02, sep="")
  
  # 5. finalize
  cat(paste("The data had been written to ", dest.dir, "\n",sep=""))
  
  write.table(mimarks.data, file=mimarks.name, 
              col.names = FALSE, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  write.table(seqset.data, file=seqset.name, 
              col.names = FALSE, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  }

