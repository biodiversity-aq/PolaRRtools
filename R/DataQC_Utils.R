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
## ideas to work out: do a QC on the Biome (use : as sep, envo vocab,...)
#    # assume words of one terms are separated by a " ", ".", and "_" => change this to space
#Metadata$env_biome <-gsub(".", " ", Metadata$env_biome, fixed=TRUE)
#Metadata$env_biome <-gsub("_", " ", Metadata$env_biome, fixed=TRUE)
#Metadata$env_biome<-gsub("\\s+", " ", Metadata$env_biome, fixed=FALSE)
# change possible hyrarchy separators to ":"
#Metadata$env_biome <-gsub(",", ":", Metadata$env_biome, fixed=TRUE)
#Metadata$env_biome<-gsub(">", ":", Metadata$env_biome, fixed=TRUE)
#Metadata$env_biome<-gsub(";", ":", Metadata$env_biome, fixed=TRUE)
#Metadata$env_biome<-gsub("|", ":", Metadata$env_biome, fixed=TRUE)
#Metadata$env_biome<-gsub(":+", ":", Metadata$env_biome, fixed=FALSE)


dataQC.dateCheck <- function(dataset, date.colnames=c("date", "Date", "collection_date")){
  #' @author Maxime Sweetlove CC-BY 4.0 2019
  #' @description dataQC.dateCheck looks in the columns of a dataset (dataframe) for a column with dates, 
  #' and transforms them to the YYYY-MM-DD format.
  #' @usage dataQC.dateCheck(dataset, date.colnames)
  #' @param dataset dataframe. The dataset where the date column should be found
  #' @param date.colnames character vector. a list of potential names for the column with the date. e.g. c("date", "Date", "collection date")
  #' @details The date column is found based on a user-provided list of possible names to look for (data.colnames argument)
  #' If a columnname is found that corresponds to a term in the list, the dates will be convered to the YYYY-MM-DD format, if the original format can be recognized.
  #' @return a list of length 2, with "$values" a dataframe with the same dimentions as the 
  #' dataset argument, and "$warningmessages" a vector with potential warning messages as character strings.
  #' @example dataQC.dateCheck(dataset, date.colnames=c("date", "Date"))
  
  warningmessages<-c()
  NAvals <- c("NA", "ND", "unnkown", "", NA, NULL, "na")
  xdate <- intersect(date.colnames, colnames(dataset))
  if(length(xdate)>1){
    warningmessages <- multi.warnings(paste("multiple date columns found, only executed QC on", xdate[1]), warningmessages)
    xdate <- xdate[1]
  }
  
  if(length(xdate)==1){
    date_values <- as.character(dataset[,xdate])
    # try to format date in YYYY-MM-DD format
    for(i in 1:length(date_values)){
      if(!is.na(date_values[i]) && !(gsub("/", "", date_values[i]) %in% NAvals) && !(gsub("-", "", date_values[i]) %in% NAvals)){
        date_values[i] <- gsub("/", "-", date_values[i], fixed=TRUE)
        date_values[i] <- gsub(" ", "", date_values[i], fixed=TRUE)
        date_values[i] <- gsub(".", "-", date_values[i], fixed=TRUE)
        date_values[i] <- gsub("_", "-", date_values[i], fixed=TRUE)
        if(grepl("[A-Za-z]", date_values[i])){
          date_values[i] <- gsub("jan.+-|jan-", "01-", tolower(date_values[i]))
          date_values[i] <- gsub("feb.+-|feb-", "02-", tolower(date_values[i]))
          date_values[i] <- gsub("mar.+-|mar-", "03-", tolower(date_values[i]))
          date_values[i] <- gsub("apr.+-|apr-", "04-", tolower(date_values[i]))
          date_values[i] <- gsub("may.+-|may-", "05-", tolower(date_values[i]))
          date_values[i] <- gsub("jun.+-|jun-", "06-", tolower(date_values[i]))
          date_values[i] <- gsub("jul.+-|jul-", "07-", tolower(date_values[i]))
          date_values[i] <- gsub("aug.+-|aug-", "08-", tolower(date_values[i]))
          date_values[i] <- gsub("sep.+-|sep-", "09-", tolower(date_values[i]))
          date_values[i] <- gsub("oct.+-|oct-", "10-", tolower(date_values[i]))
          date_values[i] <- gsub("nov.+-|nov-", "11-", tolower(date_values[i]))
          date_values[i] <- gsub("dec.+-|dec-", "12-", tolower(date_values[i]))
          #how it was written before: date_values[i] <- gsub("\\sdec.+\\s|\\sdec\\s", "-12-", tolower(date_values[i]))
        }
        date_values[i] <- strsplit(date_values[i], "T")[[1]][1]
        date_split <- strsplit(date_values[i], "-")
        
        if(length(date_split[[1]])==3 && nchar(date_split[[1]][3])==4){ #assume DD-MM-YYYY
          year <- as.character(date_split[[1]][3])
          mnth <- as.character(sprintf("%02d", as.numeric(date_split[[1]][2])))
          day <- as.character(sprintf("%02d", as.numeric(date_split[[1]][1])))
          date_values[i]<-paste(year, mnth, day, sep="-")
        }else if(length(date_split[[1]])==3 && nchar(date_split[[1]][1])==4){ #assume YYYY-MM-DD
          year <- as.character(date_split[[1]][1])
          mnth <- as.character(sprintf("%02d", as.numeric(date_split[[1]][2])))
          day <- as.character(sprintf("%02d", as.numeric(date_split[[1]][3])))
          date_values[i]<-paste(year, mnth, day, sep="-")
        }else if(length(date_split[[1]])==2 && nchar(date_split[[1]][2])==4){ #assume MM-YYYY
          year <- as.character(date_split[[1]][2])
          mnth <- as.character(sprintf("%02d", as.numeric(date_split[[1]][1])))
          date_values[i]<-paste(year, mnth, sep="-")
        }else if(length(date_split[[1]])==2 && nchar(date_split[[1]][1])==4){ #assume YYYY-MM
          year <- as.character(date_split[[1]][1])
          mnth <- as.character(sprintf("%02d", as.numeric(date_split[[1]][2])))
          date_values[i]<-paste(year, mnth, sep="-")
        }else if(!(length(date_split)==1 && nchar(date_split)==4)){ #assume YYYY, leave value as is.
          warningmessages <- multi.warnings("check the formats of the collection dates", warningmessages)
        }
      }else{
        warningmessages <- multi.warnings("some collection dates are unknown", warningmessages)
        date_values[i] <- ""
      }
    }
  }else{
    warningmessages <- multi.warnings("no collection dates found", warningmessages)
    date_values <- rep("", nrow(dataset))
  }
  return(list(values=date_values, warningmessages=warningmessages))
}

dataQC.LatitudeLongitudeCheck <- function(dataset, latlon.colnames=list(c("lat_lon"),c("latitude"), c("longitude"))){
  #' @author Maxime Sweetlove CC-BY 4.0 2019
  #' @description dataQC.LatitudeLongitudeCheck looks in the columns of a dataset (dataframe) for a column with coordinates (latitude/longitude in a single or double column system) and transforms them a standardized decimal format (see details).
  #' @usage dataQC.LatitudeLongitudeCheck(dataset, latlon.colnames)
  #' 
  #' @param dataset dataframe. The dataset where the date column should be found
  #' @param latlon.colnames a list of length 3 with character vectors. Three vectors of potential names for the columns with the latitude-longidude values. The first vector in the list are names if latitude and longitude were to be in the same column (e.g. the MIxS lat_lon format), the second and third are for when latitude and longitude, respectively, are in seperate columns. Example: list(c("lat_lon"), c("latitude"), c("longitude"))
  #' 
  #' @details The date column is found based on a user-provided list of possible names to look for (latlon.colnames argument).
  #' First, a single column is searched where latitude and longitude are noted in a single field, if this
  #' returns no result, latitude and longitude are looked for in seperate fields. When found, the coordinates
  #' are transformed to decimals and returned as a single field, with values separated by a single space 
  #' That is: (X Y), with X a numeric decimal latitude value and Y a numeric decimal longitude value.
  #' @seealso lat_lonSymbol.to.value and degree.to.decimal
  #' 
  #' @return a list of length 2, with "$values" a vector of same length as the number of rows in the
  #' dataset argument, and "$warningmessages" a vector with potential warning messages as character strings.
  
  #' @example dataQC.LatitudeLongitudeCheck(dataset, latlon.colnames=list(c("lat_lon"), c("latitude"), c("longitude")))
  
  warningmessages<-c()
  latlonName <- intersect(latlon.colnames[[1]], colnames(dataset))
  latName <- intersect(latlon.colnames[[2]], colnames(dataset))
  lonName <- intersect(latlon.colnames[[3]], colnames(dataset))
  lat_space_lon_output<-c()
  
  NAvals <- c("NA", "ND", "-", "/", "?", "unnkown", "")
  
  # check format lat_lon: 1 field
  if(length(latlonName)>=1){ #case 1: latitude and longitude are in the same field
    if(length(latlonName)>1){ #more than one column detected
      latlonName <- latlonName[1]
      warningmessages <- multi.warnings("multiple lat_lon columns found, just used the first one", warningmessages)
    }
    latlon_values <- dataset[,latlonName]
    for(i in 1:length(latlon_values)){
      if(!is.na(latlon_values[i]) && ! latlon_values[i] %in% NAvals){
        #try to find the right separator
        #first look for a tab
        latlon_val <- gsub('\t', " ", latlon_values[i])
        #ten split on for spaces
        latlon_val <- gsub("^\\s+|\\s+$", "", latlon_val) #remove leading and trailing spaces
        latlon_val <- strsplit(latlon_val," ")[[1]] #split on space
        if(length(latlon_val)==1){ #it's not separated by a space or a tab
          if(grepl(",", latlon_val)){ #comma can be decimal separator
            if(lengths(regmatches(latlon_val, gregexpr(",", latlon_val)))>1){
              # more than one comma: assume it is a decimal separator
              latlon_val <- gsub(",", ".", latlon_val)
            }else{
              #if there is just one comma, assume it is the value separator
              latlon_val <- gsub(",", " ", latlon_val)
            }
          }
          latlon_val <- gsub(",", " ", latlon_val)
          latlon_val <- gsub(";", " ", latlon_val)
          latlon_val <- gsub(":", " ", latlon_val)
          latlon_val <- gsub("|", " ", latlon_val)
          latlon_val <- gsub("\\s+", " ", latlon_val) #collapse multiple spaces
          latlon_val <- gsub("^\\s+|\\s+$", "", latlon_val) #remove leading and trailing spaces
          latlon_val <- strsplit(latlon_val," ")[[1]]
        }
        #weed out the NAs
        if(length(latlon_val)==1){#is the lat_lon is still not split, it must be NA
          if(is.na(latlon_val) | latlon_val %in% NAvals){
            lat_space_lon_output<-c(lat_space_lon_output, "")
            warningmessages <- multi.warnings("there are unknowns in the coordinates", warningmessages)
          }else{
            warningmessages <- multi.warnings("for some coordinates, the format could not be recognized", warningmessages)
            lat_space_lon_output<-c(lat_space_lon_output, latlon_val)
          }
        }else{
          #re-formatting the lat_lon values
          if(length(latlon_val)==2){
            lat<-latlon_val[1]
            lon<-latlon_val[2]
          }else if(length(latlon_val)==4){
            lat<-paste(latlon_val[1], latlon_val[2], sep="")
            lon<-paste(latlon_val[3], latlon_val[4], sep="")
          } 
          
          #converto to decimal format
          lat <- coordinate.to.decimal(lat)
          lon <- coordinate.to.decimal(lon)
          lat_space_lon_output<-c(lat_space_lon_output, paste(lat, lon))
        }
      }else{
        lat_space_lon_output<-c(lat_space_lon_output, "")
        warningmessages <- multi.warnings("there are unknowns in the coordinates", warningmessages)
      }
    }
    # check format lat_lon: 2 fields
  }else if(length(latName)>=1 & length(lonName)>=1){ # case2: latitude and longitude are be in seperate fields
    if(length(latName)>1 | length(lonName)>1){
      latName <- latName[1]
      lonName <- lonName[1]
      warningmessages <- multi.warnings("multiple latitude or longitude columns found, just used the first ones", warningmessages)
    }
    lat_values <- dataset[,latName]
    lon_values <- dataset[,lonName]
    if(length(lat_values)==length(lon_values)){
      for(i in 1:length(lat_values)){
        lat <- lat_values[i]
        lon <- lon_values[i]
        
        if(is.na(lat) | lat %in% NAvals | is.na(lon) | lon %in% NAvals){
          warningmessages <- multi.warnings("there are NAs in the coordinates", warningmessages)
          lat_space_lon_output<-c(lat_space_lon_output, "")
        }else{
          lat <- coordinate.to.decimal(lat)
          lon <- coordinate.to.decimal(lon)
          lat_space_lon_output<-c(lat_space_lon_output, paste(lat, lon))
        }
      }
    }else{
      warningmessages <- multi.warnings("the length of the latitude and longitude columns do not match, process not executed", warningmessages)
      lat_space_lon_output<-c(lat_space_lon_output, "")
    }
  }else{ #case 3: no latitude and longitude data found
    warningmessages <- multi.warnings("no latitude longitudes found", warningmessages)
    lat_space_lon_output<-""
  }
  return(list(values=lat_space_lon_output, warningmessages=warningmessages))
}

dataQC.guess.env_package.from.data <- function(dataset, pckge.colnames=c("env_package", "ScientificName")){
  #' @author Maxime Sweetlove CC-BY 4.0 2019
  #' @description dataQC.guess.env_package.from.data looks in the columns of a dataset for clues to what the most appropriate MIxS environmental could be for each sample.
  #' @usage dataQC.guess.env_package.from.data(dataset)
  #' 
  #' @param dataset dataframe. The dataset for which the MIxS environmental package should be found or guessed
  #' @param pckge.colnames a character vector. A vector with the potential names for the column where the environmental package can be found. Place the terms in order of decreasing likeliness.
  #' @details The GSC MIxS standard requires an appropriate environmental package with MIxS terms to be selected to document the data. Some data and nucleotide archives enforce their users to select such a package. This function is made to automatically either find the package in a dataset, of guess it based on the data that is present. 
  #' 
  #' @return a list of length 2, with "$values" a vector of same length as the number of rows in the dataset argument, and "$warningmessages" a vector with potential warning messages as character strings.
  
  #' @example dataQC.guess.env_package.from.data(df)
  
  warningmessages <- c()
  warningmessages<-multi.warnings("the env_package was not specified. An educated guess was made, but should be checked", warningmessages)
  
  #look if a collumn can be found
  pck <- intersect(pckge.colnames, colnames(dataset))
  if(length(pck)>0){
    env_package<-dataset[,pck[1]]
  } else{
    env_package<-NULL
    warningmessages<-multi.warnings("No env_package could be infered", warningmessages)
  }
  
  # try to convert the input to a correct package name (necessary when the package is not given directly and needs to be infered)
  if(!is.null(env_package)){
    env_package <- as.character(env_package)
    for(pk in 1:length(env_package)){
      pk_val <- tolower(env_package[pk])
      if(grepl("water", pk_val)|grepl("sea", pk_val)|grepl("ocean", pk_val)|
         grepl("pond", pk_val)|grepl("river", pk_val)|grepl("lake", pk_val)|
         grepl("aquatic", pk_val)|grepl("lagustrine", pk_val)|grepl("marine", pk_val)){
        env_package[pk]<-"water"
      }else if(grepl("soil", pk_val)|grepl("earth", pk_val)|grepl("sand", pk_val)|
               grepl("loam", pk_val)|grepl("clay", pk_val)|grepl("silt", pk_val)|
               grepl("peat", pk_val)|grepl("chalk", pk_val)|grepl(".[a-z]sol", pk_val)){
        env_package[pk]<-"soil"
      }else if(grepl("built_environment", pk_val)|grepl("built-environment", pk_val)|
               grepl("building", pk_val)|grepl("concrete", pk_val)|grepl("brick", pk_val)|
               grepl("cement", pk_val)|grepl("pavement", pk_val)|grepl("house", pk_val)|
               grepl("lobby", pk_val)|grepl("room", pk_val)){
        env_package[pk]<-"built_environment"
      }else if(grepl("air", pk_val)|grepl("wind", pk_val)){
        env_package[pk]<-"air"
      }else if(grepl("sediment", pk_val)|grepl("floor", pk_val)|grepl("bottom", pk_val)){
        env_package[pk]<-"sediment"
      }else if(grepl("microbial_mat_biofilm", pk_val)|grepl("microbial mat/biofilm", pk_val)|
               grepl("biofilm", pk_val)|grepl("microbial", pk_val)|grepl("mat", pk_val)){
        env_package[pk]<-"microbial_mat_biofilm"
      }else if(grepl("human_associated", pk_val)|grepl("human-associated", pk_val)){
        env_package[pk]<-"human_associated"
      }else if(grepl("human_gut", pk_val)|grepl("human-gut", pk_val)){
        env_package[pk]<-"human_gut"
      }else if(grepl("human_oral", pk_val)|grepl("human-oral", pk_val)){
        env_package[pk]<-"human_oral"
      }else if(grepl("human_skin", pk_val)|grepl("human-skin", pk_val)){
        env_package[pk]<-"human_skin"
      }else if(grepl("human_vaginal", pk_val)|grepl("human-vaginal", pk_val)){
        env_package[pk]<-"human_vaginal"
      }else if(grepl("wastewater_sludge", pk_val)|grepl("wastewater/sludge", pk_val)){
        env_package[pk]<-"wastewater_sludge"
      }else if(grepl("plant_associated", pk_val)|grepl("plant-associated", pk_val)){
        env_package[pk]<-"plant_associated"
      }else if(grepl("host_associated", pk_val)|grepl("host-associated", pk_val)|
               grepl("host", pk_val)|grepl("tissue", pk_val)|grepl("gut", pk_val)|
               grepl("skin", pk_val)|grepl("stomach", pk_val)|grepl("mouth", pk_val)|
               grepl("anus", pk_val)|grepl("ear", pk_val)|grepl("vagina", pk_val)|
               grepl("lungs", pk_val)|grepl("genital", pk_val)|grepl("penis", pk_val)|
               grepl("saliva", pk_val)|grepl("urine", pk_val)|grepl("feaces", pk_val)|
               grepl("oesophagus", pk_val)){
        env_package[pk]<-"host_associated"
      }else if(grepl("miscellaneous_natural_or_artificial_environment", pk_val)|
               grepl("miscellaneous natural or artificial environment", pk_val)|
               grepl("miscellaneous", pk_val)){
        env_package[pk]<-"miscellaneous_natural_or_artificial_environment"
      }else{
        env_package[pk]<-"miscellaneous_natural_or_artificial_environment"
        warningmessages<-multi.warnings("for some samples no env_package could be infered", warningmessages)
      }
    }
  }

  return(list(values=env_package, warningmessages=warningmessages))
}

dataQC.TermsCheck <- function(observed=NA, exp.standard="MIxS", exp.section=NA, fuzzy.match=TRUE,
                        out.type="full"){
  #' @author Maxime Sweetlove ccBY 4.0 2020
  #' @description checks a set of terms (e.g. columnnames) to a Standard, and flags inconsistencies, gives solutions if possible.
  #' 
  #' @usage Check.terms(observed, exp.standard="MIxS", exp.section="core")
  #' 
  #' @param observed character vector. The terms to be checked
  #' @param exp.standard character. The expected standard to which he terms should comply. Either MIxS (Minimum Information on any x sequence), DwC (DarwinCore) or INSDC (International Nucleotide Sequence Database Consortium)
  #' @param exp.section character. Optionally an specific section standard where the terms should come from. When exp.standard is MIxS, the allowed sections are: core,  air, built_environment, host_associated, human_associated, human_gut, human_oral, human_skin, human_vaginal, microbial_mat_biofilm, miscellaneous_natural_or_artificial_environment, plant_associated, sediment, soil, wastewater_sludge, water. When exp.standard is DwC, the allowed sections are: event, occurence, emof
  #' @param fuzzy.match logical. If TRUE, fuzzy matching will be done when no corresponding term is found in the Standard
  #' @param out.type character. The type of the output. Either "full" (the output is a list of three lists: terms_OK=the correct terms, terms_wrongWithSolution = wrong terms with a proposed solution, and terms_notFound = terms that had no match), "logical" (output is logical vector for exact matches or not) or "best_match" (output returns a vector with the best matching terms). default full
  
  
  require(RecordLinkage)
  
  MIxS_sections <- c("core",  "air", "built_environment", "host_associated",
                     "human_associated", "human_gut", "human_oral", "human_skin", 
                     "human_vaginal", "microbial_mat_biofilm", "sediment", "soil", 
                     "miscellaneous_natural_or_artificial_environment", "plant_associated", 
                     "wastewater_sludge", "water")
  DwC_sections <- c("event", "occurrence", "emof")
  
  #chage terms to lowercase: lowers number of possible differences
  names(observed) <- observed
  observed <- tolower(observed)
  
  # 1. check input
  if(!tolower(out.type) %in% c("full", "logical", "best_match")){
    stop("out.type should be either:\n\t \"full\": output will be a list of three lists, with:\n\t\t1) the correct terms\n\t\t2) wrong terms with a solution\n\t\t3) terms that had no match\n\t\"logical\": the output will be a TRUE/FALSE vector\n\t\"best_match\": the output will be the best matching term")
  }
  
  if(!tolower(exp.standard) %in% c("mixs", "dwc", "insdc")){
    stop("expected standard should be either \"MIxS\", \"DwC\" or \"INSDC\".")
  }
  if(tolower(exp.standard) == "mixs"){
    exp.standard<-"MIxS"
    if(!(is.na(exp.section)) && !tolower(exp.section) %in% MIxS_sections || length(exp.section)>1){
      stop(paste("expected section for MIxS should be NA or ONE of the following:\n\t", paste(MIxS_sections, collapse="\n\t"), sep=""))
    }
  }
  if(tolower(exp.standard) == "dwc"){
    exp.standard<-"DwC"
    if(!(is.na(exp.section)) && !tolower(exp.section) %in% DwC_sections || length(exp.section)>1){
      stop(paste("expected section for DwC should be NA or ONE of the following:\n\t", paste(DwC_sections, collapse="\n\t"), sep=""))
    }
  }
  if(tolower(exp.standard) == "insdc"){
    exp.standard<-"INSDC"
    exp.section <- NA
  }
  
  # 2. get the right synonym list
  if(is.na(exp.section)){
    synonymList <- sapply(as.character(TermsLib[TermsLib$name_origin %in% exp.standard,]$synonyms), function(x){strsplit(x, ";")})
    names(synonymList) <- TermsLib[TermsLib$name_origin %in% exp.standard,]$name 
  }else{
    exp.section <- gsub("^core$", "MIxS_core", exp.section, fixed=FALSE)
    exp.section <- gsub("^event$", "DwC_Event", exp.section, fixed=FALSE)
    exp.section <- gsub("^occurrence$", "DwC_Occurrence", exp.section, fixed=FALSE)
    exp.section <- gsub("^emof$", "DwC_eMoF", exp.section, fixed=FALSE)
    synonymList <- sapply(as.character(TermsLib[TermsLib[,exp.section] > 0,]$synonyms), function(x){strsplit(x, ";")})
    names(synonymList) <- TermsLib[TermsLib[,exp.section] > 0,]$name 
  }
  
  # 3. check terms
  terms_OK <- c()
  terms_wrongWithSolution <- c()
  terms_notFound <- c()
  for(obs_term in observed){
    obs_termName <- names(observed[grepl(paste("\\b", obs_term, "\\b", sep=""), observed, fixed=FALSE)])
    # first look for complete matches
    term_match <- synonymList[grepl(paste("\\b", obs_term, "\\b", sep=""), synonymList, fixed=FALSE)]
    if(length(term_match)==1){
      # one complete match
      if(obs_termName==names(term_match)){
        # terms is completely OK
        terms_OK <- c(terms_OK, obs_termName)
      }else{
        # term is one of the synonyms
        solutionNames <- names(terms_wrongWithSolution) #the names are the obsrved terms
        terms_wrongWithSolution <- c(terms_wrongWithSolution, names(term_match)) #the values are the matching term from the standard
        names(terms_wrongWithSolution) <- c(solutionNames, obs_termName)
      }
    }else{
      # no match found: check for spelling errors with a fuzzy match
      if(fuzzy.match){
        distances <- RecordLinkage::jarowinkler(obs_term, unlist(synonymList)) # caluculate all word distances
        if(max(distances)>0.85){ #treshold step, not bothering with worse matches
          dist_mtch <- match(max(distances),distances) #get index of best match
          if(length(dist_mtch)==1){ #only one match can be best, otherwise need human input
            dist_mtch <- unlist(synonymList)[dist_mtch] #get term that matched best
            term_match <- synonymList[grepl(paste("\\b", dist_mtch, "\\b", sep=""), synonymList, fixed=FALSE)]
            # return best solution
            solutionNames <- names(terms_wrongWithSolution)
            terms_wrongWithSolution <- c(terms_wrongWithSolution, names(term_match))
            names(terms_wrongWithSolution) <- c(solutionNames, obs_termName)
          }else{
            terms_notFound <- c(terms_notFound, obs_termName)
          }
        }else{
          terms_notFound <- c(terms_notFound, obs_termName)
        }
      }else{
        terms_notFound <- c(terms_notFound, obs_termName)
      }
    }
  }
  
  
  if(out.type=="full"){
    out<-list(terms_OK = terms_OK,
              terms_wrongWithSolution = terms_wrongWithSolution,
              terms_notFound = terms_notFound)
  }else if(out.type=="logical"){
    out <- sapply(observed, function(w){
      if(w %in% terms_OK){w<-TRUE}else{w<-FALSE}
    })
  }else if(out.type=="best_match"){
    out <- sapply(observed, function(w){
      if(w %in% terms_wrongWithSolution){w<-names(terms_wrongWithSolution[terms_wrongWithSolution %in% w])
      }else if(w %in% terms_notFound){w<-NA
      }else{w}})
  }
  
  
  return(out)
}

dataQC.generate.footprintWKT <- function(dataset, NA.val=""){
  #' @author Maxime Sweetlove ccBY 4.0 2020
  #' @description generate a footprintWKT column for a dataset (preferably DarwinCore eventCore). Lowest level events must be point locations.
  #' @param dataset data.frame. A data.frame with Event data, formatted as a DarwinCore EventCore file. Must include decimalLatitude, decimalLongitude, eventID and parentEventID.
  #' @param NA.val character. What the fill in when there is no data. Default ""
  #' @usage dataQC.generate.footprintWKT(dataset)
  #' @details This function will format the geographic coordinates of an event (in the fields decimalLatitude and decimalLongitude) into the WellKnownText format. For events that are at the lowest level of the hierarchy, it assumes a single point location. Higher level events are desribed by polygons based on the coordinates of the child events.
  require(mapview)
  if(!"decimalLatitude" %in% colnames(dataset) | 
     !"decimalLongitude" %in% colnames(dataset) |
     !"eventID" %in% colnames(dataset)){
    stop("The data must include the following columns:\n\tdecimalLatitude, decimalLongitude, eventID")
  }
  footprintWKT_vec <- c()
  if("parentEventID" %in% colnames(dataset)){
    for(ev in dataset$eventID){
      if(!ev %in% dataset$parentEventID){
        # in this case the event is at the lowest level, assume it is a point
        lat<-dataset[dataset$eventID==ev,]$decimalLatitude
        lon<-dataset[dataset$eventID==ev,]$decimalLongitude
        if(!is.na(lat) & !is.na(lon) & lat!="" & lon!=""){
          footprintWKT_vec <- c(footprintWKT_vec, paste("POINT (", as.character(lat), " ", as.character(lon), ")", sep=""))
        } else{
          footprintWKT_vec <- c(footprintWKT_vec, NA.val)
        }
      }else{
        # this event is the parent of one or more child events. we need to find all the coordinates of the lowest events under this parent
        child_ev <- dataset[dataset$parentEventID == ev,]$eventID #the direct descendants
        doneList<-c()
        n_added <- length(child_ev)
        while(! n_added == 0){ #loop through the events to add all the other descendants
          n_added = 0
          for(d in child_ev){
            if(!d %in% doneList){
              d2 <- dataset[dataset$parentEventID == d,]$eventID
              doneList <- c(doneList,d)
              if(length(d2)>0){
                child_ev<- c(child_ev, d2)
                n_added <- n_added +1
              }
            }
          }
        }
        child_ev <- setdiff(child_ev, unique(dataset$parentEventID))
        
        lat<-dataset[dataset$eventID %in% child_ev,]$decimalLatitude
        lon<-dataset[dataset$eventID %in% child_ev,]$decimalLongitude
        if(!all(is.na(lat)) & !all(is.na(lon)) & all(lat!="") & all(lon!="")){
          coords <- as.matrix(data.frame(x=as.numeric(lat), y=as.numeric(lon)))
          polygx <- mapview::coords2Polygons(coords, hole=FALSE, ID=ev)
          footprintWKT_vec <- c(footprintWKT_vec, writeWKT(polygx))
        } else{
          footprintWKT_vec <- c(footprintWKT_vec, NA.val)
        }
      }
    }
  }else{
    footprintWKT_vec <- paste("POINT (", as.character(dataset$decimalLatitude),
                              ", ", as.character(dataset$decimalLongitude), ")", 
                              sep="")
  }
  
  footprintWKT_vec <- gsub("POINT (NA NA)", NA.val, footprintWKT_vec, fixed=TRUE)
  return(footprintWKT_vec)
  
}


dataQC.eventStructure <- function(dataset, check.parentEventID=TRUE, 
                                  complete.hierarchy=FALSE, ask.input=TRUE){
  #' @author Maxime Sweetlove ccBY 4.0 2020
  #' @description checks if an eventID is present in a (QC'd) dataset, and generates one if not. Does the same of parentEventID if check.parentEventID is TRUE, and check the hierarchical relationships between eventID and parentEventID
  #' @param dataset data.frame. A data.frame
  #' @param check.parentEventID logical. If TRUE, the function will check for the presence of a parentEventID, and generate one if needed. The hierarchy between the eventIDs and the parentEventIDs is then also checked. Default TRUE.
  #' @param complete.hierarchy logical. If TRUE, the event-parentEvent hierarchy structure will be completed upt to a root. Any parentEvents that are not listed among the eventIDs will also be created. Default FALSE.
  #' @param ask.input logical. If TRUE, console input will be requested to the user when a problem occurs. Default TRUE
  
  #' @usage dataQC.eventStructure(dataset, check.parentEventID=TRUE, complete.hierarchy=FALSE, ask.input=TRUE)
  
  require(stringr)
  eventIDs<-c()
  parentEventIDs<-c()
  
  eventTable <- data.frame(eventID=rep(NA, nrow(dataset)), 
                           original_event=rep(TRUE, nrow(dataset)),
                           project=rep(NA, nrow(dataset)))
  rownames(eventTable)<-rownames(dataset)
  
  if(complete.hierarchy & !check.parentEventID){
    warning("complete.hierarchy is meaningless without checking the parentEventID column.")
  }
  
  # 1. check for the presence of eventID, if not, try to generate one
  if("eventID" %in% colnames(dataset)){
    eventTable$eventID <- dataset$eventID
    # QC the eventIDs: they have to be unique
    if(length(unique(eventTable$eventID)) != nrow(dataset)){
      if(dealWithErrors){
        # make the eventIDs unique
        eventTable$eventID <- make.names(eventTable$eventID, unique=TRUE)
      }else{
        stop("Some eventIDs were not unique.")
      }
    }
  }else{ ## there is no eventID
    # try to generate an eventID from sample names (but these have to be unique with non-missing values)
    for(attempt in c("original_name", "INSDC_SampleID", "sra_run_number", 
                     "genbank_accession_numbers")){
      if(attempt %in% colnames(dataset)){
        if(length(unique(dataset[,attempt]))==nrow(dataset) & 
           !any(is.na(dataset[,attempt])) & 
           !any(dataset[,attempt]=="")){
          eventTable$eventID <- dataset[,attempt]
          if(attempt!="original_name"){
            eventTable$eventID <- paste("INSDC:", eventTable$eventID, sep="")
          }
          break
        }
      }
    }
    if(all(is.na(eventTable$eventID))){
      # if no decent eventID could be found, generate one based on the project name
      for(attempt in c("biosample", "project_name","bioproject")){
        if(attempt %in% colnames(dataset)){
          eventTable$eventID <- dataset[,attempt]
          eventTable$eventID <- paste(eventTable$eventID, ":E",
                            stringr::str_pad(1:nrow(dataset), nchar(nrow(dataset)), pad = "0"), 
                            sep="")
          eventTable$eventID<-gsub("^:","", eventTable$eventID)
          eventTable$eventID<-gsub(" ","_", eventTable$eventID)
          eventTable$eventID<-gsub("-","_", eventTable$eventID)
          break
          }
        }
    }
    if(all(is.na(eventTable$eventID))){
      # still no eventID has been generated, just give it a number
      if(ask.input){
        cat("No eventIDs were found.\n\tPlease provide an eventID prefix, or type n to ignore.") 
        doNext <- readline() 
        if(!tolower(doNext) %in% c("n")){
          eventTable$eventID <- paste(doNext, ":E",
                            stringr::str_pad(1:nrow(dataset), nchar(nrow(dataset)), pad = "0"), 
                            sep="")
        }else{
          eventTable$eventID <- paste("unnamed_event:E",
                            stringr::str_pad(1:nrow(dataset), nchar(nrow(dataset)), pad = "0"), 
                            sep="")
        }
      }
    }
  }
  
  # 2. check for the presence of parentEventID, if not present, generate one
  if(check.parentEventID){
    if("parentEventID" %in% colnames(dataset)){
      eventTable$parentEventID <- dataset$parentEventID
    }else{
      # generate a parentEventID based on the project name
      if("project_name" %in% colnames(dataset)){
        eventTable$parentEventID <- dataset$project_name
        eventTable$project <- dataset$project_name
        eventTable$parentEventID<-gsub(" ","_", eventTable$parentEventID)
        eventTable$parentEventID<-gsub("-","_", eventTable$parentEventID)
      }
    }
  }
  
  if(all(is.na(eventTable$parentEventID))){
    if(ask.input){
      cat("No parentEventID was found.\n\tPlease provide an one that is equivalent to the scientific project in which the data was generated, or type n to ignore.") 
      doNext <- readline() 
      if(!tolower(doNext) %in% c("n")){
        eventTable$parentEventID<-rep(doNext, nrow(dataset))
        eventTable$project <- doNext
      }
    }else{
      eventTable$parentEventID<-rep("unnamed_project", nrow(dataset))
      eventTable$project <- "unnamed_project"
    }
  }
  
  #try to find a project name, or ask for one
  if(all(is.na(eventTable$project))){
    if("project_name" %in% colnames(dataset)){
      eventTable$project <- dataset$project_name
    }else if(ask.input){
      cat("No project name was found.\n\tPlease provide a project name, or type n to ignore.") 
      doNext <- readline() 
      if(!tolower(doNext) %in% c("n")){
        eventTable$project <- doNext
      }else{
        eventTable$project <- "unnamed_project"
      }
    }
  }
  
  #if there is just one project: this will be the root project
  #if there are multiple projects: create a root project
  if(complete.hierarchy){
    if(length(unique(eventTable$project)==1) && !unique(eventTable$project)==""){
      PRJ <- unique(eventTable$project)
    }else if(ask.input){
      if("" %in% unique(eventTable$project) & length(unique(eventTable$project))>1){
        cat("Multiple projects found or some sample had no project.\n\tPlease provide an overarching project name that can be used to unite all samples, or type n to ignore.") 
      }else if("" %in% unique(eventTable$project) & length(unique(eventTable$project))==1){
        cat("No project found.\n\tPlease provide an overarching project name that can be used to unite all samples, or type n to ignore.") 
      }else{
        cat("Multiple projects found.\n\tPlease provide an overarching project name that can be used to unite all samples, or type n to ignore.") 
      }
      doNext <- readline() 
      if(!tolower(doNext) %in% c("n")){
        PRJ <- doNext
      }else{
        PRJ <- "unnamed_project"
      }
    }
    
    ## complete the event hierarchy
    # first
    
    # then the overarching project
    PRJ
    
    
    # 6. include higer level events
    # add any parentEventIDs that are missing from the EventIDs
    if(length(missing_parents)>0){
      for(mp in missing_parents){
        dataset[nrow(dataset)+1,]<-NA
        rownames(dataset)[nrow(dataset)] <- mp
        dataset[nrow(dataset),]$eventID <- mp
        dataset[nrow(dataset),]$parentEventID <- ""
      }
    }
    
    # use the project as a root to all EventIDs and parentEventIDs
    if(length(missing_project)>0){
      if(length(PRJ)>1){
        PRJ <- "unnamed_project"
      }
      dataset[dataset$parentEventID=="","parentEventID"] <- PRJ
      dataset[nrow(dataset)+1,]<-NA
      rownames(dataset)[nrow(dataset)] <- "project"
      dataset[nrow(dataset),]$eventID <- PRJ
      dataset[nrow(dataset),]$parentEventID <- ""
    }
    
    

    
    
  }
  
  
  # 5. check if the EventIDs match the parentEventIDs
  if(check.parentEventID){
    missing_project<-c()
    missing_parents<-c()
    u_ev <- unique(dataset$eventID)
    u_pev <- unique(dataset$parentEventID)
    if(!"" %in% u_pev){
      missing_project<-c("") #the allows a project-level parent event to be added
    }else{
      #check if eventIDs with no parentEvent is also the parent of an event
      if(!all(dataset[dataset$parentEventID=="",]$eventID %in% dataset$parentEventID)){
        # if not, a project-level parent can be added
        missing_project<-c("")
      }else{
        if(length(dataset[dataset$parentEventID=="",]$eventID)>1){
          missing_project<-c("") # multiple high-ranking parents can be united by adding project-level parent
        }
      }
    }
    if(!all(u_pev[u_pev!=""] %in% u_ev)){
      missing_parents <- c(setdiff(u_pev, u_ev))
    }
  }
  

  

  
  if(check.parentEventID){
    data.out <- dataset[,colnames(dataset) %in% c("eventID", "parentEventID")]
  }else{
    data.out <- dataset[,colnames(dataset) %in% c("eventID")]
  }
  return(data.out)
}






### still need to finish
dataQC.measurement.to.meter <- function(dataset, meter.colnames=c("depth", "elev", "alt_elev", "filter_size", "tot_depth_water_col")){
  #' @author Maxime Sweetlove CC-BY 4.0 2019
  #' @description dataQC.measurement.to.meter looks in the columns of a dataset for given size/depth/length related measurements, and converts cells with value_unit format into a single value, saving the unit in a separate vector. some less used units are convrted into their nearest standard (See details).
  #' @usage dataQC.measurement.to.meter(dataset, meter.colnames=c("depth", "elev", "alt_elev", "filter_size", "tot_depth_water_col"))
  #' 
  #' @param dataset dataframe. The dataset where the measurement columns should be found
  #' @param meter.colnames a character vector. A vector with the names of columns names to target
  #' @details This function is mainly aimed at formatting cells, and removing the unit from the cell. If multiple decimal indicators (e.g. nano-, centi-, kilo-,...) of a unit are found, everything is converted to the smalest decimal indicator. Any non-metric units are converted to meter. If multiple units are found, nothing will be done.
  #' 
  #' @return a list of length 3, with "$values" a dataframe of same length as the number of rows in the dataset argument and the same number of columns as columns that were chacked, "$units" a named character vector with the unit associated with each column, and "$warningmessages" a vector with potential warning messages as character strings.
  #' @example dataQC.measurement.to.meter(df, c("depth"))

  QC_units<-c()
  items_shared <- intersect(meter.colnames, colnames(dataset))

  if(length(items_shared)>0){
    alternative_units<-data.frame(name=items_shared, unit_full=rep(NA, length(items_shared)))
    possible_units <- c("nm", "um", "??m", "mm", "cm", "dm", "m", "km")
    for(name_unit in items_shared){
      vals <- as.vector(as.character(dataset[,colnames(dataset) %in% name_unit]))
      names(vals)<-row.names(dataset)
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
  
  
}

### ideas for functions:
#dataQC.footprintWKT <- function(){
#to check if footprintWKT in DwC is correctly formatted
#POINT, LINESTRING, POLYGON, MULTIPOINT
#}

