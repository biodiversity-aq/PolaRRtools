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
### general use

multi.warnings <- function(message_text, warningmessages){
  #' @usage multi.warnings(message_text, warningmessages)
  #' @param message_text a character string. The error message to be added
  #' @param warningmessages a vector with one ore more character strings. The previous error messages.
  if(! message_text %in% warningmessages){
    warningmessages <- c(warningmessages, message_text)
  }
  return(warningmessages)
}

#==============================================================
### Data table formatting ultis

combine.data.frame <- function(df1, df2, fill=NA, merge.cols=TRUE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description combine.data.frame merges two dataframes, optionally merging common columns.
  #' 
  #' @usage combine.data.frame(df1, df2, fill=NA, merge.cols=TRUE)
  #' 
  #' @param df1 a dataframe
  #' @param df2 a dataframe
  #' @param fill character or NA. A value to put into the cells that have no data.
  #' @param merge.cols boolean. Merge colomns with common name. default TRUE
  #' 
  #' @details Columns with matching names can or cannot be merged, rows are automatically
  #' bound (a wrapper of rbind), not merged.
  #' Any missing data as a result of the non-matchng columns will be filled by the fill argument.
  #' @example
  #' > df1 <- data.frame(col1=c(9,7,5), col2=c("ts",4,3))
  #' > df2 <- data.frame(col2=c(8,6), col3=c(8,7))
  #' > combine.data.frame(df1, df2, fill=NA, merge.cols=TRUE, merge.rows=FALSE)
  #'     col1 col2 rowNames col3
  #'1     9    ts        1   NA
  #'2     7    4        2   NA
  #'3     5    3        3   NA
  #'11   NA    8        1    8
  #'21   NA    6        2    7
  
  if(!class(df1) %in% c("data.frame", "matrix") | 
     !class(df2) %in% c("data.frame", "matrix") |
     ncol(df1) == 0 | nrow(df1) == 0 | 
     ncol(df2) == 0 | nrow(df2) == 0){
    stop("invalid input, must be a dataframe of matrix with at leaste 1 column and row.")
  }
  for(nc in 1:ncol(df1)){ #merging columns with factors is always trouble, so convert to character
    if(class(df1[,nc])=="factor"){
      df1[,nc] <- as.character(df1[,nc])
    }
  }
  for(nc in 1:ncol(df2)){ #idem
    if(class(df2[,nc])=="factor"){
      df2[,nc] <- as.character(df2[,nc])
    }
  }
  shared_cols <- intersect(colnames(df1), colnames(df2))
  if(merge.cols){
    df11 <- df1[,!colnames(df1) %in% shared_cols, drop=FALSE]
    df22 <- df2[,!colnames(df2) %in% shared_cols, drop=FALSE]
  }else{df22 <- df2; df11 <- df1}
  shared_rows <- intersect(rownames(df1), rownames(df2))
  
  if(length(shared_rows)>0){ #new column with original rownames
    df1$original_rowName <- row.names(df1)
    df2$original_rowName <- row.names(df2)
    shared_cols<-c(shared_cols, "original_rowName")
  }
  
  if(!ncol(df22)==0){ #if ncol(df22)==0, then df2 has no other colnames than df1
    df_UpRight <- data.frame(matrix(nrow = nrow(df1), ncol = ncol(df22), data=fill))
    colnames(df_UpRight) <- colnames(df22)
    rownames(df_UpRight) <- rownames(df1)
    df1_b <- cbind(df1, df_UpRight)
  }else{
    df1_b <- df1
  }
  
  if(!ncol(df11)==0){
    df11 <- df1[,!colnames(df1) %in% shared_cols, drop=FALSE]
    df_DownLeft <- data.frame(matrix(nrow = nrow(df2), ncol = ncol(df11), data=fill))
    rownames(df_DownLeft) <- rownames(df2)
    colnames(df_DownLeft) <- colnames(df11)
    df2_b <- cbind(df_DownLeft, df2)
  }else{
    df2_b <- df2
  }
  
  df_out <- data.frame(rbind(df1_b, df2_b))
  colnames(df_out)<-colnames(df1_b)
  
  if(merge.cols){
    df_out[c((nrow(df1_b)+1):(nrow(df_out))),shared_cols]<-df2[,shared_cols, drop=FALSE]
  }else{
    if(length(shared_rows)>0){
      df_out[c((nrow(df1_b)+1):(nrow(df_out))),"rowNames"]<-df2[,"rowNames", drop=FALSE]
    }
  }
  
  return(df_out)
}

combine.data <- function(d1, d2, fill=NA, variables.as.cols=TRUE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description combine.data combines two metadata.MIxS objects into one, merging all common variables (columns as default).
  #' 
  #' @usage combine.data(d1, d2, fill=NA, variables.as.cols=TRUE)
  #' 
  #' @param d1 a metadata.MIxS object
  #' @param d2 a metadata.MIxS object
  #' @param fill character or NA. A value to put into the cells that have no data.
  #' @param variables.as.cols boolean. If TRUE, the input data is assumed to have rows as samples and variables/parameters as columns. If FALSE the data is formatted the other was around. default TRUE
  #' 
  #' @details Variables with matching names are merged, if variables.as.cols=TRUE this means columns are merged,
  #' if FALSE rows are merged
  #' Any missing data as a result of the non-matchng variable-name will be filled by the fill argument.
  #' @example
  
  
  if(class(d1) %in% c("data.frame", "matrix") & 
     class(d2) %in% c("data.frame", "matrix")){
    if(ncol(d1) == 0 | nrow(d1) == 0 | 
       ncol(d2) == 0 | nrow(d2) == 0){
      stop("invalid input, each dataset must at leaste have 1 column and row.")
    }
    class_out <- "data.frame"
    if(variables.as.cols){
      df1<-d1
      df2<-d2
    } else{
      df1<-data.frame(t(d1))
      df2<-data.frame(t(d2)) 
    }
  } else if(check.valid.metadata.MIxS(d1) & check.valid.metadata.MIxS(d2)){
    class_out <- "metadata.MIxS"
    variables.as.cols<-TRUE
    df1<-d1@data
    df2<-d2@data
  } else{
    stop("invalid input, both datasets must be a dataframe, matrix or metadata.MIxS class object.
         possible causes of this error include:
         - a dataframe/matrix or metadata.MIxS@data object had 0 columns or rows. 
         - the two input dataset differed in their object class.
         - the metadata.MIxS object was invalid.")
  }
  
  for(nc in 1:ncol(df1)){ #merging columns with factors is always trouble, so convert to character
    if(class(df1[,nc])=="factor"){
      df1[,nc] <- as.character(df1[,nc])
    }
  }
  for(nc in 1:ncol(df2)){ #idem
    if(class(df2[,nc])=="factor"){
      df2[,nc] <- as.character(df2[,nc])
    }
  }
  shared_cols <- intersect(colnames(df1), colnames(df2))
  df11 <- df1[,!colnames(df1) %in% shared_cols, drop=FALSE]
  df22 <- df2[,!colnames(df2) %in% shared_cols, drop=FALSE]
  
  shared_rows <- intersect(rownames(df1), rownames(df2))
  
  if(length(shared_rows)>0){ #new column with original rownames
    df1$original_rowName <- row.names(df1)
    df2$original_rowName <- row.names(df2)
    shared_cols<-c(shared_cols, "original_rowName")
  }
  
  if(!ncol(df22)==0){ #if ncol(df22)==0, then df2 has no other colnames than df1
    df_UpRight <- data.frame(matrix(nrow = nrow(df1), ncol = ncol(df22), data=fill))
    colnames(df_UpRight) <- colnames(df22)
    rownames(df_UpRight) <- rownames(df1)
    df1_b <- cbind(df1, df_UpRight)
  }else{
    df1_b <- df1
  }
  
  if(!ncol(df11)==0){
    df11 <- df1[,!colnames(df1) %in% shared_cols, drop=FALSE]
    df_DownLeft <- data.frame(matrix(nrow = nrow(df2), ncol = ncol(df11), data=fill))
    rownames(df_DownLeft) <- rownames(df2)
    colnames(df_DownLeft) <- colnames(df11)
    df2_b <- cbind(df_DownLeft, df2)
  }else{
    df2_b <- df2
  }
  
  df_out <- data.frame(rbind(df1_b, df2_b))
  colnames(df_out)<-colnames(df1_b)
  df_out[c((nrow(df1_b)+1):(nrow(df_out))),shared_cols]<-df2[,shared_cols, drop=FALSE]
  
  if(!variables.as.cols){
    df_out <- data.frame(t(df_out))
  }
  
  if(class_out=="metadata.MIxS"){
    df_units <- c()
    df_section <- c()
    for(cname in colnames(df_out)){
      df_u1 <- d1@units[cname]
      df_u2 <- d2@units[cname]
      if(!is.na(df_u1) & !is.na(df_u2)){
        if(df_u1 != df_u2){
          stop("some variable units do not match. Correct and try again")
        }
      }
      if(is.na(df_u1)){
        df_u1 <- df_u2
      }
      df_units[cname] <- as.character(df_u1)
      
      df_s1 <- d1@section[cname]
      if(is.na(df_s1)){
        df_s1 <- d2@section[cname]
      }
      df_section[cname] <- as.character(df_s1)
    }
    
    df_t1<-d1@type
    df_t2<-d2@type
    if(df_t1=="flexible.MIxS" | df_t2=="flexible.MIxS"){
      df_t<-"flexible.MIxS"
    }else{
      df_t<-"strict.MIxS"
    }
    
    df_env1<-d1@env_package
    df_env2<-d2@env_package
    if(df_env1==df_env2){
      df_env<-df_env1
    }else{
      df_env<-"multiple_packages"
    }
    
    df_out <- new("metadata.MIxS",
                  data = df_out,
                  units = df_units,
                  section = df_section,
                  type = df_t,
                  env_package = df_env
    )
    
  }
  
  return(df_out)
  }

wideTab.to.hierarchicalTab <- function(dataTab, col_hierarchy){
  #' @author Maxime Sweetlove CC-BY 4.0 2019
  #' @description wideTab.to.hierarchicalTab turns a regular wide table into a hierarchical recursive table.
  #' @usage wideTab.to.hierarchicalTab(dataTab, col_hierarchy)
  #' 
  #' @param dataTab a data.frame. The wide table to be transformed, with columns as hierarchical variables.
  #' @param col_hierarchy a data.frame with two columns named "child" and "parent". The hierarchical relations between the columns formatted in a child-parent table listing the column names. Use "root" of NA for columns with nu perent.
  #' @details The input data will be transformed into a data.frame with three columns: one for the value (cell content in the origical data.frame), it's rank (the original column name), and parent (cell content of the column that is one up in the hierarchy)
  #' @example 
  #' Table <- data.frame(children=c("Agust", "Benny"), their_parents=c("Bernadette", "Rosa"), and_their_grandparents=c("Odille", "Jean-Pierre")) 
  #' TabHierarch <- data.frame(child=c("children", "their_parents", "and_their_grandparents"), parent=c("their_parents", "and_their_grandparents", "root"))
  #' wideTab.to.hierarchicalTab(Table, TabHierarch)
  
  if(colnames(col_hierarchy) != c("child", "parent")){
    stop("invalid input of the col_hierarchy argument.\nMust be a data.frame with two columns, names \"child\" and \"parent\".")
  }
  if(nrow(col_hierarchy) != ncol(dataTab)){
    stop("invalid input of the col_hierarchy argument.\n Every column in the dataTab input must be represented as a row in the col_hierarchy data.frame.")
  }
  
  
  recs_out <- data.frame(matrix(nrow=0, ncol=3))
  colnames(recs_out) <- c("value", "rank", "parent")
  
  dataTab <- dataTab[!duplicated(dataTab),]
  
  for(i in 1:nrow(dataTab)){
    for(j in 1:ncol(dataTab)){
      val <- as.character(dataTab[i,j])
      rank <- colnames(dataTab[,j, drop=FALSE])
      rankParent <- as.character(col_hierarchy[col_hierarchy$child == rank,]$parent)
      if(rankParent %in% colnames(dataTab)){
        valParent <- dataTab[i, rankParent]
      } else{
        valParent <- NA
      }
      recs_out <- rbind(recs_out,data.frame(value=val, rank=rank, parent=valParent))
    }
  }
  return(recs_out)
}

#==============================================================
### data content formating utils

coordinate.to.decimal<-function(val){
  #' @author Maxime Sweetlove CC-BY 4.0 2019
  #' @description degree.to.decimal turns a latutude or longitude value in a degrees-minutes-seconds (DMS) 
  #' format into a decimal value
  #' @usage degree.to.decimal(val)
  #' 
  #' @param val a character string. A single latitude or longitude value to be transformed.
  #' @details NSWE as well as degrees, minutes and seconds are recognized and turned into a numeric decimal coordinate value
  
  degreeSym <- c('°', '\302\260', '\241', "<U+00B0>", "\u00b0", "<c2><b0>", "<U+00C2>", "\u00c2", "\u00c2\u00b0")
  minSym <- c('\'', '\342\200\262', "'")
  secSym <- c('\"', '\342\200\263')
  
  val<-as.character(val)
  Encoding(val)<-"UTF-8"
  
  if(grepl("S|W", val)){
    s=-1
    val <- gsub("S|W", "", val)[[1]]
  }else{
    s=1
    val <- gsub("N|E", "", val)[[1]]
    }
  val <- gsub(paste(degreeSym, collapse="|"), "DDD", val)[[1]]
  val <- gsub(paste(minSym, collapse="|"), "MMM", val)[[1]]
  val <- gsub(paste(secSym, collapse="|"), "SSS", val)[[1]]
  val <- gsub(" ", "", val)[[1]]

  if(grepl("DDD", val)){
    degrees <- as.numeric(strsplit(val, "DDD")[[1]][1])
    degrees2 <- strsplit(val, "DDD")[[1]][2]
    if(grepl("MMM", degrees2)){
      minutes <- as.numeric(strsplit(degrees2, "MMM")[[1]][1])
      minutes2 <- strsplit(degrees2, "MMM")[[1]][2]
      if(grepl('SSS', minutes2)){
        seconds <- as.numeric(strsplit(minutes2, 'SSS')[[1]][1])
      }else{
        seconds <- 0
      }
    } else{
      minutes <- seconds <- 0
    }
    decimal <- s*(degrees + minutes/60 + seconds/3600)
  }else{
    decimal <- s*(as.numeric(val))
  }
  return(decimal)
}

parse.primer.text<-function(primerstring){
  #' @author Maxime Sweetlove CC-BY 4.0 2019
  #' @description parse.primer.text splits a text string into a forward and reverse primer sequence and name (depending on the content of the string)
  #' @usage parse.primer.text(primerstring)
  #' 
  #' @param primerstring a character string. A single sting to be parsed
  #' @details 
  # assumptions: DNA has no numbers and is of length >5 
  
  primerstring<-"AATGTACCTAGTGGTA GGTAGTAAYGTAG"
  
  DNAchars <- c("atcg", "ryswkmbdhvn")
  
  # 1. find the right separator, remove non-usefull chars
  #first look for a tab
  primerstring <- gsub('\t', " ", primerstring)
  primerstring <- gsub(",", " ", primerstring)
  primerstring <- gsub(";", " ", primerstring)
  primerstring <- gsub(":", " ", primerstring)
  primerstring <- gsub("|", " ", primerstring)
  primerstring <- gsub("\\s+", " ", primerstring) #collapse multiple spaces
  primerstring <- gsub("^\\s+|\\s+$", "", primerstring) #remove leading and trailing spaces
  
  primerstring <- gsub("(", "", primerstring)
  primerstring <- gsub(")", "", primerstring)
  primerstring <- gsub("[", "", primerstring)
  primerstring <- gsub("]", "", primerstring)
  primerstring <- gsub("{", "", primerstring)
  primerstring <- gsub("}", "", primerstring)
  
  primerstring <- strsplit(primerstring," ")[[1]] #split on space
  
  # 2. look if non-dna is present (e.g. primer names)
  
  # 3. spliting forward and everse primer
  

}

commonTax.to.NCBI.TaxID<-function(taxon, fill.unknown="12908"){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description commonTax.to.NCBI.TaxID converts taxon names of common taxa (superkingdom and phylum level) to it's NCBI taxID using an internal library. For taxa not in the internal library, please see https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
  #' @param taxon character or character vector. The taxon names to be converted to NCBI tax IDs.
  #' @param fill.unknown character. The string to return when a taxon was not found in the list with common taxa. Default is NCBI:txid12908 (the ID for "unknown sequences"), other options include NA or ""
  #' @details This function makes use of a limited internal library, not all the taxIDs are present
  #' @return a character or character vector with the matching taxon IDs, with fill.unknown for those not found. 
  #' @example commonTax.to.NCBI.TaxID(c("Bacteria", "Alveolata"), fill.unknown="")
  
  #note: ID 12908 stands for "unclassified sequences"
  taxIDs<-c()
  failed_taxa<-c()
  for(tx in taxon){
    if(tolower(tx) %in% rownames(TaxIDLib)){
      if(tolower(tx) %in% c("eukaryotes", "eukarya", "eukaryote")){
        tx<-"eukaryota"
      }
      taxIDs <- c(taxIDs, as.character(TaxIDLib[tolower(tx),]$NCBItaxID))
    }else{
      failed_taxa <- tx
      taxIDs <- c(taxIDs, fill.unknown)
    }
  }
  if(length(failed_taxa)>0){
    messagex <- paste("The following taxa were not found under the most common taxa: ", 
                      paste(unique(failed_taxa), collapse=", "), "\n an NCBI taxon ID for these taxa can be found at https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi",
                      sep="")
    warning(messagex)
  }
  return(taxIDs)
}

get.ENAName <- function(variable){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description get.ENAName gets the ENA variant a MIxS term
  #' @param variable character a MIxS term.
  ENAName <- as.character(TermsLib[TermsLib$name==variable,]$name_variant_ENA)
  return(ENAName)
}

parse.citation <- function(citation){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description parse.citation splits a character string with a scientific citation into the following components: autors, year, title, journal, issue, volume, pages, dio
  #' @param citation character string. A bibliographic citation, scientific reference
  #' @return a list with the autors, year, title, journal, issue, volume, pages and dio
  
  ref <- "Wiles, T. J., & Guillemin, K. J. (2020). Zebrafish as a Model for Investigating Animal–Microbe Interactions. In The Zebrafish in Biomedical Research (pp. 627-635). Academic Press."

  test <- bibentry(ref, bibtype="Article")
  format(test, "text")
  print(test, style = "citation")
  
  ## MLA
  # Wiles, Travis J., and Karen J. Guillemin. "Zebrafish as a Model for Investigating Animal–Microbe Interactions." The Zebrafish in Biomedical Research. Academic Press, 2020. 627-635.
  
  ## APA
  # Wiles, T. J., & Guillemin, K. J. (2020). Zebrafish as a Model for Investigating Animal–Microbe Interactions. In The Zebrafish in Biomedical Research (pp. 627-635). Academic Press.
  
  ## Chicago
  # Wiles, Travis J., and Karen J. Guillemin. "Zebrafish as a Model for Investigating Animal–Microbe Interactions." In The Zebrafish in Biomedical Research, pp. 627-635. Academic Press, 2020.
  
  ## Harvard
  # Wiles, T.J. and Guillemin, K.J., 2020. Zebrafish as a Model for Investigating Animal–Microbe Interactions. In The Zebrafish in Biomedical Research (pp. 627-635). Academic Press.
  
  ## Vancouver
  # Wiles TJ, Guillemin KJ. Zebrafish as a Model for Investigating Animal–Microbe Interactions. InThe Zebrafish in Biomedical Research 2020 Jan 1 (pp. 627-635). Academic Press.
  

  grepl('([12][0-9]{3})', "bbb(1999)p")
  
  
  ### split into author list, title and rest:
  
  
  ref <- strsplit(citation, ",")[[1]]
  
  # in author list
  gsub("and ", "", ref)
  gsub("& ", "", ref)
  
  ref <- unname(unlist(sapply(ref, function(x){strsplit(x, ".", fixed=TRUE)})))
  
  #remove leading and trailing spaces
  ref <- unname(unlist(sapply(ref, function(x){trimws(x, "both")})))
  
  ### guess author format:
  autors <- c()

  for(i in 1:length(ref)){
    
    if(ref[i] )
    
    print(ref[i])
    
    #replace &, and
    
  }
  
  
  if(nchar(ref[1])==1){
    # first initial(s), than name
  } else{
    # first name, (then initial(s))
    authors <- c(authors, ref[1]=ref[2])
  }
}






