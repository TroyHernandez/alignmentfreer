# Sequence_Extraction.R

path <- "/home/troy/Dropbox/alignmentfreer/data/abalone.gbk"

GetSeqGbk <- function(path, upper = FALSE){
  
  # Reading in lines
  templines <- readLines(path)
  
  # Document reference features are extracted here 
  temp2 <- 0
  
  for (i in 1:(length(templines) - 1)) {
    if (strsplit(templines[i], "  ")[[1]][1] == "ORIGIN") {
      temp0 <- i + 1
    }
    if (strsplit(templines[i], "  ")[[1]][2] == "ORGANISM") {
      temp1 <- i
    }
    if (temp2 == 0 & strsplit(templines[i], "  ")[[1]][1] == "REFERENCE") {
      temp2 <- i - 1
    }
    if (strsplit(templines[i], "  ")[[1]][1] == "VERSION") {
      temp3 <- i
    }
  }
  
  ########################################################
  # Finding and writing dna sequence  
  templine <- character(0)
  
  for (i in temp0:(length(templines) - 1)) {
    dnaseq <- strsplit(paste(templines[i]), " ")[[1]]
    dnaseq <- paste(dnaseq[(max(which(dnaseq == "")) + 2):length(dnaseq)],
                    collapse = "")
    templine <- paste(templine, dnaseq, sep = "")
  }
  
  dnaseq <- templine
  
  ########################################################
  # Attributes assigned here
  
  # Number of basepairs
  temp <- strsplit(templines[1], " ")
  temp <- temp[[1]][temp[[1]] != ""]
  bp <- as.numeric(temp[3])
  
  ########################################################
  # Finding and writing gi/accession
  temp <- strsplit(templines[temp3], " ")
  temp <- temp[[1]][temp[[1]] != ""]
  accession <- temp[2]
  gi <- as.integer(substring(temp[3], 4, nchar(temp[3])))
  
  ########################################################
  # organism extraction
  organism <- substr(templines[temp1], 13, nchar(templines[temp1]))
  
  ########################################################
  # UnCapitalizing each letter
  if(upper == TRUE) {
    ltrs=list(LETTERS)[[1]]
    names(ltrs)=letters
    gbk <- paste(ltrs[strsplit(dnaseq, "")[[1]]], collapse = "")
  } else {
    gbk <- dnaseq
  }
  
  attributes(gbk) <- list(gi = gi, accession = accession,
                          bp = bp, organism = organism)
  class(gbk) <- "GBK"
  gbk
}

print.GBK <- function(x, ...) {
  #   cat("ROC curve: ")
  print(as.character(x))
}

path <- "/home/troy/Dropbox/alignmentfreer/data/abalone.fasta"

GetSeqFasta <- function(path, lower = FALSE){
  
  # Reading in lines
  templines <- readLines(path)
  
  info <- strsplit(templines[1], "|", fixed = TRUE)
  gi <- as.integer(info[[1]][2])
  accession <- info[[1]][4]
  # Removing white-space
  organism <- substring(info[[1]][5], 2, nchar(info[[1]][5]))
  
  ########################################################
  # Finding and writing dna sequence  
  templine <- character(0)
  
  for (i in 2:(length(templines) - 1)) {
    dnaseq <- strsplit(paste(templines[i]), " ")[[1]]
    templine <- paste(templine, dnaseq, sep = "")
  }
  
  dnaseq <- templine
  
  bp <- nchar(dnaseq)

  ########################################################
  # UnCapitalizing each letter
  
  if(lower == TRUE) {
    ltrs=list(letters)[[1]]
    names(ltrs)=LETTERS
    fasta <- paste(ltrs[strsplit(dnaseq, "")[[1]]], collapse = "")
  } else {
    fasta <- dnaseq
  }
  
  attributes(fasta) <- list(organism = organism, accession = accession,
                            gi = gi, bp = bp)
  class(fasta) <- "FASTA"
  fasta
}

print.FASTA <- function(x) {
#   cat(attributes(x)$organism)
  print(as.character(x))
}
