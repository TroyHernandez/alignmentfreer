# Sequence_Extraction.R
#=========================================================

# path <- "/home/troy/Dropbox/alignmentfreer/data/abalone.gbk"

#=========================================================

#' Obtain the DNA sequence, and other information, from a gbk file.
#'
#' This function extracts a DNA sequence from a gbk file in addition to other 
#' information; e.g. gi number, accession number, number of base pairs,
#' and organism name.  A \code{gbk} class is output.  
#'
#' @param path the path to a gbk file
#' @param upper logical indicating if lower-case nucleotides should be converted
#' to upper-case
#' @param phylo string or \code{NA} indicating if phylogenetic information should be returned.
#' Currently, only \code{"virus"} option is working.
#' @keywords vectorizer
#' @export
#' 
gbk <- function(path, upper = FALSE, phylo = NA){
  
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
  
  #---------------------------------------------------------
  # Finding and writing dna sequence  
  templine <- character(0)
  
  for (i in temp0:(length(templines) - 1)) {
    dnaseq <- strsplit(paste(templines[i]), " ")[[1]]
    dnaseq <- paste(dnaseq[(max(which(dnaseq == "")) + 2):length(dnaseq)],
                    collapse = "")
    templine <- paste(templine, dnaseq, sep = "")
  }
  
  dnaseq <- templine
  
  #---------------------------------------------------------
  # Attributes assigned here
  
  # Number of basepairs
  temp <- strsplit(templines[1], " ")
  temp <- temp[[1]][temp[[1]] != ""]
  bp <- as.numeric(temp[3])
  
  #---------------------------------------------------------
  # Finding and writing gi/accession
  temp <- strsplit(templines[temp3], " ")
  temp <- temp[[1]][temp[[1]] != ""]
  accession <- temp[2]
  gi <- as.integer(substring(temp[3], 4, nchar(temp[3])))
  
  #---------------------------------------------------------
  # organism extraction
  organism <- substr(templines[temp1], 13, nchar(templines[temp1]))
  
  #---------------------------------------------------------
  # UnCapitalizing each letter
  if(upper == TRUE) {
    ltrs <- list(LETTERS)[[1]]
    names(ltrs) <- letters
    gbk <- paste(ltrs[strsplit(dnaseq, "")[[1]]], collapse = "")
  } else {
    gbk <- dnaseq
  }
  
  #---------------------------------------------------------
  # Extracting ORGANISM/phylogenetic info
  if(is.na(phylo)) {
    phylo <- list()
  } else {
    if (phylo == "virus") {
      # Get first part of baltimore predictor
      temp <- strsplit(templines[1], " ")
      temp <- temp[[1]][temp[[1]] != ""]
      balt1 <- temp[5]
      
      viralphylo <- .GetVirusPhylo(templines, temp1, temp2)
      
      baltimore <- .GetBaltimore(balt1, balt2 = viralphylo$balt2,
                                 family = viralphylo$family)
      phylo <- list(baltimore = baltimore, order = viralphylo$order,
                    family = viralphylo$family, subfamily = viralphylo$subfamily,
                    genus = viralphylo$genus)
    }
  }
  
  attributes(gbk) <- list(gi = gi, accession = accession,
                          bp = bp, organism = organism, phylo = phylo)
  class(gbk) <- "gbk"
  gbk
}

#=========================================================

print.gbk <- function(x, ...) {
  print(as.character(x))
}

#=========================================================

# path <- "/home/troy/Dropbox/alignmentfreer/data/abalone.fasta"

#=========================================================
#' Obtain the DNA sequence, and other information, from a fasta file.
#'
#' This function extracts a DNA sequence from a fasta file in addition to other 
#' information; e.g. gi number, accession number, number of base pairs,
#' and organism name.  A \code{fasta} class is output.  
#'
#' @param path the path to a gbk file
#' @param lower logical indicating if upper-case nucleotides should be converted
#' to lower-case
#' @keywords vectorizer
#' @export
#' 
fasta <- function(path, lower = FALSE){
  
  # Reading in lines
  templines <- readLines(path)
  
  info <- strsplit(templines[1], "|", fixed = TRUE)
  gi <- as.integer(info[[1]][2])
  accession <- info[[1]][4]
  # Removing white-space
  organism <- substring(info[[1]][5], 2, nchar(info[[1]][5]))
  
  #---------------------------------------------------------
  # Finding and writing dna sequence  
  templine <- character(0)
  
  for (i in 2:(length(templines) - 1)) {
    dnaseq <- strsplit(paste(templines[i]), " ")[[1]]
    templine <- paste(templine, dnaseq, sep = "")
  }
  
  dnaseq <- templine
  
  bp <- nchar(dnaseq)

  #---------------------------------------------------------
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
  class(fasta) <- "fasta"
  fasta
}

#=========================================================

print.fasta <- function(x) {
  print(as.character(x))
}
