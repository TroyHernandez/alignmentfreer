# Sequence_Extraction3.R
#=========================================================

# path <- "/home/troy/Dropbox/alignmentfreer/data/abalone.gbk"

#=========================================================
# Remove spaces and periods from virusphlyo info
.RemoveSpacePeriod <- function(templine){
  for (i in 1:length(templine)) {
    temp <- templine[i]
    
    #---------------------------------------------------------
    # Remove period
    if (substring(temp, nchar(temp), nchar(temp)) == ".") {
      temp <- substring(temp, 1, nchar(temp) - 1)
    }
    
    #---------------------------------------------------------
    # Remove white spaces at beginning
    while (substring(temp, 1, 1) == " ") {
      temp <- substring(temp, 2, nchar(temp))
    }
    templine[i] <- temp
  }
  templine
}

#=========================================================
# Obtains phylogenetic info from templines
.GetPhylo <- function(templine2, level = c("order", "family",
                                          "subfamily", "order")) {
  # phylo scrape
  phylo <- ""
  .kPhyloSplit <- c(order = "virales", family = "viridae",
                   subfamily = "virinae", genus = "virus")
  split <- .kPhyloSplit[level]
  for(j in 1:length(templine2)) {
    if (templine2[j] != "" & phylo == "") {
      temp.split <- strsplit(templine2[j], split)[[1]]
      # Checking to see if temp was split
      if (nchar(templine2[j]) > nchar(temp.split)[1]) {
        phylo <- templine2[j]
      }
    }
  }
  # Removes unnecessary detail
  if (substr(phylo, 1, 2) == "un") {
    phylo <- strsplit(phylo, " ")[[1]][2]
  }
  
  # Fixes bug resulting from gbk file format
  if(level == "genus" && phylo == "ssDNA" | phylo == "dsDNA"){
    phylo <- ""
  }
  phylo
}
