# Sequence_Extraction2.R
#=========================================================

# path <- "/home/troy/Dropbox/alignmentfreer/data/abalone.gbk"

#=========================================================

# Obtain the phylogenetic information from a virus' gbk file.
.GetVirusPhylo <- function(templines, temp1, temp2) {
  templine <- paste(templines[temp1:temp2], collapse = "; ")
  templine <- substring(templine, 13, nchar(templine))
  templine <- strsplit(templine,";")[[1]]
  
  templine <- .RemoveSpacePeriod(templine)
  
  order <- .GetPhylo(templine[4:length(templine)], level = "order")
  family <- .GetPhylo(templine[4:length(templine)], level = "family")
  subfamily <- .GetPhylo(templine[4:length(templine)], level = "subfamily")
  genus <- .GetPhylo(templine2 = templine[4:length(templine)], level = "genus")
  balt2 <- templine[3]
  
  list(order = order, family = family, subfamily = subfamily,
       genus = genus, balt2 = balt2)
}

#=========================================================

# Obtain the phylogenetic information from a virus' gbk file.
.GetBaltimore <- function(balt1, balt2 = viralphylo$balt2,
                          family = viralphylo$family) {
  Balt=""
  
  if (balt2 == "Deltavirus") Balt <- "V"
  if (balt2=="dsDNA viruses, no RNA stage") Balt <- "I"
  if (balt2=="dsRNA viruses") Balt <- "III"
  if (balt2=="Satellites") Balt <- NA
  if (balt2=="ssDNA viruses") Balt <- "II"
  if (balt2=="ssRNA negative-strand viruses") Balt <- "V"
  if (balt2=="ssRNA positive-strand viruses, no DNA stage") Balt <- "IV"
  
  #Heterocapsa_circularisquama_RNA_virus_uid16157
  # J Gen Virol. 2011 Aug;92(Pt 8):1960-70. Epub 2011 May 11.
  # Three-dimensional reconstruction of Heterocapsa circularisquama RNA virus by electron cryo-microscopy.
  # Miller JL, Woodward J, Chen S, Jaffer M, Weber B, Nagasaki K, Tomaru Y, Wepf R, Roseman A, Varsani A, Sewell T.
  # http://www.ncbi.nlm.nih.gov/pubmed/21562120
  if (balt2=="unassigned ssRNA viruses") Balt <- "IV"
  if (balt2=="unassigned viruses") Balt <- ""
  if (balt2=="unclassified archaeal viruses") Balt <- ""
  if (balt2=="unclassified phages") Balt <- ""
  if (balt2=="unclassified virophages") Balt <- ""
  if (balt2=="unclassified viruses") Balt <- ""
  
  if (balt2=="Retro-transcribing viruses" & balt1=="DNA") Balt <- "VII"
  if (balt2=="Retro-transcribing viruses" & balt1=="ms-DNA") Balt <- "VII"
  if (balt2=="Retro-transcribing viruses" & balt1=="RNA") Balt <- "VI"
  if (balt2=="Retro-transcribing viruses" & balt1=="ss-RNA") Balt <- "VI"
  
  # 1078 Mouse_mammary_tumor_virus_uid14435 NC_001503.
  # http://en.wikipedia.org/wiki/Mouse_mammary_tumor_virus
  if (balt2=="Retro-transcribing viruses" & balt1=="ss-DNA") Balt <- "VI"
  
  #Redoing Retroviruses
  if (family == "Retroviridae") Balt <- "VI"
  if (family == "Metaviridae") Balt <- "VI"
  if (family == "Pseudoviridae") Balt <- "VI"  
  if (family == "Hepadnaviridae") Balt <- "VII"
  if (family == "Caulimoviridae") Balt <- "VII"

  Balt
}
