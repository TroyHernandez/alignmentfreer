# Vectorizer5.R

#########################################################################
SufficientStats <- function(kmer.list = kmer.list, j = j){
  if (j == 3) {
    ss <- unlist(lapply(kmer.list, var))
  }
  if (j == 4) {
    ss <- unlist(lapply(kmer.list, skewness))
  }
  if (j == 5) {
    ss <- unlist(lapply(kmer.list, kurtosis))
  }
  if (j > 5) {
    cout("WARNING!!! Method Sufficient not able to handle
         statistics greater than 5; i.e. kurtosis.\n")
    stop()
  }
  ss
}

#########################################################################
# LmomentStats=function(kmer.list=kmer.list,j=j){
#   if(j>5){
#     cout("WARNING!!! Method Lmoments not currently able to handle dimension greater than 5; i.e. kurtosis.\n")
#   }else{
#     ss=unlist(lapply(kmer.list,Lmoments))[seq((j-1),4*length(kmer.list),4)]
#   }
#   ss
# }

#########################################################################

AllocateIrregKmers <- function(irreg.name, kmers, reg.tbl,
                               tbl, kmer.seq, kmer.list) {
  #Check each letter of the name
  num.loc.irreg <- CalcNumLocIrreg(irreg.name)
  num.perms <- num.loc.irreg$num.perms
#   extra.letters.vec <- num.loc.irreg$extra.letters.vec
  
  #substring(irreg.names[i],1:k,1:k)
  perm.mat <- CalcPermutationMat(irreg.name,num.perms)
  
  add.to.tbl.kmer <- AddToTblKmer(perm.mat, kmers, reg.tbl, tbl,
                                  kmer.seq, kmer.list, irreg.name)
  kmer.list <- add.to.tbl.kmer$kmer.list

  reg.tbl <- add.to.tbl.kmer$reg.tbl
  
  list(kmer.list = kmer.list, reg.tbl = reg.tbl)
}
