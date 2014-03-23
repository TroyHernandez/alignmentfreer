# Vectorizer5.R
#=========================================================

.SufficientStats <- function(kmer.list = kmer.list, j = j){
  if (j == 3) {
    ss <- unlist(lapply(kmer.list, var))
  }
  if (j == 4) {
    ss <- unlist(lapply(kmer.list, e1071::skewness))
  }
  if (j == 5) {
    ss <- unlist(lapply(kmer.list, e1071::kurtosis))
  }
  if (j > 5) {
    stop("WARNING!!! Method Sufficient not able to handle
         statistics greater than 5; i.e. kurtosis.\n")
  }
  ss
}

#=========================================================

# .LmomentStats=function(kmer.list=kmer.list,j=j){
#   if(j>5){
#     cout("WARNING!!! Method Lmoments not currently able to handle dimension greater than 5; i.e. kurtosis.\n")
#   }else{
#     ss=unlist(lapply(kmer.list,Lmoments))[seq((j-1),4*length(kmer.list),4)]
#   }
#   ss
# }

#=========================================================

.AllocateAmbigKmers <- function(ambig.name, kmers, reg.tbl,
                               tbl, kmer.seq, kmer.list, kmer.wt.list) {
  #Check each letter of the name
  num.loc.ambig <- .CalcNumLocAmbig(ambig.name)
  num.perms <- num.loc.ambig$num.perms
#   extra.letters.vec <- num.loc.ambig$extra.letters.vec
  
  #substring(ambig.names[i],1:k,1:k)
  perm.mat <- .CalcPermutationMat(ambig.name,num.perms)
  
  add.to.tbl.kmer <- .AddToTblKmer(perm.mat, kmers, reg.tbl, tbl,
                                  kmer.seq, kmer.list, kmer.wt.list, ambig.name)
  kmer.list <- add.to.tbl.kmer$kmer.list
  kmer.wt.list <- add.to.tbl.kmer$kmer.wt.list
  reg.tbl <- add.to.tbl.kmer$reg.tbl
  
  list(kmer.list = kmer.list, kmer.wt.list = kmer.wt.list, reg.tbl = reg.tbl)
}

#=========================================================

.AmbigSufficientStats <- function(kmer.list, kmer.wt.list, statistic){
  if (statistic == 3) {
    ss <- unlist(lapply(1:length(kmer.list),
                        function(i, x, w) weighted.var(x[[i]],w[[i]]),
                        x = kmer.list, w = kmer.wt.list))
  }
  if (statistic == 4) {
    ss <- unlist(lapply(1:length(kmer.list),
                        function(i, x, w) weighted.skew(x[[i]],w[[i]]),
                        x = kmer.list, w = kmer.wt.list))
  }
  if (statistic == 5) {
    ss <- unlist(lapply(1:length(kmer.list),
                        function(i, x, w) weighted.kur(x[[i]],w[[i]]),
                        x = kmer.list, w = kmer.wt.list))
  }
  if (statistic > 5) {
    stop("WARNING!!! Method Sufficient not able to handle
         statistics greater than 5.\n")
  }
  ss
}
