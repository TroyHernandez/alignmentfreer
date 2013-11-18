# Vectorizer6.R

##############################################################

CalcNumLocIrreg <- function(irreg.name) {
  
  num.perms <- rep(1, nchar(irreg.name))
  extra.letters.vec <- rep(0, nchar(irreg.name))
  
  for (j in 1:nchar(irreg.name)) {
    # Iterate through kExtraLetters
    for (h in 1:3) {
      num.char <- nchar(irreg.name)
      extra.letter.pos <- which(substring(irreg.name,
                                          1:num.char,
                                          1:num.char)[j] == kExtraLetters[[h]])
      # Check to see letter is irregular;
      # create vector to determine nrow of matrix
      if (length(extra.letter.pos) > 0) {
        num.perms[j] <- h + 1
        extra.letters.vec[j] <- extra.letter.pos
      }
    }
  }
  list(num.perms = num.perms, extra.letters.vec = extra.letters.vec)
}
##############################################################
kTwoLetters <- c("R", "Y", "M", "K", "S", "W")
kThreeLetters <- c("B", "D", "H", "V")
kFourLetters <- c("N")
kExtraLetters <- list(kTwoLetters, kThreeLetters, kFourLetters)
# The x.letters.list correspond element-wise to x.letters.
kTwoLettersList <- list(c("G", "A"), c("T", "C"), c("A", "C"),
                        c("G", "T"), c("G", "C"), c("A", "T"))
kThreeLettersList <- list(c("A", "C", "T"), c("C", "T", "G"),
                          c("A", "C", "G"), c("A", "T", "G"))
kFourLettersList <- list(c("A", "C", "T", "G"))
kExtraLettersList <- list(kTwoLettersList, kThreeLettersList, kFourLettersList)

###############################################################

CalcPermutationMat <- function(irreg.name, num.perms) {
  num.char <- nchar(irreg.name)
  perm.mat <- matrix(substring(irreg.name, 1:num.char, 1:num.char),
                     nrow = prod(num.perms),
                     ncol = nchar(irreg.name), byrow = T)
  irreg.ind <- which(num.perms!=1)
  #We use nrow.perm.mat to carve up matrix into rows later
  nrow.perm.mat <- nrow(perm.mat)
  #Going through each extra letter
  for(j in 1:length(irreg.ind)){
    ir.ind.j <- irreg.ind[j]
    irreg.loc <-
      which(kExtraLetters[[num.perms[ir.ind.j] - 1]] == perm.mat[1, ir.ind.j])
    replacement.letters <-
      kExtraLettersList[[num.perms[ir.ind.j] - 1]][[irreg.loc]]
    #The letters that get reinserted are permuted here.
    real.replace <- c()
    for (l in 1:length(replacement.letters)) {
      real.replace <- c(real.replace,
                        rep(replacement.letters[l],
                            (nrow.perm.mat / length(replacement.letters))))
    }
    
    # Adjust the nrow.perm.mat to shrink with every replacement.
    nrow.perm.mat <- (nrow.perm.mat / length(replacement.letters))
    #insert real.replace into matrix
    perm.mat[, ir.ind.j] <- real.replace
  }
  perm.mat
}

###################################################################

AddToTblKmer <- function(perm.mat, kmers, reg.tbl, tbl,
                         kmer.seq, kmer.list, kmer.wt.list, irreg.name) {
  #paste letters back together and overwrite xac matrix
  overwrite <- apply(perm.mat, 1, paste, collapse = "")
  for (j in 1:length(overwrite)) {
    #Add to totals in reg.tbl
    overwrite.ind <- which(kmers == overwrite[j])
    reg.tbl[overwrite.ind] <-
      reg.tbl[overwrite.ind] + tbl[irreg.name] / length(overwrite)
    
    #Add positions to kmer.list 
    kmer.ind <- which(kmer.seq == irreg.name)
    kmer.list[[kmers[overwrite.ind]]] <- c(kmer.list[[kmers[overwrite.ind]]],
                                           kmer.ind)
    kmer.wt.list[[kmers[overwrite.ind]]] <-
      c(kmer.wt.list[[kmers[overwrite.ind]]],
        rep(1 / length(overwrite), length(kmer.ind)))
  }
  list(kmer.list = kmer.list, kmer.wt.list = kmer.wt.list, reg.tbl = reg.tbl)
}

weighted.var <- function(x, w, na.rm = FALSE) {
  # https://stat.ethz.ch/pipermail/r-help/2008-July/168762.html
  # By Gavin Simpson
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                                       na.rm)
}