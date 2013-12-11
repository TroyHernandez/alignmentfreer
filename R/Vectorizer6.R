# Vectorizer6.R
#=========================================================

.CalcNumLocAmbig <- function(ambig.name) {
  
  num.perms <- rep(1, nchar(ambig.name))
  extra.letters.vec <- rep(0, nchar(ambig.name))
  
  for (j in 1:nchar(ambig.name)) {
    # Iterate through .kExtraLetters
    for (h in 1:3) {
      num.char <- nchar(ambig.name)
      extra.letter.pos <- which(substring(ambig.name,
                                          1:num.char,
                                          1:num.char)[j] == .kExtraLetters[[h]])
      # Check to see letter is ambigular;
      # create vector to determine nrow of matrix
      if (length(extra.letter.pos) > 0) {
        num.perms[j] <- h + 1
        extra.letters.vec[j] <- extra.letter.pos
      }
    }
  }
  list(num.perms = num.perms, extra.letters.vec = extra.letters.vec)
}
#=========================================================

.kTwoLetters <- c("R", "Y", "M", "K", "S", "W")
.kThreeLetters <- c("B", "D", "H", "V")
.kFourLetters <- c("N")
.kExtraLetters <- list(.kTwoLetters, .kThreeLetters, .kFourLetters)
# The x.letters.list correspond element-wise to x.letters.
.kTwoLettersList <- list(c("G", "A"), c("T", "C"), c("A", "C"),
                        c("G", "T"), c("G", "C"), c("A", "T"))
.kThreeLettersList <- list(c("A", "C", "T"), c("C", "T", "G"),
                          c("A", "C", "G"), c("A", "T", "G"))
.kFourLettersList <- list(c("A", "C", "T", "G"))
.kExtraLettersList <- list(.kTwoLettersList, .kThreeLettersList, .kFourLettersList)

#=========================================================

.CalcPermutationMat <- function(ambig.name, num.perms) {
  num.char <- nchar(ambig.name)
  perm.mat <- matrix(substring(ambig.name, 1:num.char, 1:num.char),
                     nrow = prod(num.perms),
                     ncol = nchar(ambig.name), byrow = T)
  ambig.ind <- which(num.perms!=1)
  #We use nrow.perm.mat to carve up matrix into rows later
  nrow.perm.mat <- nrow(perm.mat)
  #Going through each extra letter
  for(j in 1:length(ambig.ind)){
    ir.ind.j <- ambig.ind[j]
    ambig.loc <-
      which(.kExtraLetters[[num.perms[ir.ind.j] - 1]] == perm.mat[1, ir.ind.j])
    replacement.letters <-
      .kExtraLettersList[[num.perms[ir.ind.j] - 1]][[ambig.loc]]
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

#=========================================================

.AddToTblKmer <- function(perm.mat, kmers, reg.tbl, tbl,
                         kmer.seq, kmer.list, kmer.wt.list, ambig.name) {
  #paste letters back together and overwrite xac matrix
  overwrite <- apply(perm.mat, 1, paste, collapse = "")
  for (j in 1:length(overwrite)) {
    #Add to totals in reg.tbl
    overwrite.ind <- which(kmers == overwrite[j])
    reg.tbl[overwrite.ind] <-
      reg.tbl[overwrite.ind] + tbl[ambig.name] / length(overwrite)
    
    #Add positions to kmer.list 
    kmer.ind <- which(kmer.seq == ambig.name)
    kmer.list[[kmers[overwrite.ind]]] <- c(kmer.list[[kmers[overwrite.ind]]],
                                           kmer.ind)
    kmer.wt.list[[kmers[overwrite.ind]]] <-
      c(kmer.wt.list[[kmers[overwrite.ind]]],
        rep(1 / length(overwrite), length(kmer.ind)))
  }
  list(kmer.list = kmer.list, kmer.wt.list = kmer.wt.list, reg.tbl = reg.tbl)
}

#=========================================================

#' Compute a weighted variance.
#'
#' @param x an object containing the values whose weighted variance is to be
#' computed. 
#' @param w a numerical vector of weights the same length as x giving the
#' weights to use for elements of x. 
#' @param na.rm a logical value indicating whether NA values in x should be
#' stripped before the computation proceeds. 
#' 
#' @keywords vectorizer
#' @export
#'
#' @examples
#' 
#' 
#' ## GPA from Siegel 1994
#' wt <- c(5,  5,  4,  1)/15
#' x <- c(3.7,3.3,3.5,2.8)
#' xv <- weighted.var(x, wt)
#' 

#=========================================================

weighted.var <- function(x, w, na.rm = FALSE) {
  # https://stat.ethz.ch/pipermail/r-help/2008-July/168762.html
  # By Gavin Simpson; accessed 11.18.13
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  v1 <- sum(w)
  v2 <- sum(w ^ 2)
  x.bar <- sum(x * w) / v1
  m2 <- sum(w * (x - x.bar) ^ 2) / v1
  (v1 ^ 2 / (v1 ^ 2 - v2)) * m2
}

#=========================================================

#' Compute a weighted skewness.
#'
#' @param x an object containing the values whose weighted skewness is to be
#' computed. 
#' @param w a numerical vector of weights the same length as x giving the
#' weights to use for elements of x. 
#' @param na.rm a logical value indicating whether NA values in x should be
#' stripped before the computation proceeds. 
#' 
#' @keywords vectorizer
#' @export
#'
#' @examples
#' 
#' 
#' ## GPA from Siegel 1994
#' wt <- c(5,  5,  4,  1)/15
#' x <- c(3.7,3.3,3.5,2.8)
#' xs <- weighted.skew(x, wt)
#'
weighted.skew <- function(x, w, na.rm = FALSE) {
  # http://arxiv.org/abs/1304.6564
  # Lorenzo Rimoldini; accessed 11.19.13
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  v1 <- sum(w)
  v2 <- sum(w ^ 2)
  v3 <- sum(w ^ 3)
  x.bar <- sum(x * w) / v1
  m3 <- sum(w * (x - x.bar) ^ 3) / v1
  
  (v1 ^ 3 / (v1 ^ 3 - 3 * v1 * v2 + 2 * v3)) * m3
}

#=========================================================

#' Compute a weighted kurtosis.
#'
#' @param x an object containing the values whose weighted kurtosis is to be
#' computed. 
#' @param w a numerical vector of weights the same length as x giving the
#' weights to use for elements of x. 
#' @param na.rm a logical value indicating whether NA values in x should be
#' stripped before the computation proceeds. 
#' 
#' @keywords vectorizer
#' @export
#'
#' @examples
#' 
#' 
#' ## GPA from Siegel 1994
#' wt <- c(5,  5,  4,  1)/15
#' x <- c(3.7,3.3,3.5,2.8)
#' xk <- weighted.kur(x, wt)
#' 
weighted.kur <- function(x, w, na.rm = FALSE) {
  # http://arxiv.org/abs/1304.6564
  # Lorenzo Rimoldini; accessed 11.19.13
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  v1 <- sum(w)
  v2 <- sum(w ^ 2)
  v3 <- sum(w ^ 3)
  v4 <- sum(w ^ 4)
  x.bar <- sum(x * w) / v1
  m2 <- sum(w * (x - x.bar) ^ 2) / v1
  m4 <- sum(w * (x - x.bar) ^ 4) / v1
  denom <- (v1 ^ 2 - v2) *
    (v1 ^ 4 - 6 * v1 ^ 2 * v2 + 8 * v1 * v3 + 3 * v2 ^ 2 - 6 * v4)
  num1 <- v1 ^ 2 * (v1 ^ 4 - 4 * v1 * v3 + 3 * v2 ^ 2)
  num2 <- 3 * v1 ^ 2 * (v1 ^ 4 - 2 * v1 ^ 2 * v2 + 4 * v1 * v3 - 3 * v2 ^ 2)
  
  (num1 / denom * m4 - num2 / denom * m2 ^ 2)
}
