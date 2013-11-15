# Vectorizer4.R

#######################################################

CalcMeanKmerList <- function(kmer.seq, tbl, statistic, kmers) {
  kmer.list <- as.list(rep(0, length(kmers)))
  ans.temp <- rep(0, length(kmers))
  if (statistic > 1) {
    for (i in 1:length(kmers)) {
      temp <- which(kmer.seq == kmers[i]) - 1
      ans.temp[i] <- mean(temp) / length(kmer.seq)
      #Handle kmer.list in uniform(0, 1) for computational efficiency
      kmer.list[[i]] <- ((which(kmer.seq == kmers[i]) - 1) / length(kmer.seq))
      - ans.temp[i]
    }
  }
  
  list(mean = ans.temp, list = kmer.list)
}

#########################################################################

CalcDescriptiveStats <- function(ans, kmer.list, statistic, method) {
  if (statistic > 2) {
    for (j in 3:statistic) {
      if (method == "Moment") {
        ans[j, ] <- unlist(lapply(lapply(kmer.list, "^", (j - 1)), 
                                  sum)) / ans[1, ]
      } else if (method == "Sufficient") {
        ans[j, ] <- SufficientStats(kmer.list, j)
      } else if (method == "Lmoment") {
        ans[j, ] <- LmomentStats(kmer.list, j)
      }
    }
  }
  
  ans
}

##########################################################################
CorrectZeroCases <- function(ans, statistic, method = "Moment") {
  # Changes mean and moment entries with zero counts to
  # halfway mean and 0 moment
  if (sum(ans[1, ] == 0) > 0) {
    ind <- which(ans[1, ] == 0)
    ans[2, ind] <- .5
    if (statistic > 2 & method == "Moment") {
      ans[3:statistic, ind] <- 0
    }
    if (statistic > 2 & method == "Sufficient") {
      ans[3, ind] <- 0
      if (statistic > 3) {
        ans[4, ind] <- 0
      }
      if (statistic > 4) {
        ans[5, ind] <- -6/5
      }
    }
  }
  
  ans
}

#########################################################################

CorrectSingletonCases <- function(ans, statistic, method = "Moment") {
  #Changes mean and moment entries with single counts to halfway mean and 0 moment
  if (sum(ans[1, ] > 0 & ans[1, ] <= 1) > 0) {
    ind <- which(ans[1, ] > 0 & ans[1, ] <= 1)
    if (statistic > 2 & method == "Moment") {
      ans[3:statistic, ind] <- 0
    }
    if (statistic > 2 & method == "Sufficient") {
      ans[3, ind] <- 0
      if (statistic > 3) {
        ans[4, ind] <- 0
      }
      if (statistic > 4) {
        ans[5, ind] <- -6 / 5
      }
    }
  }
  
  ans
}

#######################################################################

CalcIrregMeanKmerList <- function(kmer.seq, tbl, irreg.names, kmers, ans) {
  
  # Initialize with zeroes as place holders
  kmer.list <- as.list(rep(0, length(kmers)))
  names(kmer.list) <- kmers
  reg.tbl <- rep(0, length(kmers))
  names(reg.tbl) <- kmers
  
  #Check each name not in list
  for (i in 1:length(irreg.names)) {
    irreg.name <- irreg.names[i]
    ind <- which(kmers == irreg.name)
    #if the name DNE
    if (length(ind) < 1) {
      allocate.irreg.kmers <- AllocateIrregKmers(irreg.name, kmers, reg.tbl,
                                                 tbl, kmer.seq, kmer.list)
      reg.tbl <- allocate.irreg.kmers$reg.tbl
      kmer.list <- allocate.irreg.kmers$kmer.list
    }else{
      reg.tbl[ind] <- reg.tbl[ind] + tbl[irreg.name]
      kmer.list[[ind]] <- c(kmer.list[[ind]], which(kmer.seq == kmers[ind]))
    }
  }
  
  #Remove initial zeroes
  for(i in 1:length(kmers)){
    kmer.list[[i]] <- kmer.list[[i]][-1]
  }
  
  tbl <- reg.tbl
  ans[1, ] <- tbl
  ans[2, ] <- unlist(lapply(kmer.list,sum))/reg.tbl

  list(kmer.list = kmer.list, tbl = reg.tbl, ans = ans)
}
