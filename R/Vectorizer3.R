# Vectorizer3.R


#######################################################
#k-mer generator in caps
KmerGenerator <- function(kmer) {
  
  dna <- c("A","C","G","T")
  
  #kmercolnames is final output
  kmer.list <- dna
  
  if (kmer > 1) {
    for (j in 2:kmer) {
      kmer.list <- c(t(outer(kmer.list, dna, FUN = "paste", sep = "")))
    }
  }
  
  kmer.list
}

#######################################################
#Appends statistic name to kmer 
StatAppender <- function(kmers, statistic) {
  for (i in 1:statistic) {
    if (i == 1) {
      col.names <- paste("n", kmers, sep = "")
    }
    if (i == 2) {
      col.names <- c(col.names,paste("mu", kmers, sep = ""))
    }
    if (i == 3) {
      col.names <- c(col.names,paste("var", kmers, sep = ""))
    }
    if (i == 4) {
      col.names <- c(col.names,paste("skew", kmers, sep = ""))
    }
    if (i == 5) {
      col.names <- c(col.names,paste("kur", kmers, sep = ""))
    }
  }
  col.names
}

##############################################################
CalcRegLetters <- function(kmer.seq, tbl, statistic, kmers, method = "Moment") {
  ans <- matrix(0, nrow = statistic, ncol = length(kmers))
  ans[1, ] <- tbl
  
  mean.kmer.list <- CalcMeanKmerList(kmer.seq, tbl, statistic, kmers)
  ans[2, ] <- mean.kmer.list$mean
  kmer.list <- mean.kmer.list$list
  
  ans <- CalcDescriptiveStats(ans, kmer.list, statistic, method)
  
  if (statistic > 1) {
    ans <- CorrectZeroCases(ans, statistic)
    ans <- CorrectSingletonCases(ans, statistic)
  }
  
  # Changes counts to frequencies
  ans[1, ] <- ans[1, ] / sum(ans[1, ])
  
  ans <- c(t(ans))
}



################################################################

CalcIrregLetters <- function (kmer.seq, statistic, kmers, method = "Moment") {
  ans <- matrix(0, nrow = statistic, ncol = length(kmers))
  tbl <- table(kmer.seq)
  irreg.names <- names(tbl)
  
  irreg.kmer.list <- CalcIrregMeanKmerList(kmer.seq, tbl,
                                           irreg.names, kmers, ans)
  tbl <- irreg.kmer.list$tbl
  kmer.list <- irreg.kmer.list$kmer.list
  ans <- irreg.kmer.list$ans

  ans <- CalcDescriptiveStats(ans, kmer.list, statistic, method)
  
  if (statistic > 1) {
    ans <- CorrectZeroCases(ans, statistic)
    ans <- CorrectSingletonCases(ans, statistic)
  }
  
  ans[1, ] <- ans[1, ] / sum(ans[1, ])
  
  ans <- c(t(ans))
  ans
}

##############################################################
CalcEmptyLetters=function(x,tbl,d,kmers,KMERS,method="Moment"){
  ans=matrix(0,nrow=d,ncol=length(kmers))
  
  if(d>1){
    if(sum(ans[1,]==0)>0){
      ind=which(ans[1,]==0)
      ans[2,ind]=.5
      if(d>2 & method=="Moment"){
        ans[3:d,ind]=0
      }
      if(d>2 & method=="Sufficient"){
        ans[3,ind]=0
        if(d>3){
          ans[4,ind]=0
        }
        if(d>4){
          ans[5,ind]=-6/5
        }
      }
    }
    
    #Changes mean and moment entries with single counts to halfway mean and 0 moment
    ind=1:ncol(ans)
    if(d>2 & method=="Moment"){
      ans[3:d,ind]=0
    }
    if(d>2 & method=="Sufficient"){
      ans[3,ind]=0
      if(d>3){
        ans[4,ind]=0
      }
      if(d>4){
        ans[5,ind]=-6/5
      }
    }
    
  }
  
  ans=c(t(ans))
  ans
}
