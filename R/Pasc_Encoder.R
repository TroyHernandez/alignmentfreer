# Pasc_Encoder.R

#=========================================================

#' Convert a Pasc similarity file into a matrix. 
#'
#' This function converts a similarity file as obtained from the PASC tool and
#' converts it into a matrix and a corresponding data frame.
#'
#' @param path the path to a pasc similarity file
#' @param similarity logical determining if similarity matrix or dissimilarity
#' matrix should be returned
#' @keywords vectorizer
#' @references http://www.ncbi.nlm.nih.gov/sutils/pasc/viridty.cgi
#' @export
#' 
pasc <- function(path, similarity = TRUE) {
  # Reading in lines
  templines <- read.csv(path, sep = "\t", header = FALSE)
  name <- strsplit(path, "/")[[1]]
  name <- name[length(name)]
  name <- substring(name, 1, nchar(name) - 4)
  
  #gi1 is the first element in the column
  gi1 <- strsplit(as.character(templines[, 2]), "|", fixed = T)
  gi1 <- unlist(gi1)
  gi1 <- gi1[seq(2, length(gi1), 2)]
  gi1 <- as.numeric(gi1)
  
  #gi2 is the second element in the column
  gi2 <- strsplit(as.character(templines[, 9]), "|", fixed = T)
  gi2 <- unlist(gi2)
  gi2 <- gi2[seq(2, length(gi2), 2)]
  gi2 <- as.numeric(gi2)
  
  #orders according to first element and then according to the second element
  ord <- order(gi1, gi2)
  
  #Number of viruses
  nvirus <- length(unique(gi1)) + 1
  #Similarity matrix
  pasc <- dist(matrix(0, nrow = nvirus, ncol = 1))
  
  #Reordered templines matrix
  templines2=templines[ord, ]
  
  #dist objects work with this weird, but very helpful, syntax
  pasc[1:length(pasc)] <- templines2[ , 1]
  
  #Fills in the lower triangular section of the similarity matrix
  pasc <- as.matrix(as.dist(pasc, upper = T, diag = T), nrow = nvirus)
  #Fills in the diagonal of the matrix
  diag(pasc) <- 1
  
  if (similarity == F) {
    pasc = 1 - pasc
  }
  rownames(pasc) <- colnames(pasc) <- c(gi1[ord][1], gi2[ord][1:(nvirus - 1)])
  
  #indx is the classification information of the first virus.
  indx <- templines2[1, 3:7]
  
  #This must be done to add rest of viruses to indx
  colnames(indx) <- colnames(templines2)[10:14] <- c("Isolate", "Species",
                                                     "Genus", "Subfamily",
                                                     "Family")
  #The rest of the viruses are added in here from the pairing with virus 1.
  indx <- rbind(indx, templines2[1:(nvirus - 1), 10:14])
  
  indx <- cbind(rownames(pasc), indx, name)
  
  colnames(indx)[c(1, 7)] <- c("GI", "Family/Genus")
  rownames(indx) <- 1:nvirus
  
  attributes(pasc) <- list(index = indx, similarity = similarity)
  class(pasc) <- "pasc"
  pasc
}

#=========================================================

print.pasc <- function(x, ...) {
  row.names <- attributes(x)$index[,1]
  x <- matrix(x, nrow = sqrt(length(x)))
  rownames(x) <- colnames(x) <- row.names
  print(x)
}
