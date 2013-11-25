alignmentfreer
==============

This is an R package that calculates alignment-free vectors, in particular the generalized vector, for dna sequences.  Functions also exist to extract relevent information from fasta and gbk files.

An existing Bioconductor package, Biostrings, can extract kmer counts, but leaves you to calculate the frequencies or other statistics.  Biostrings also ignores ambiguous nucleotides such as R, B, or N.  Additionally, there are no existing functions to extract sequence information from gbk files.

The functions GetSeqGbk and GetSeqFasta extract sequence and other information from gbk and fasta files, respectively.  The Vectorizer function will convert the DNA or RNA character string into the generalized vector described here:
http://arxiv.org/abs/1309.0408

The proper configuration of the options can allow the Vectorizer funtion to extract only kmer frequencies or the natural vector.  The xxx function can convert this to the composition vector.

