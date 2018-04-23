# Calculate kmer frequencies of
# the sequences in a given fasta file.
# Uses k = a list a numbers passed on by snakemake.

# Saves the kmer frequencies in .csv files and
# distance matrices as RDS files.

# Author: Sam Nooij
# Date: 6 Apr 2018

library('seqinr')
library('stringr')

kmers <- snakemake@params[["kmers"]]
input_file <- snakemake@input[[1]]

sequences <- read.fasta(file = input_file)

# Calculate for all given k sizes...
for (kmer in kmers) {
  # ...the frequencies of each kmer.
  frequencies <- function(sequence) {
    return(count(seq = sequence, wordsize = kmer))
  }
  # And save them in a dataframe.
  kmer_df <- as.data.frame(do.call(rbind, lapply(X = sequences, FUN = frequencies)))
  # Add a column with all the sequence IDs.
  kmer_df$name <- rownames(kmer_df)
  # And reorder to put the 'name' column first (before all the kmers).
  kmer_df <- kmer_df[c(ncol(kmer_df), 1:ncol(kmer_df) - 1)]
  
  # Save the kmer-count dataframe to a file:
  index <- which(kmers == kmer)
  output_kmers <- snakemake@output[[1]][index]
  print(paste("Now writing to file:", output_file))
  write.table(x = kmer_df, file = output_file, sep = ",", row.names = FALSE, col.names = TRUE)
}