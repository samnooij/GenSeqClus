# Simple script that opens a csv file, uses the values in it
# to calculate euclidean distances between samples and
# stores those distances (as an R "dist object") as RDS file,
# and as a dataframe in a .csv file.

# Sam Nooij, 13 April 2018

library(reshape2)

frequencies_df_list <- snakemake@input[[1]]

for (frequencies_df in frequencies_df_list) {
  frequencies_df <- read.csv(frequencies_df, header = TRUE, sep = ",")
  distances <- as.matrix(dist(x = frequencies_df, method = "euclidean"))
  rownames(distances) <- frequencies_df$name
  colnames(distances) <- frequencies_df$name
  
  # Extract the index by checking which value (index) in the list
  # matches the current dataframe:
  index <- which(frequencies_df_list == frequencies_df)
  
  output_distances_RDS <- snakemake@output[["matrix"]][index]
  saveRDS(object = distances, file = output_distances_RDS)
  
  output_distances_csv <- snakemake@output[["dataframe"]][index]
  distances_df <- melt(distances)
  write.table(x = distances_df,
    file = output_distances_csv, 
    sep=",", 
    row.names = FALSE)
}