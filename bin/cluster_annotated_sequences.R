# Cluster calculations for annotated sequences (i.e. those for which
# genogroup/-type metadata is present). This script should prepare
# all the data that is necessary for the visualisations: histograms +
# cluster analysis information, distance maps and heatmaps + geno-
# group or -type annotation.

# Sam Nooij, 23 April 2018

## To do:
#   - fix the input for multiple kmers (1mers, 2mers, etc.)
#   - fix cluster calculations - either per df, or merge them
#   - write cluster statistics to output files

# Output:
# Table with results of cluster analysis calculations =
# clustering accordance, clustering cost, clustering quality(?),
# dunn index, which sequences were put in which cluster

### REQUIRED LIBRARIES ----------------------------------------------
library(reshape2) #for converting distance matrices - dataframes
library(dplyr) #for merging multiple dataframes


### READ INPUT FILES ------------------------------------------------

blast_df <- read.csv(snakemake@input[["blast"]],
  header = FALSE, sep = "\t")
column.names <-  c("query", "target", "identity", 
  "length", "mismatches", "gap.openings", "q.start", 
  "q.end", "t.start", "t.end", "evalue", "bit.score")
colnames(blast_df) <- column.names
blast_df <- within(blast_df, blast.distance <- 100 - identity)

nvr_dist_mat <- readRDS(snakemake@input[["nvr_mat"]])
nvr_dist_df <- read.csv(snakemake@input[["nvr_df"]],
  header = TRUE, sep = ",")
nvr_df <- melt(nvr_dist_mat)
colnames(nvr_df) <- c("query", "target", "nvr.distance")

## Note that there are multiple files for kmers, one for each k!
kmer_dist_mat <- readRDS(snakemake@input[["kmer_mat"]])
kmer_dist_df <- read.csv(snakemake@input[["kmer_df"]], 
  header = TRUE, sep = ",")

aa_ml_mat <- read.csv(snakemake@input[["aa_ml"]], 
  header = FALSE, sep = "", skip = 1, row.names = 1)
aa_ml_df <- melt(as.matrix(aa_ml_mat))
colnames(aa_ml_df) <- c("query", "target", "aa.ml.distance")

nt_ml_mat <- read.csv(snakemake@input[["nt_ml"]], 
  header = FALSE, sep = "", skip = 1, row.names = 1)
nt_ml_df <- melt(as.matrix(nt_ml_mat))
colnames(nt_ml_df) <- c("query", "target", "nt.ml.distance")

aa_ml_clean_mat <- read.csv(snakemake@input[["aa_ml_clean"]], 
  header = FALSE, sep = "", skip = 1, row.names = 1)
aa_ml_clean_df <- melt(as.matrix(aa_ml_clean_mat))
colnames(aa_ml_clean_df) <- c("query", "target", "aa.ml.clean.distance")

metadata_df <- read.csv(snakemake@input[["metadata"]],
  header = TRUE, sep = "\t")

###TESTING/DEBUGGING: ===============================================
setwd("/data/sapo/experiments/compare_clustering_methods/")
blast_df <- read.csv("tmp/SaV_genomes_tblastn_aa.blast", header = FALSE, sep = "\t")
column.names <-  c("query", "target", "identity", 
  "length", "mismatches", "gap_openings", "q.start", 
  "q.end", "t.start", "t.end", "evalue", "bitscore")
colnames(blast_df) <- column.names
blast_df <- within(blast_df, blast.distance <- 100 - identity)

nvr_dist_mat <- readRDS("tmp/SaV_genomes_nt_nvr_distances.RDS")
nvr_df <- melt(nvr_dist_mat)
colnames(nvr_df) <- c("query", "target", "nvr_distance")

### make a for-loop for kmers?

aa_ml_mat <- read.csv("tmp/SaV_genomes-mafft_aa.fas.mldist",
  header = FALSE, sep = "", skip = 1, row.names = 1)
colnames(aa_ml_mat) <- rownames(aa_ml_mat)
aa_ml_df <- melt(as.matrix(aa_ml_mat))
colnames(aa_ml_df) <- c("query", "target", "aa_ml_distance")

nt_ml_mat <- read.csv("tmp/SaV_genomes-mafft-RevTrans_aa.fas.mldist",
  header = FALSE, sep = "", skip = 1, row.names = 1)
colnames(nt_ml_mat) <- rownames(nt_ml_mat)
nt_ml_df <- melt(as.matrix(nt_ml_mat))
colnames(nt_ml_df) <- c("query", "target", "nt_ml_distance")

aa_ml_clean_mat <- read.csv("tmp/SaV_genomes-mafft-Gblocks-mafft_aa.fas.mldist",
  header = FALSE, sep = "", skip = 1, row.names = 1)
colnames(aa_ml_clean_mat) <- rownames(aa_ml_clean_mat)
aa_ml_clean_df <- melt(as.matrix(aa_ml_clean_mat))
colnames(aa_ml_clean_df) <- c("query", "target", "aa_ml_clean_distance")

metadata_df <- read.csv("data/SaV_genomes_metadata.tsv",
  header = TRUE, sep = "\t")

sample <- "SaV_genomes"

### END OF DATA IMPORTS =============================================

distances.df <- list(blast_df, nvr_df, aa_ml_df, nt_ml_df, aa_ml_clean_df) %>%
  Reduce(function(df1,df2) full_join(df1, df2, by = c("query", "target")), .)
distances.df <- merge.data.frame(x = distances.df, y = metadata_df, by.x = "query", by.y = "accession_id")

#Distances are all columns that contain "distance"
distances <- Filter(function(x) grepl("distance", x), colnames(distances.df))

#For all distances:
for (d in distances){
  print(d)
  #the distance measure is the name before "distance"
  measure <- sub("*.distance", "", d)
  #perform the cluster analysis
  sl_clustering <- single.linkage.clustering(df = distances.df, dist = d)
  
  #and save the results
  write.table(x = sl_clustering, file = paste("tmp/", sample, "_", measure, "_clustering_results.csv", sep = ""),
    row.names = FALSE)
}
### CLUSTER ANALYSIS FUNCTIONS --------------------------------------


### Copy from Cluster_analysis.Rmd
#{r calculate.clusters}
library(fpc) #for cluster.stats()

convert.to.matrix <- function(df, dist) {
  dist.mat <- acast(data = df[, c("query", "target", dist)],
    formula = query ~ target, fun.aggregate = mean)
  
  return(dist.mat)
}

#Find 'plateaus' or optima in the 'clusters/cutoff-curve'
find.consecutives <- function(x, k = 5) {
  #Find the length of vector x
  n <- length(x)
  #Create an empty object for end points
  plateau.values <- NULL
  #Loop over x, from start+5 to k
  for (i in k:n) {
    #if 5 or more elements before i are equal to i, save the endpoint
    if (all(x[(i-k):i] == x[i])) {
      plateau.values <- c(plateau.values, x[(i-k):i])
    }
  }
  
  return(plateau.values)
}

single.linkage.clustering <- function(df, dist) {
  #input: dataframe, distance measure
  #output: plot with marked plateaus/optima, dataframe with cluster statistics
  dist.mat <- convert.to.matrix(df, dist)
  
  sl.clustering <- hclust(dist(dist.mat), method = "single")
  
  #Set the maximum distance to be used as cutoff
  max.dist <- ceiling(max(df[dist]))
  
  #Select step sizes accordingly
  if (max.dist > 1000) {
    step.size <- 100
  } else if (max.dist > 40) {
    step.size <- 1
  } else if (max.dist > 4) {
    step.size <- .1
  } else {
    step.size <- .05
  }
  
  number.of.steps <- max.dist / step.size
  
  #Check the number of clusters for each selected cutoff value
  step <- 1
  
  #Prepare a matrix to hold the results
  cluster.analysis.mat <- matrix(nrow = number.of.steps + 1, ncol = 10)
  
  for (cutoff in seq(from = 0, to = max.dist, by = step.size)) {
    tree.cut <- cutree(sl.clustering, h = cutoff) 
    
    clustering.accordance <- calculate.clustering.accordance(tree.cut)
    
    #Note that quality is *per cluster* and cost *sums* the cost over all clusters,
    # so this function should be calculated per cluster, and summed.
    # Besides, recording qualities makes little sense if you report only 1 number.
    clustering.quality.cost <- calculate.cluster.quality.cost(dist = dist, clusters = tree.cut, cutoff = cutoff)
    #Make sure to pass "distances" as a vector of numeric values!
    clustering.cost <- clustering.quality.cost[[3]] #a single number
    clustering.quality <- clustering.quality.cost[c("cluster", "quality")] #a data frame of clusternumber - quality
    
    clustering.statistics <- calculate.cluster.statistics(dist.mat = dist.mat, cl.numbers = tree.cut)
    avg.size <- clustering.statistics$avg.size
    median.size <- clustering.statistics$median.size
    average.between <- clustering.statistics$average.between
    average.within <- clustering.statistics$average.within
    dunn.index <- clustering.statistics$dunn.index
    dunn2.index <- clustering.statistics$dunn2.index
    number.of.clusters <- max(tree.cut)
    
    #Write the information to the result table
    cluster.analysis.mat[step, 1] <- cutoff
    cluster.analysis.mat[step, 2] <- number.of.clusters
    cluster.analysis.mat[step, 3] <- avg.size
    cluster.analysis.mat[step, 4] <- median.size
    cluster.analysis.mat[step, 5] <- average.between
    cluster.analysis.mat[step, 6] <- average.within
    cluster.analysis.mat[step, 7] <- dunn.index
    cluster.analysis.mat[step, 8] <- dunn2.index
    cluster.analysis.mat[step, 9] <- clustering.accordance
    cluster.analysis.mat[step, 10] <- clustering.cost
    
    step <- step + 1
  }
  
  #Write the results to a dataframe
  cluster.analysis.df <- data.frame(cluster.analysis.mat)
  colnames(cluster.analysis.df) <- c("cutoff", "number.of.clusters",
    "average.size", "median.size", "average.between.distance", "average.within.distance",
    "dunn.index", "dunn2.index","clustering.accordance", "clustering.cost")
  
  return(cluster.analysis.df)
}

calculate.cluster.statistics <- function(dist.mat, cl.numbers) {
  #input: distance matrix (as used in the clustering)
  #       cluster numbers (e.g. the output of cutree, or the 'cluster' column from a kmeans-output dataframe)
  #output: list (dataframe)  of statistics:
  #       avg cluster size, median size, avg between dist, avg within dist, dunn index (x2)
  
  #Work-around for 'unclusterable' numbers:
  # check if cl.numbers is as long as dist.mat
  # (i.e. all sequences form their own cluster)
  # If it is, use cluster.stats(silhouette = FALSE)
  # (this will yield the same results, except "average.within" will be NaN)
  # else: proceed with the normal function
  
  if (length(cl.numbers) == length(rownames(dist.mat))) {
    statistics <- cluster.stats(d = dist.mat, clustering = cl.numbers, silhouette = FALSE)
  } else {
    statistics <- cluster.stats(d = dist.mat, clustering = cl.numbers)
  }
  
  # Calculate the statistics:
  avg.size <- mean(statistics$cluster.size)
  median.size <- median(statistics$cluster.size)
  
  stats.df <- as.data.frame(c(avg.size, median.size, 
    statistics[c("average.between", "average.within", "dunn", "dunn2")]))
  
  colnames(stats.df) <- c("avg.size", "median.size", "average.between", "average.within", "dunn.index", "dunn2.index")
  
  return(stats.df)
}

calculate.clustering.accordance <- function(tree.cut) {
  df <- as.data.frame(tree.cut)
  
  df$query <- as.factor(rownames(df))
  df <- merge(df, distances.df[c("query", "genogroup")], by = "query")
  colnames(df)[colnames(df) == "tree.cut"] <- "cluster"
  
  no.clusters <- max(df$cluster)
  no.groups <- length(levels(df$genogroup))
  
  if (no.groups <= no.clusters) {
    agg <- aggregate(cluster ~ genogroup, df, paste)
    X <- sum(rapply(agg$cluster, check.identical))
    Y <- no.groups - X
    Z <- no.clusters - X
  } else {
    agg <- aggregate(genogroup ~ cluster, df, paste)
    X <- sum(rapply(agg$genogroup, check.identical))
    Y <- no.groups - X
    Z <- no.clusters - X
  }
  
  clustering.accordance <- X / (Y + X + Z)
  
  return(clustering.accordance)
}

check.identical <- function(list) {
  if (length(levels(as.factor(list))) == 1) {
    identical <- TRUE
  } else {
    identical <- FALSE
  }
  
  return(identical)
}

fetch.cluster.distances <- function(tree.cut, distance) {
  df <- as.data.frame(tree.cut)
  colnames(df) <- "cluster"
  df$accession <- rownames(df)
  aggregated.df <- aggregate(accession ~ cluster, df, paste)
  
  grab.distances <- function(accessions) {
    distances <- convert.to.matrix(distances.df, distance)[accessions, accessions]
    return(distances)
  }
  
  distance.list <- lapply(aggregated.df$accession, FUN = grab.distances)
  # aggregate the rownames (IDs) per cluster number (look at above function)
  # use the IDs to extract distances from the original dataframe (using the distance metric provided to the function)
  # return: a list of distance lists(?) or matrices?
  return(distance.list)
}

#Note: For cluster quality, these distances should be intragroup distances!
# (i.e. within distances) This score counts *per cluster*, so reporting
# a single CQ per cutoff is useless.
calculate.cluster.quality.cost <- function(dist, clusters, cutoff) {
  ##dist = column required from overall dataframe
  ##clusters = cluster numbers as assigned by hclust/cutree
  ##cutoff = cutoff value between clusters, used for counting violations
  number.of.clusters <- max(clusters)
  
  #Failsafe for when cutoff = 0 (which is meaningless)
  if (cutoff == 0){
    return(c(data.frame(cluster = 1:number.of.clusters, quality = number.of.clusters * 0), NaN))
  } else {
    #Initiate a dataframe for storing quality scores
    cluster.qualities <- data.frame(cluster = 1:number.of.clusters, quality = number.of.clusters * NaN)
    
    #Start the calculations per cluster
    for (cluster in 1:number.of.clusters) {
      accession.ids <- attr(which(clusters == cluster), "names")
      #List all accession IDs in a cluster by finding ('which') all
      # "names" attributes that match cluster number 'cluster'
      within.distances <- subset(distances.df, query %in% accession.ids & target %in% accession.ids, dist)[,1]
      #Within-cluster distances is the subset of distances 
      # for which "query" and "target" match the accession IDs.
      #The object is turned to number only by appending [,1]
      # (instead of a dataframe)
      number.of.comparisons <- length(within.distances)
      
      ### Algorithm explained: ###
      #Sort the distances from large to small
      #Check if the number is greater than the cutoff (= violation)
      #  if not: break the loop (don't waste time)
      #  if greater: found a violation
      sorted.distances <- sort(within.distances, decreasing = TRUE)
      violations <- 0
      cluster.cost <- 0
      for (n in sorted.distances) {
        if (n < cutoff) {
          break
        } else {
          violations <- violations + 1
          cluster.cost <- cluster.cost + ((n - cutoff) / cutoff)
        }
      }
      #Calculate per-cluster quality
      cluster.quality <- 1 - violations / number.of.comparisons
      #And write it to the dataframe of collected qualities
      cluster.qualities$quality[cluster] <- cluster.quality
    }
  }
  
  return(c(cluster.qualities, cluster.cost))
  #Output = dataframe of qualities (per cluster), summed up cost
}
