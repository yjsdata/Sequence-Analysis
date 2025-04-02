library(haven)
library(TraMineR)
library(mice)
library(dplyr)
library(fpc)
library(cluster)

# Run multiple imputation on sequence columns
imputed_data <- mice(data[, 16:21], m = 6, method = "pmm", seed = 123)

# Get the complete (imputed) data for the first imputation
data_imputed <- complete(imputed_data, 1)

# Integrate imputed columns back into the original data frame
data[, 16:21] <- data_imputed

(data.alph <- seqstatl(data[, 16:21]))

##Rename
library(dplyr)
data <- data %>%
  rename("yr17" = F_3_leisuremp_17)
data <- data %>%
  rename("yr15" = F_3_leisuremp_15)
data <- data %>%
  rename("yr13" = F_3_leisuremp_13)
data <- data %>%
  rename("yr11" = F_3_leisuremp_11)
data <- data %>%
  rename("yr09" = F_3_leisuremp_09)
data <- data %>%
  rename("yr07" = F_3_leisuremp_07)

data.lab <- c("No LA", "Household Only", "Social Only", 
              "Social and Household", "Mental/PA Only", 
              "Mental/PA and Household", 
              "Mental/PA and Social", "ALL LA")
data.shortlab <- c("1", "2", "3", "4", "5", "6", "7", "8")
data.seq <- seqdef(data[, 16:21], alphabet = data.alph,
                   labels = data.lab, states = data.shortlab)

seqstatl(data[, 16:21])

data.alph <- c("1", "2", "3", "4", "5", "6", "7", "8")
data.seq <- seqdef(data[, 16:21], states = data.shortlab,
                   labels = data.lab, xtstep = 1,
                   alphabet = data.alph)

seqIplot(data.seq, border = NA)

seqIplot(data.seq, xtstep = 1, sortv = "from.start")

seqfplot(data.seq, border = NA)

seqdplot(data.seq, border = NA)

seqmsplot(data.seq, with.legend = "right")

seqmtplot(data.seq, ylim = c(0, 1))

selected_data_young <- subset(data, agegroup == 0)
selected_data_old <- subset(data, agegroup == 1)

###Matching###
# Compute the optimal matching dissimilarity matrix
om_dist <- seqdist(data.seq, method = "OM", indel = 1, sm = "CONSTANT", 
                   with.missing = TRUE)

# Perform hierarchical clustering
hclust_om <- hclust(as.dist(om_dist), method = "ward.D2")
plot(hclust_om)

# MDS plot to visualize sequence dissimilarities
mds_om <- cmdscale(as.dist(om_dist), k = 2)
plot(mds_om, xlab = "Dimension 1", ylab = "Dimension 2", 
     main = "MDS Plot of OM Distances")

#######
# Calculate ASW for each number of clusters (3 to 8)
asw_results <- data.frame(Clusters = 3:8, ASW = NA)  # Create a data frame to store results

for (k in 3:8) {
  # Assign clusters for k clusters
  clusters <- cutree(hclust_om, k = k)
  
  # Calculate silhouette values
  sil <- silhouette(clusters, dist = as.dist(om_dist))
  
  # Compute silhouette scores for different cuts
  sil_scores <- sapply(3:8, function(k) {
    clusters_k <- cutree(hclust_om, k = k)
    silhouette_score <- cluster.stats(as.dist(om_dist), clusters_k)$avg.silwidth
    return(silhouette_score)
  })
  
  # Compute the average silhouette width (ASW)
  asw_results$ASW[asw_results$Clusters == k] <- mean(sil[, 3])
}

# Print ASW results
print(asw_results)

# Plot ASW vs. number of clusters
plot(asw_results$Clusters, asw_results$ASW, type = "b", pch = 19, col = "blue",
     xlab = "Number of Clusters", ylab = "Average Silhouette Width (ASW)",
     main = "ASW for Different Numbers of Clusters")

# Define a custom clustering function that returns the appropriate structure
custom_cutree <- function(hclust_obj, k) {
  cluster_assignments <- cutree(hclust_obj, k = k)
  # Return a list with the cluster assignments as part of the object
  return(list(cluster = cluster_assignments))
}

# Compute the Gap Statistic using the hierarchical clustering result
gap_stat <- clusGap(om_dist, FUNcluster = function(x, k) custom_cutree(hclust_om, k), K.max = 8, B = 50)

# Plot the Gap Statistic
plot(gap_stat, main = "Gap Statistic for Different Cluster Numbers")
#######
#######
# Calculate Hubert's Gamma for each number of clusters (3 to 8)
hg_results <- data.frame(Clusters = 3:8, Hubert_Gamma = NA)  # Data frame to store results

for (k in 3:8) {
  # Assign clusters for k clusters
  clusters <- cutree(hclust_om, k = k)
  
  # Compute the dissimilarity matrix for cluster assignments
  # Convert cluster assignments to a binary membership matrix
  cluster_diss <- outer(clusters, clusters, FUN = function(x, y) as.numeric(x != y))
  
  # Flatten both dissimilarity matrices for correlation calculation
  om_diss_flat <- as.vector(as.dist(om_dist))
  cluster_diss_flat <- as.vector(as.dist(cluster_diss))
  
  # Compute Hubert's Gamma (correlation between dissimilarities)
hg_results$Hubert_Gamma[hg_results$Clusters == k] <- cor(om_diss_flat, cluster_diss_flat, method = "pearson")
}

# Print Hubert's Gamma results
print(hg_results)

# Plot Hubert's Gamma vs. number of clusters
plot(hg_results$Clusters, hg_results$Hubert_Gamma, type = "b", pch = 19, col = "red",
     xlab = "Number of Clusters", ylab = "Hubert's Gamma (HG)",
     main = "Hubert's Gamma for Different Numbers of Clusters")
#######
#######
# Calculate PBC for each number of clusters (3 to 8)
pbc_results <- data.frame(Clusters = 3:8, PBC = NA)  # Data frame to store results

for (k in 3:8) {
  # Assign clusters for k clusters
  clusters <- cutree(hclust_om, k = k)
  
  # Create a binary membership matrix (same cluster = 1, different = 0)
  cluster_membership <- outer(clusters, clusters, FUN = function(x, y) as.numeric(x == y))
  
  # Flatten both matrices for correlation calculation
  om_diss_flat <- as.vector(as.dist(om_dist))
  membership_flat <- as.vector(as.dist(cluster_membership))
  
  # Calculate the Point Biserial Correlation (PBC)
  pbc_results$PBC[pbc_results$Clusters == k] <- cor(om_diss_flat, membership_flat, method = "pearson")
}

# Print PBC results
print(pbc_results)

# Plot PBC vs. number of clusters
plot(pbc_results$Clusters, pbc_results$PBC, type = "b", pch = 19, col = "green",
     xlab = "Number of Clusters", ylab = "Point Biserial Correlation (PBC)",
     main = "PBC for Different Numbers of Clusters")
#######
#######
# Create a data frame to store metrics for comparison
comparison_results <- data.frame(
  Clusters = 3:8,
  ASW = asw_results$ASW,
  HG = hg_results$Hubert_Gamma,
  PBC = pbc_results$PBC
)

# Normalize metrics to range [0, 1] (Optional, for better visualization)
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
comparison_results$ASW_norm <- normalize(comparison_results$ASW)
comparison_results$HG_norm <- normalize(comparison_results$HG)
comparison_results$PBC_norm <- normalize(comparison_results$PBC)

# Print the table of results
print(comparison_results)

# Plot the metrics for comparison
plot(comparison_results$Clusters, comparison_results$ASW_norm, type = "b", pch = 19, col = "blue",
     xlab = "Number of Clusters", ylab = "Normalized Metric Values",
     main = "Comparison of Clustering Metrics")
lines(comparison_results$Clusters, comparison_results$HG_norm, type = "b", pch = 19, col = "red")
lines(comparison_results$Clusters, comparison_results$PBC_norm, type = "b", pch = 19, col = "green")

# Add a legend to the plot
legend("right", legend = c("ASW", "Hubert's Gamma", "PBC"), 
       col = c("blue", "red", "green"), lty = 1, pch = 19)

#######
#######
# Cut the dendrogram into a chosen number of clusters (5 clusters)
clusters_5 <- cutree(hclust_om, k = 5)
print(clusters_5)
table(clusters_5)

seqdplot(data.seq, group = data$clusters_5, border = NA, main = "Cluster", 
         with.legend = "down", cex.legend = 1.5)

seqIplot(data.seq, group = data$clusters_5, xtstep = 1, border = NA)

table(data$clusters_5)

# Sequence distribution plot by cluster
seqdplot(data.seq, group = clusters_5, border = NA, main = "Cluster", 
         with.legend = FALSE)

seqdplot(data.seq, group = clusters_5, border = NA, main = "Cluster", 
         with.legend = "right", cex.legend = 1.5)

data$clusters_5 <- clusters_5

# Print cluster size table
table(data$clusters_5)


# Cut the dendrogram into a chosen number of clusters (6 clusters)
clusters_6 <- cutree(hclust_om, k = 6)
print(clusters_6)
table(clusters_6)

seqdplot(data.seq, group = data$clusters_6, border = NA, main = "Cluster", 
         with.legend = "down", cex.legend = 1.5)

seqIplot(data.seq, group = data$clusters_6, xtstep = 1, border = NA)

table(data$clusters_5)

# Sequence distribution plot by cluster
seqdplot(data.seq, group = clusters_6, border = NA, main = "Cluster", 
         with.legend = FALSE)

seqdplot(data.seq, group = clusters_6, border = NA, main = "Cluster", 
         with.legend = "right", cex.legend = 1.5)

data$clusters_6 <- clusters_6

# Print cluster size table
table(data$clusters_6)

#######
#######
# Choose the best number of clusters based on the highest silhouette score or gap statistic
optimal_k <- which.max(sil_scores)  # or use the gap statistic result

# Cut the dendrogram at the optimal number of clusters
optimal_clusters <- cutree(hclust_om, k = optimal_k)

# Add cluster assignments to the data
data$optimal_clusters <- optimal_clusters

# Visualize the sequence data by the optimal clusters
seqdplot(data.seq, group = data$optimal_clusters, border = NA, main = "Cluster", with.legend = FALSE, cex.legend = 1.5)
seqdplot(data.seq, group = data$optimal_clusters, border = NA, main = "Cluster", with.legend = "right", cex.legend = 1.5)

# Print cluster size table
table(data$optimal_clusters)