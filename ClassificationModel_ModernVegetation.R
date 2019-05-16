# When using modern pollen  dataset for Southern Africa please acknowledge and cite: Sobol, M. K., Scott, L., & Finkelstein, S. A. (2019). Reconstructing past biomes states using machine learning and modern pollen assemblages: A case study from Southern Africa. Quaternary Science Reviews, 212, 1-17.
# This code was used for classification of both biomes and bioregions of Southern Africa by loading either either SAModern_Biome.csv or SAModern_Bioregion.csv

rm(list=ls())
setwd("") # set directory 

# Call required package
library(randomForest)

SAmod<- read.csv("") # Load modern pollen dataset
Pollen_only <- SAmod[,-(1:7)] # Pollen data only; gets rid of unwanted columns (row, column)
Polsum = SAmod[,7, drop=FALSE]
Pollen_matrix <- data.matrix(Pollen_only) # Converting pollen data into matrix
Abundances <- Pollen_matrix/Polsum$POLSUM # Calculating abundances by multiplying each taxon by the pollen sum
Abundances <- data.frame(Abundances)

# Function for excluding/including rare taxa 
exclude_rare_taxa = FALSE # Set to TRUE if rare taxa to be excluded
if (exclude_rare_taxa) {
    presence <- Pollen_matrix > 0
    Taxa.sum <- colSums(presence) # Count how many sites each taxa appears in
    Taxa.sum # Showing the list of taxa
    taxa.sum.mat <- data.matrix(lapply(Taxa.sum, unique)) # Iterate over the columns of Taxa.sum to find all unique values
    taxa.sum.mat # Showing the number of taxa
    taxa.sum.mat >=3 # Converting into binary again
    Abundances <- Abundances[,which(taxa.sum.mat >= 3)] # Replacing original pollen abundances data with pollen taxa present in at least 3 samples
    Pollen_filtered <- Pollen_matrix[,which(taxa.sum.mat >= 3)] # New matrix with pollen taxa present in at least 3 samples
    # Finding taxa with 3% in at least 1 sample
    max_abs <- apply(Abundances, 2, max) # Finding taxa with maximum abundaces across sites
    max_abs <- data.matrix(lapply(max_abs, unique))
    Abundances <- Abundances[,which(max_abs >= 0.0300)] # Filtering out abundances taxa without 3% in at least 1 site
    Pollen_filtered <- Pollen_filtered[,which(max_abs >= 0.0300)] # Filtering out Pollen_filtered taxa without 3% in at least 1 site
    Poll.sqr <- sqrt(Pollen_filtered)
} 
bio_n_data = Abundances
bio_n_data$BIO_N = SAmod$BIO_N

# Find which biomes/bioregions to remove
table(bio_n_data$BIO_N) # Show the number of site in each biome/bioregion
threshold = 4 # Set threshold for X number of pollen sites per biome/bioregions
counts = table(bio_n_data$BIO_N) # Show which biomes are kept and number of sites in each class

# Function for removing biomes/bioregions in which the site count is less than the threshold "threshol=X"
for (b in levels(bio_n_data$BIO_N)) {
    if (counts[b] < threshold) {
        bio_n_data = bio_n_data[bio_n_data$BIO_N != b,]
    }
}
bio_n_data = droplevels(bio_n_data) # Drop biomes/bioregions with occurances of "threshol=X" 

# Make a matrix to hold the class errors of each run
n_runs = 100 # Choose a number of models to run
n_levels = nlevels(bio_n_data$BIO_N)
bio_errors = matrix(, nrow=n_runs, ncol=n_levels)
colnames(bio_errors) = levels(bio_n_data$BIO_N)
best_bio_accuracy = 0
all_bio_accuracies = c() # For storing accuracies of each run
best_bio_accuracies = c() # For storing the best accuracie at the end of each run

# The biome/bioregion class labels are the last column of the data frame
bio_idx = ncol(bio_n_data)
tuna = tuneRF(bio_n_data[,-bio_idx], bio_n_data[,bio_idx], ntreeTry=500, plot=FALSE)
best_mtry = tuna[which.min(tuna[,2])] # Picks best number of pollen taxa
sampsize_accuracies = c()

for (i in c(1:n_runs)) {
    bio_n_RF = randomForest(BIO_N ~ ., data=bio_n_data, importance=T)
    con_mat = bio_n_RF$confusion
    bio_errors[i,] = con_mat[,(n_levels + 1)]
    
    # Calculate the out-of-bag (OOB) error
    con_mat = con_mat[,1:n_levels] # Get rid of the last column (class.error)
    
    num_correct_predictions = sum(diag(con_mat))
    num_total_predictions = sum(con_mat)
    oob_bio_accuracy = num_correct_predictions / num_total_predictions
    all_bio_accuracies = append(all_bio_accuracies, oob_bio_accuracy)
    
    if (oob_bio_accuracy > best_bio_accuracy) {
        best_bio_RF = bio_n_RF
        best_bio_accuracy = oob_bio_accuracy
    }
    best_bio_accuracies = append(best_bio_accuracies, best_bio_accuracy)
    print(best_bio_accuracy)
}

summary(all_bio_accuracies)
sd(all_bio_accuracies)
summary(bio_errors) # Provides basic stats (min, max, quantiles, mean)
best_bio_RF # Print the best Random Forest model statistics
print(best_bio_RF)# Print table
save(best_bio_RF, file="NAME.RData") # Save the Random Forest object
write.csv(best_bio_RF$confusion, file="NAME.csv")

# Importance variables
best_bio_RF$importance # Show importance variables for pollen taxa
write.csv(best_bio_RF$importance, file="NAME.csv")
importance_bio = varImpPlot(best_bio_RF)
dev.copy2eps(file="NAME.csv") 
save(importance_bio, file="NAME.RData")