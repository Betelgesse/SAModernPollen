rm(list=ls())
setwd ("")# set directory

# call required packages
library(randomForest)

Best_biome = load(file="") # Load the Random Forest model (.RData) for biomes or bioregions 
WK<- read.csv("WK_R.csv", row.names=NULL) # Load WK fossil data 
WKpollen <- WK[,-(1:6)] # Exclude unwanted columns
Polsum = WK[,6, drop=FALSE]
Pollen_matrix <- data.matrix(WKpollen) # Convert pollen data into matrix
WKAbundances <- Pollen_matrix/Polsum[row(Pollen_matrix)] # Calculate abundances by multiplying each taxon by the pollen sum
WKAbundances <- data.frame(WKAbundances)
WK_RF = WKAbundances

# Predicting past biomes/bioregions from fossil sequence
PredictedBiomes = WK$pred_bio = predict(best_bio_RF, WK_RF)
WK$Biomes = as.numeric(WK$pred_bio)
WK$prob_bio = predict(best_bio_RF, WK_RF, type="prob")
write.csv(WK, "NAME.csv")
write.csv(PredictedBiomes, "NAME.csv")
