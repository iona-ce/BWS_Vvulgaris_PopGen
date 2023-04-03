# Working with the 01_BWS_DATA_COMPLETE.csv dataset.
# This file contains all data for every single trap and associated microsatellite information.

# The aim of this script is to create a file with genotype data to use in the rest of the study.
# Each sample must have at least 50% complete genotype to be used later on in the study.

# 1. Initial tidying ----

# Load data
bws_data <- read.csv("Data/01_BWS_DATA_COMPLETE.csv")

# Aim is to remove samples with fewer than 7 successfully amplified loci

# Firstly: replace any NAs that are not associated with loci with something other than NA

# Regions (2018 samples do not have regions as unecessary for later analyses)

bws_data$Region[is.na(bws_data$Region)] <- 1
bws_data$cluster[is.na(bws_data$cluster)] <- 1

bws_data$count_na_row <- (rowSums(is.na(bws_data)))/2 

# 2. Sample sizes

# For each method

nrow(subset(bws_data, bws_data$Method == "Dneasy")) #81
nrow(subset(bws_data, bws_data$Method == "Chelex")) #302

# 2. Removing samples with fewer than 7 successfully amplified microsatellites 

# Create new object to work with
bws_data_50 <- bws_data

# Remove of genotypes with fewer 7 loci
bws_data_50 <- subset(bws_data_50, bws_data_50$count_na_row <= 7) 

# Sample sizes

# Methods 
nrow(subset(bws_data_50, bws_data_50$Method == "Dneasy")) #79
nrow(subset(bws_data_50, bws_data_50$Method == "Chelex")) #222

# Year
nrow(subset(bws_data_50, bws_data_50$Year == "2017")) #105
nrow(subset(bws_data_50, bws_data_50$Year == "2018")) #196

# Looking at PCR success rates 

# Across all samples
pcr_success <- 1 - colSums(is.na(bws_data[,7:34]))/nrow(bws_data) 

# Depending on method 
bws_dneasy <- subset(bws_data, bws_data$Method == "Dneasy")
bws_chelex <- subset(bws_data, bws_data$Method == "Chelex")

pcr_success_dneasy <- 1 - colSums(is.na(bws_dneasy[,7:34]))/nrow(bws_dneasy)
pcr_success_chelex <- 1 - colSums(is.na(bws_chelex[,7:34]))/nrow(bws_chelex)

# Depending on year 
bws17_datacheck <- subset(bws_data, bws_data$Year == 2017) 
bws18_datacheck <- subset(bws_data, bws_data$Year == 2018)

pcr_success_17_datacheck <- 1 - colSums(is.na(bws17_datacheck[,7:34]))/nrow(bws17_datacheck)
pcr_success_18_datacheck <- 1 - colSums(is.na(bws17_datacheck[,7:34]))/nrow(bws17_datacheck)

# This is repeated after removal of incomplete genotypes to check which loci are suitable to keep for future analyses
bws17 <- subset(bws_data_50, bws_data_50$Year == 2017) 
bws18 <- subset(bws_data_50, bws_data_50$Year == 2018)

pcr_success_17 <- 1 - colSums(is.na(bws17[,7:34]))/nrow(bws17)
pcr_success_18 <- 1 - colSums(is.na(bws18[,7:34]))/nrow(bws18)

# In 2017 data, D3.13 and VMA4 have very low success rates - this is because these loci were not amplified 
# with DNeasy samples.
# In 2018 data, LIST2004 has only 0.4% success rate. This must be removed for future analyses.
# All other loci retained.

# Reformat data to be used later on 
bws_data_50[is.na(bws_data_50)] <- 0

#Remove this column 
bws_data_50$count_na_row <- NULL

#Save data
write.csv(bws_data_50, "Data/02_BWS_DATA_50.csv", row.names = FALSE)

# If need to read this in
bws_data_50 <- read.csv("Data/02_BWS_DATA_50.csv")





