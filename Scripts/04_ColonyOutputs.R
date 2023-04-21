# Looking at the outputs of COLONY

# Aim is to: 
# 1. Determine whether there is any evidence for siblings within/between traps
# We can expect that there will be more siblings within traps
# This is also an indication of whether the markers are appropriate for sibship analyses 

# 2. Remove any siblings from future analyses as the presence of siblings could interfere 
# with the results 

# 3. If there is any evidence for siblings between traps, determine the distance (in km) 
# between them

# 1. Load data ----

# Full siblings 2017
FScluster_bws2017 <- read.table("Data/05_COLONY/BWS17_colony_results.FullSibDyad", sep = ",", header = TRUE)
nrow(FScluster_bws2017) #32

# Half siblings 2017
HScluster_bws2017 <- read.table("Data/05_COLONY/BWS17_colony_results.HalfSibDyad", sep = ",", header = TRUE)
nrow(HScluster_bws2017) #201

# Full siblings 2018
FScluster_bws2018 <- read.table("Data/05_COLONY/BWS18_colony_results.FullSibDyad", sep = ",", header = TRUE)
nrow(FScluster_bws2018) #63

# Half siblings 2018
HScluster_bws2018 <- read.table("Data/05_COLONY/BWS18_colony_results.HalfSibDyad", sep = ",", header = TRUE)
nrow(HScluster_bws2018) #330

# 2. Tidy data for later analyses 

# Add columns indicating full/half sibship
FScluster_bws2017$SibType <- "Full"
FScluster_bws2018$SibType <- "Full"
HScluster_bws2017$SibType <- "Half"
HScluster_bws2018$SibType <- "Half"

# Adding year 
FScluster_bws2017$Year <- "2017"
FScluster_bws2018$Year <- "2018"
HScluster_bws2017$Year <- "2017"
HScluster_bws2018$Year <- "2018"

# Remove unecessary columns (here, we can remove maternal/paternal column from the half siblings tables: 
# vespine wasps can only be half siblings by the mother, not the father)

HScluster_bws2017$Paternal.Maternal <- NULL
HScluster_bws2018$Paternal.Maternal <- NULL

# Combine these data frames 
BWS_Sibs <- rbind(FScluster_bws2017, FScluster_bws2018, HScluster_bws2017, HScluster_bws2018)

# Add geographical information 

# Latitude and Longitude
match1 <- match(BWS_Sibs$OffspringID1, bws_data$Sample)
BWS_Sibs$Lat1 <- bws_data[match1,]$Latitude
BWS_Sibs$Long1 <- bws_data[match1,]$Longitude
match2 <- match(BWS_Sibs$OffspringID2, bws_data$Sample)
BWS_Sibs$Lat2 <- bws_data[match2,]$Latitude
BWS_Sibs$Long2 <- bws_data[match2,]$Longitude

# Cluster information (relevant to 2018 samples only)

#Adding cluster to the BWS Sibs 
BWS_Sibs$Cluster1 <- bws_data[match1,]$cluster
BWS_Sibs$Cluster2 <- bws_data[match2,]$cluster

# Adding trap number 
BWS_Sibs$Trap1 <- bws_data[match1,]$Trap_ID
BWS_Sibs$Trap2 <- bws_data[match2,]$Trap_ID

# Add completedness of genotypes - this will serve to later on make a decision on which
# samples to remove 
BWS_Sibs$count_na_row1 <- bws_data[match1,]$count_na_row
BWS_Sibs$count_na_row2 <- bws_data[match2,]$count_na_row

# Remove sib pairs if inferred probability of sibship is less than 0.8
BWS_Sibs <- subset(BWS_Sibs, BWS_Sibs$Probability > 0.85)

## 3. Pairwise distances: calculating distance between siblings in different traps ----

# Libraries 
library(dplyr)
library(geosphere)

# Calculate distance in km 
BWS_Sibs <- BWS_Sibs %>% rowwise() %>% 
  mutate(Distance = distGeo(c(Long1, Lat1), c(Long2, Lat2)))

#Divde this by 1000 to transform to metres 
BWS_Sibs$Distance <- BWS_Sibs$Distance/1000

# Make this into dataframe
BWS_Sibs <- data.frame(BWS_Sibs)

# Get data
nrow(subset(BWS_Sibs, BWS_Sibs$SibType == "Full")) #24
nrow(subset(BWS_Sibs, BWS_Sibs$SibType == "Half")) #1
length(unique(c(BWS_Sibs$OffspringID1, BWS_Sibs$OffspringID2))) # 34

# There are 31 full sib pairs, 1 half sib pair comprising 39 individuals 

# Save file 
write.csv(BWS_Sibs, "Data/06_BWS_Sibling_Results.csv", row.names = FALSE)
BWS_Sibs <- read.csv("Data/06_BWS_Sibling_Results.csv")

# Edit this document for Supplementary Material 
BWS_Sibs_Supplementary <- BWS_Sibs[c("Year", "OffspringID1", "Trap1", "Cluster1", "OffspringID2", "Trap2", "Cluster2", "Probability", "SibType", "Distance")]

# Replace cluster numbers 

BWS_Sibs_Supplementary$Cluster1[BWS_Sibs_Supplementary$Cluster1 == 1] <- "-"
BWS_Sibs_Supplementary$Cluster2[BWS_Sibs_Supplementary$Cluster2 == 1] <- "-"
BWS_Sibs_Supplementary[BWS_Sibs_Supplementary == 22] <- "Hastings"
BWS_Sibs_Supplementary[BWS_Sibs_Supplementary == 29] <- "Poole"
BWS_Sibs_Supplementary[BWS_Sibs_Supplementary == 205] <- "Shawford"
BWS_Sibs_Supplementary[BWS_Sibs_Supplementary == 286] <- "Wrecclesham"
BWS_Sibs_Supplementary[BWS_Sibs_Supplementary == 328] <- "Walton"

# Change colnames
colnames(BWS_Sibs_Supplementary) <- c("Year", "Sample 1", "Trap Number 1", "Cluster 1 (if applicable)", "Sample 2", "Trap Number 2", "Cluster 2 (if applicable)", 
                                      "Probability of Sibship",	"Full/Half Sibship",	"Distance between traps (km)")

# Save file 
write.csv(BWS_Sibs_Supplementary, "Results/BWS_Sibling_Results.csv", row.names = FALSE)

# 4. Create datasets for subsequent analyses with this information

# All siblings from the 2017 dataset need to be removed
# Siblings from the 2018 dataset need to be removed to check the data/loci. 
# However these need to be retained for population structure analyses. 

# Create vector of which individuals to remove for locus descriptor analyses 
inds_to_remove_17 <- c("1285A", "1285H", "3129A", "3129F", "3129G", "3129J", "3129K", "3129Q") 

# Removal of 8 individuals 

inds_to_remove_18 <- c("0645H", "0645E", "1187A", "1715A", 
                       "2159D", "2448C", "2448G", "2448I", 
                       "5637G", "5912E")

# If performing this analysis with 0.8, then also remove"0645B", "3026B".

# Removal of 10 individuals 
# Although indications for 3253A and 5597B being full sibs, these have been retained as distance >135km - unlikely
# This leaves 184 individuals from 2018 (check this)

# Remove these individuals 
# Create function to remove 
`%nin%` <- Negate(`%in%`)

# Remove siblings 
bws_nosibs <- bws_data_50[bws_data_50$Sample %nin% c(inds_to_remove_17, inds_to_remove_18),]

# Save this 
write.csv(bws_nosibs, "Data/07_BWS_Data_50_NoSibs.csv", row.names = FALSE)

# To read this in
bws_nosibs <- read.csv("Data/07_BWS_Data_50_NoSibs.csv")

# These data need to be analysed separately as they are from different years 
# Subset these 

bws17_nosibs <- subset(bws_nosibs, bws_nosibs$Year == "2017") #nrow 97
bws18_nosibs <- subset(bws_nosibs, bws_nosibs$Year == "2018") #nrow 186

# Save files 

write.csv(bws17_nosibs, "Data/07a_BWS_Data_50_NoSibs_2017.csv", row.names = FALSE)
write.csv(bws18_nosibs, "Data/07b_BWS_Data_50_NoSibs_2018.csv", row.names = FALSE)

# 5. Data formatting for future analyses 
# These datasets will be used for the following: Genepop analyses and GenAlEx (Excel plugin)

# GenAlEx format requires: sample, population information, genotype data

# Remove unecessary columns

bws17_nosibs_genalex <- bws17_nosibs[ , ! names(bws17_nosibs) %in% c("Trap_ID", "Method",
                                                                     "cluster", "Year", 
                                                                     "Latitude", "Longitude", 
                                                                     "Vespula_vulgaris",
                                                                     "D3.15", "D3.15.1",
                                                                     "VMA4", "VMA4.1")]

bws18_nosibs_genalex <- bws18_nosibs[ , ! names(bws18_nosibs) %in% c("Trap_ID","Method",
                                                                      "Region", "Year", "Latitude",
                                                                      "Longitude", "Vespula_vulgaris", 
                                                                      "LIST2004", "LIST2004.1")]

# Save these for analysis with GenAlEx
write.csv(bws17_nosibs_genalex, "Data/08a_BWS_Data_50_NoSibs_2017_GenAlEx.csv", row.names = FALSE)
write.csv(bws18_nosibs_genalex, "Data/08b_BWS_Data_50_NoSibs_2018_GenAlEx.csv", row.names = FALSE)



