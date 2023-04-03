# Preparing data for use in COLONY

# Load libraries 
library(dplyr)

# Using file "Data/02_BWS_DATA_50.csv" (made in 01_TidyingData_50Complete.R script)

# 1. 2017 data (national scale) ----

#Remove columns that are not required in COLONY analyses 
bws17_colony <- bws17[ , ! names(bws17) %in% c( "Trap_ID", "Method", "Region", "cluster", 
                                         "Year", "Latitude",
                                         "Longitude", "Vespula_vulgaris")]



#Write this as a .txt file for later analyses 
write.table(bws17_colony, "Data/03_COLONY_BWS17.txt", col.names = FALSE, 
            row.names = FALSE, quote = FALSE)

#Note that Colony folder should have COLONY executable

## 2. 2018 data (regional scale) ----
#Repeat the same for 2018 data 

# Re-order columns
# Note that here, LIST2004 is removed due to low amplification success 
bws18_colony <- bws18[ , ! names(bws18) %in% c( "Trap_ID", "Method", "Region", "cluster", 
                                                "Year", "Latitude",
                                                "Longitude", "Vespula_vulgaris", "LIST2004", "LIST2004.1")]

#Write this as a .txt file for later analyses 
write.table(bws18_colony, 
            "Data/04_COLONY_BWS18.txt", col.names = FALSE, 
            row.names = FALSE, quote = FALSE)

# Note that to be used for COLONY analyses, these files need to be manually edited 
# to add in necessary information (see COLONY manual)

# Download COLONY here: 
# https://www.zsl.org/science/software/colony


