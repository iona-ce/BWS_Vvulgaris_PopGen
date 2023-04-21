#### Data preparation for use with STRUCTURE ####

#Datasets for STRUCTURE analyses have following format:
# - Individual identifier (here trap ID)
# - Putative population identifier (here 'region'). Note that when doing this in terminal, must be numeric value
# - Genotype data
# - Missing values should be -9

# For 2017 ----

## Remove all unecessary columns
bws17_nosibs_structure <- bws17_nosibs[ , ! names(bws17_nosibs) %in% c("Trap_ID", "Method",
                                                                     "cluster", "Year", 
                                                                     "Latitude", "Longitude", 
                                                                     "Vespula_vulgaris",
                                                                     "D3.15", "D3.15.1",
                                                                     "VMA4", "VMA4.1")] #loci to remove

## Change region names to numeric values (if not does not work using terminal; 
# if using structure GUI, this does not need to be changed)
bws17_nosibs_structure[bws17_nosibs_structure == "NE"] <- 1
bws17_nosibs_structure[bws17_nosibs_structure == "SE"] <- 2
bws17_nosibs_structure[bws17_nosibs_structure == "SC"] <- 3
bws17_nosibs_structure[bws17_nosibs_structure == "EE"] <- 4
bws17_nosibs_structure[bws17_nosibs_structure == "WE"] <- 5
bws17_nosibs_structure[bws17_nosibs_structure == "NI"] <- 6

#Replace missing values from NA to -9
bws17_nosibs_structure[bws17_nosibs_structure == 0] <- -9

#Write in .txt format. 
write.table(bws17_nosibs_structure, "Data/11_Structure/bws17_structure.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# For 2018 ----
# Siblings are retained for BWS 2018 structure analyses; using data set bws18

bws18_structure <- bws18_nosibs[ , ! names(bws18_nosibs) %in% c("Trap_ID", "Method",
                                                  "Region", "Year", 
                                                  "Latitude", "Longitude", 
                                                  "Vespula_vulgaris", "count_na_row",
                                                  "LIST2004", "LIST2004.1")] #loci to remove

#Replace missing values from NA to -9
bws18_structure[is.na(bws18_structure)] <- -9

# Remove individuals to only have those that belong to a cluster 
bws18_structure <- bws18_structure[bws18_structure$cluster %in% c("17", "22", "29", "62", "122", "193", "205", "286", "328"),]

#Write in .txt format. 
write.table(bws18_structure, "Data/11_Structure/bws18_structure.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
