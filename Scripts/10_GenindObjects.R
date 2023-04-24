###Transforming the data into a Genind object

# 1. Load libraries ----
basic_libraries <- c("testthat",
                     "ape",
                     "pegas",
                     "ggplot2",
                     "adegenet",
                     "dplyr",
                     "poppr",
                     "geosphere",
                     "sp",
                     "vegan",
                     "hierfstat",
                     "devtools",
                     "pegas",
                     "genepop",
                     "hierfstat",
                     "adegenet",
                     "ade4")

for (lib in basic_libraries) {
  if (require(package = lib, character.only = TRUE)) {
    print("Successful")
  } else {
    print("Installing")
    install.packages(lib)
    library(lib, character.only = TRUE )
  }
}

## 2018 data
# 2. Formatting data  ----

# 2018 data ----

# Cluster level data 
# Remove columns that are not needed

bws18_genind_table <- bws18[ , ! names(bws18) %in% c("Method", "Region", "Year", 
                                                     "Code", "Postcode", 
                                                     "Vespula_vulgaris", "count_na_row",
                                                     "LIST2004", "LIST2004.1")] #loci to remove

bws18_genind_table <- bws18_genind_table[bws18_genind_table$cluster %in% c("17", "22", "29", "57", "62", "122", "193", "205", "286", "328"),]


# Function to rename colnames in dataframe 
renaming_locus_name <- function(my_dataframe){
  
  # stop the function if its input is neither a data frame nor a matrix
  stopifnot(is.data.frame(my_dataframe) || is.matrix(my_dataframe))
  
  # obtain locus names
  locus_name_vec <- unique(gsub(x = colnames(my_dataframe)[4:29],
                                pattern = "\\.1$",
                                replacement = ""))
  
  # remove full stops
  locus_name_vec <- gsub(x = locus_name_vec,
                         pattern = "\\.",
                         replacement = "")
  
  # create a empty vectore for column names
  column_names_subset <- c()
  
  # for each element in the vector
  for(locus_name in locus_name_vec){
    
    # repeat the name of the locus with "_1" and "_2"
    column_names_subset <- append(column_names_subset,
                                  c(paste(locus_name, "_1", sep = ""),
                                    paste(locus_name, "_2", sep = "")))
  }
  
  
  # Rename columns 
  colnames(my_dataframe) <- c("Sample",
                              "Trap_ID",
                              "Pop",
                              column_names_subset,
                              "Lat",
                              "Long")
  
  return(my_dataframe)
}

bws18_genind_table <- renaming_locus_name(bws18_genind_table)

# Modify format to fit adegenet genind object

# Change loci to single column separated by a "/"
bws18_reformated <- as.data.frame(cbind(paste(bws18_genind_table$LIST2007_1,
                                              bws18_genind_table$LIST2007_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$LIST2003_1,
                                              bws18_genind_table$LIST2003_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$LIST2013_1,
                                              bws18_genind_table$LIST2013_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$LIST2018_1,
                                              bws18_genind_table$LIST2018_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$VMA3_1,
                                              bws18_genind_table$VMA3_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$LIST2001_1,
                                              bws18_genind_table$LIST2001_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$LIST2011_1,
                                              bws18_genind_table$LIST2011_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$VMA6_1,
                                              bws18_genind_table$VMA6_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$R1169_1,
                                              bws18_genind_table$R1169_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$LIST2018_1,
                                              bws18_genind_table$LIST2018_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$Rufa19_1,
                                              bws18_genind_table$Rufa19_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$D315_1,
                                              bws18_genind_table$D315_2,
                                              sep = "/"),
                                        paste(bws18_genind_table$VMA4_1,
                                              bws18_genind_table$VMA4_2,
                                              sep = "/")),
                                  stringsAsFactors = FALSE)

# obtain locus names, without the _1 or _2
locus_name_vec <- unique(gsub(x = colnames(bws18_genind_table)[4:29],
                              pattern = "\\_[1,2]$",
                              replacement = ""))

# Rename columns and row with loci and sample names 
colnames(bws18_reformated) <- locus_name_vec
rownames(bws18_reformated) <- bws18_reformated$Sample

# Create genind object
bws18_genind <- df2genind(X      = bws18_reformated,
                          sep    = "/",
                          ind.names = bws18_genind_table$Sample,
                          loc.names = locus_name_vec,
                          pop    = bws18_genind_table$Pop,
                          NA.char = "0",
                          ploidy = 2,
                          type   = "codom",
                          strata = data.frame(bws18_genind_table$Trap_ID)) #adding in subpop info

# Add xy data x = long, y = lat 
# Create lat/long dataframe 
xydat <- as.data.frame(cbind(bws18_genind_table$Long,
                             bws18_genind_table$Lat))
# Insert row and column names 
rownames(xydat) <- bws18_genind_table$Sample
colnames(xydat) <- c("x", "y")

# Add xy dataframe to genind object 
bws18_genind@other$xy <- xydat

# Subset to cluster data ----

# Aim: subset data by population 

bws18_genind_cl17 <- popsub(bws18_genind,
                            sublist = "17")

bws18_genind_cl17@pop <- bws18_genind_cl17@strata$bws18_genind_table.Trap_ID

bws18_genind_cl22 <- popsub(bws18_genind,
                            sublist = "22")

bws18_genind_cl22@pop <- bws18_genind_cl22@strata$bws18_genind_table.Trap_ID

bws18_genind_cl29 <- popsub(bws18_genind,
                            sublist = "29")

bws18_genind_cl29@pop <- bws18_genind_cl29@strata$bws18_genind_table.Trap_ID

bws18_genind_cl57 <- popsub(bws18_genind,
                            sublist = "57")

bws18_genind_cl57@pop <- bws18_genind_cl57@strata$bws18_genind_table.Trap_ID

bws18_genind_cl62 <- popsub(bws18_genind,
                            sublist = "62")

bws18_genind_cl62@pop <- bws18_genind_cl62@strata$bws18_genind_table.Trap_ID

bws18_genind_cl122 <- popsub(bws18_genind,
                             sublist = "122")

bws18_genind_cl122@pop <- bws18_genind_cl122@strata$bws18_genind_table.Trap_ID

bws18_genind_cl193 <- popsub(bws18_genind,
                             sublist = "193")

bws18_genind_cl193@pop <- bws18_genind_cl193@strata$bws18_genind_table.Trap_ID

bws18_genind_cl205 <- popsub(bws18_genind,
                             sublist = "205")

bws18_genind_cl205@pop <- bws18_genind_cl205@strata$bws18_genind_table.Trap_ID

bws18_genind_cl286 <- popsub(bws18_genind,
                             sublist = "286")

bws18_genind_cl286@pop <- bws18_genind_cl286@strata$bws18_genind_table.Trap_ID

bws18_genind_cl328 <- popsub(bws18_genind,
                             sublist = "328")

bws18_genind_cl328@pop <- bws18_genind_cl328@strata$bws18_genind_table.Trap_ID

# Regional level data (without siblings) ----

bws18_genind_table_regional <- bws18_nosibs[ , ! names(bws18_nosibs) %in% c("Method", "Region", "Year", 
                                                                            "Code", "Postcode", 
                                                                            "Vespula_vulgaris", "count_na_row",
                                                                            "LIST2004", "LIST2004.1")] #loci to remove

bws18_genind_table_regional <- bws18_genind_table_regional[bws18_genind_table_regional$cluster %in% c("17", "22", "29", "62", "122", "193", "205", "286", "328"),]

bws18_genind_table_regional <- renaming_locus_name(bws18_genind_table_regional)

bws18_regional_reformated <- as.data.frame(cbind(paste(bws18_genind_table_regional$LIST2007_1,
                                                       bws18_genind_table_regional$LIST2007_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$LIST2003_1,
                                              bws18_genind_table_regional$LIST2003_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$LIST2013_1,
                                              bws18_genind_table_regional$LIST2013_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$LIST2018_1,
                                              bws18_genind_table_regional$LIST2018_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$VMA3_1,
                                              bws18_genind_table_regional$VMA3_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$LIST2001_1,
                                              bws18_genind_table_regional$LIST2001_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$LIST2011_1,
                                              bws18_genind_table_regional$LIST2011_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$VMA6_1,
                                              bws18_genind_table_regional$VMA6_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$R1169_1,
                                              bws18_genind_table_regional$R1169_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$LIST2018_1,
                                              bws18_genind_table_regional$LIST2018_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$Rufa19_1,
                                              bws18_genind_table_regional$Rufa19_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$D315_1,
                                              bws18_genind_table_regional$D315_2,
                                              sep = "/"),
                                        paste(bws18_genind_table_regional$VMA4_1,
                                              bws18_genind_table_regional$VMA4_2,
                                              sep = "/")),
                                  stringsAsFactors = FALSE)

# obtain locus names, without the _1 or _2
locus_name_vec <- unique(gsub(x = colnames(bws18_genind_table_regional)[4:29],
                              pattern = "\\_[1,2]$",
                              replacement = ""))

# Rename columns and row with loci and sample names 
colnames(bws18_regional_reformated) <- locus_name_vec
rownames(bws18_regional_reformated) <- bws18_genind_table_regional$Sample

# Create genind object
bws18_regional_genind <- df2genind(X      = bws18_regional_reformated,
                          sep    = "/",
                          ind.names = bws18_genind_table_regional$Sample,
                          loc.names = locus_name_vec,
                          pop    = bws18_genind_table_regional$Pop,
                          NA.char = "0",
                          ploidy = 2,
                          type   = "codom",
                          strata = data.frame(bws18_genind_table_regional$Trap_ID)) #adding in subpop info

# Add xy data x = long, y = lat 
# Create lat/long dataframe 
xydat <- as.data.frame(cbind(bws18_genind_table_regional$Long,
                             bws18_genind_table_regional$Lat))
# Insert row and column names 
rownames(xydat) <- bws18_genind_table_regional$Sample
colnames(xydat) <- c("x", "y")

# Add xy dataframe to genind object 
bws18_regional_genind@other$xy <- xydat


# Save the genind objects 
save(bws18_regional_genind, file = "Data/12_GenindObjects/12_BWS_50_NoSibs_2018_Regional_Genind.RData")

save(bws18_genind_cl17, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl17_Genind.RData")
save(bws18_genind_cl22, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl22_Genind.RData")
save(bws18_genind_cl29, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl29_Genind.RData")
save(bws18_genind_cl62, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl62_Genind.RData")
save(bws18_genind_cl122, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl122_Genind.RData")
save(bws18_genind_cl193, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl193_Genind.RData")
save(bws18_genind_cl205, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl205_Genind.RData")
save(bws18_genind_cl286, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl286_Genind.RData")
save(bws18_genind_cl328, file = "Data/12_GenindObjects/12_BWS_50_WithSibs_2018_cl328_Genind.RData")

# 2017 data ----

# Cleaning data 
bws17_genind_table <- bws17_nosibs[ , ! names(bws17_nosibs) %in% c("Trap_ID", "Method", "cluster", 
                                                                   "Year",
                                                                   "D3.15", "D3.15.1", "VMA4", "VMA4.1",
                                                                   "Vespula_vulgaris", "count_na_row")]

# Renaming locus name for 2017 data
renaming_locus_name <- function(my_dataframe){
  
  # stop the function if its input is neither a data frame nor a matrix
  stopifnot(is.data.frame(my_dataframe) || is.matrix(my_dataframe))
  
  # obtain locus names
  locus_name_vec <- unique(gsub(x = colnames(my_dataframe)[3:26],
                                pattern = "\\.1$",
                                replacement = ""))
  
  # remove full stops
  locus_name_vec <- gsub(x = locus_name_vec,
                         pattern = "\\.",
                         replacement = "")
  
  # create a empty vectore for column names
  column_names_subset <- c()
  
  # for each element in the vector
  for(locus_name in locus_name_vec){
    
    # repeat the name of the locus with "_1" and "_2"
    column_names_subset <- append(column_names_subset,
                                  c(paste(locus_name, "_1", sep = ""),
                                    paste(locus_name, "_2", sep = "")))
  }
  
  # Rename columns 
  colnames(my_dataframe) <- c("Sample",
                              "Pop",
                              column_names_subset,
                              "Lat",
                              "Long")
  
  return(my_dataframe)
}

#Apply this function to BWS dataframe 
bws17_genind_table <- renaming_locus_name(bws17_genind_table)

# Change loci to single column separated by a "/" (as per genind format)
#Here it is necessary to make sure that all loci are accounted for in the correct dataset 

bws17_genind_reformated <- as.data.frame(cbind(paste(bws17_genind_table$LIST2007_1,
                                                     bws17_genind_table$LIST2007_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$LIST2003_1,
                                                     bws17_genind_table$LIST2003_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$LIST2013_1,
                                                     bws17_genind_table$LIST2013_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$LIST2018_1,
                                                     bws17_genind_table$LIST2018_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$VMA3_1,
                                                     bws17_genind_table$VMA3_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$LIST2004_1,
                                                     bws17_genind_table$LIST2004_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$LIST2001_1,
                                                     bws17_genind_table$LIST2001_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$LIST2011_1,
                                                     bws17_genind_table$LIST2011_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$VMA6_1,
                                                     bws17_genind_table$VMA6_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$R1169_1,
                                                     bws17_genind_table$R1169_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$LIST2017_1,
                                                     bws17_genind_table$LIST2017_2,
                                                     sep = "/"),
                                               paste(bws17_genind_table$Rufa19_1,
                                                     bws17_genind_table$Rufa19_2,
                                                     sep = "/")),
                                         stringsAsFactors = FALSE)

# Obtain locus names, without the _1 or _2 (as these are now in one single column)
locus_name_vec <- unique(gsub(x = colnames(bws17_genind_table)[3:26],
                              pattern = "\\_[1,2]$",
                              replacement = ""))

# Rename columns and row with loci and sample names 
colnames(bws17_genind_reformated) <- locus_name_vec
rownames(bws17_genind_reformated) <- bws17_genind_table$Sample

## 3.b. Create Genind Object with LatLong data ####

bws17_genind <- df2genind(X      = bws17_genind_reformated,
                          sep    = "/",
                          ind.names = bws17_genind_table$Sample,
                          loc.names = locus_name_vec,
                          pop    = bws17_genind_table$Pop,
                          NA.char = "0",
                          ploidy = 2,
                          type   = "codom")

# Add xy data x = long, y = lat 
# Create lat/long dataframe 
xydat <- as.data.frame(cbind(bws17_genind_table$Long,
                             bws17_genind_table$Lat))
# Insert row and column names 
rownames(xydat) <- bws17_genind_table$Sample
colnames(xydat) <- c("x", "y")

# Add xy dataframe to genind object 
bws17_genind@other$xy <- xydat

# Save this 
save(bws17_genind, file = "Data/12_GenindObjects/12_BWS_50_2017_National_Genind.RData")

# 2017 without Northern Ireland data 
bws17_genind_NI <- popsub(bws17_genind,
                             sublist = c("EE", "SE", "NE", "SC", "WE"))
