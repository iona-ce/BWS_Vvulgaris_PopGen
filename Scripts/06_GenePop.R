# Analysing data with GenPop to look for

# Null alleles (NA) 
# Deviation from Hardy-Weinberg (HWE)
# Linkage Disequilibrium (LD)
# Inbreeding Coefficient (Fis)

# Load libraries 
library(genepop)

# 1. Genetic tests ----
# Get basic information

# For 2017
basic_info("Data/09a_BWS_Data_50_NoSibs_2017_GenePop.txt", 
           outputFile = "Results/GenePop/01a_BWS_Data_50_NoSibs_2017_BasicInfo",
           verbose = interactive())

# For 2018
basic_info("Data/09b_BWS_Data_50_NoSibs_2018_GenePop.txt", 
           outputFile = "Results/GenePop/01b_BWS_Data_50_NoSibs_2018_BasicInfo",
           verbose = interactive())

# Deviations from Hardy-Weinberg and Fis

# For 2017
test_HW(
  "Data/09a_BWS_Data_50_NoSibs_2017_GenePop.txt",
  which = "Proba",
  outputFile = "Results/GenePop/02a_BWS_Data_50_NoSibs_2017_HardyWeinberg",
  settingsFile = "",
  enumeration = FALSE,
  dememorization = 10000,
  batches = 20,
  iterations = 5000,
  verbose = interactive()
)

# For 2018
test_HW(
  "Data/09b_BWS_Data_50_NoSibs_2018_GenePop.txt",
  which = "Proba",
  outputFile = "Results/GenePop/02b_BWS_Data_50_NoSibs_2018_HardyWeinberg",
  settingsFile = "",
  enumeration = FALSE,
  dememorization = 10000,
  batches = 20,
  iterations = 10000,
  verbose = interactive()
)

# Linkage disequilibrium 
# For 2017
test_LD(
  "Data/09a_BWS_Data_50_NoSibs_2017_GenePop.txt",
  outputFile = "Results/GenePop/03b_BWS_Data_50_NoSibs_2017_LinkageDisequilibrium",
  settingsFile = "",
  dememorization = 10000,
  batches = 100,
  iterations = 5000,
  verbose = interactive()
)

# For 2018
test_LD(
  "Data/09b_BWS_Data_50_NoSibs_2018_GenePop.txt",
  outputFile = "Results/GenePop/03b_BWS_Data_50_NoSibs_2018_LinkageDisequilibrium",
  settingsFile = "",
  dememorization = 10000,
  batches = 100,
  iterations = 5000,
  verbose = interactive()
)

# Null alleles 

# For 2017
nulls(
  "Data/09a_BWS_Data_50_NoSibs_2017_GenePop.txt",
  outputFile = "Results/GenePop/04a_BWS_Data_50_NoSibs_2017_NullAlleles",
  settingsFile = "",
  nullAlleleMethod = "",
  CIcoverage = 0.95,
  verbose = interactive()
)

# For 2018
nulls(
  "Data/09b_BWS_Data_50_NoSibs_2018_GenePop.txt",
  outputFile = "Results/GenePop/04b_BWS_Data_50_NoSibs_2018_NullAlleles",
  settingsFile = "",
  nullAlleleMethod = "",
  CIcoverage = 0.95,
  verbose = interactive()
)

# 2. Bonferroni for HWE ----

HW_all_result_vec_2017 <- readLines("Results/GenePop/02a_BWS_Data_50_NoSibs_2017_HardyWeinberg")[31:42]
HW_all_result_vec_2018 <- readLines("Results/GenePop/02b_BWS_Data_50_NoSibs_2018_HardyWeinberg")[31:43]

HW_all_result_vec_updated_2017 <- gsub(pattern = "^[A-Z].*    ",
                                  x = HW_all_result_vec_2017,
                                  replacement = "")

HW_all_result_vec_updated_2018 <- gsub(pattern = "^[A-Z].*    ",
                                       x = HW_all_result_vec_2018,
                                       replacement = "")

HW_all_p_values_vec_2017 <- as.numeric(gsub(pattern = " .*",
                                       x = HW_all_result_vec_updated_2017,
                                       replacement = ""))

HW_all_p_values_vec_2018 <- as.numeric(gsub(pattern = " .*",
                                            x = HW_all_result_vec_updated_2018,
                                            replacement = ""))

# Apply Bonferroni correction on P-value from HWE
HW_all_adj_p_values_vec_2017 <- p.adjust(HW_all_p_values_vec_2017,
                                    method = "bonferroni")

HW_all_adj_p_values_vec_2018 <- p.adjust(HW_all_p_values_vec_2018,
                                         method = "bonferroni")

# Create table with these results

# Create vectors of loci names 
# For 2017
loci_2017 <- c("LIST2007", "LIST2003", "LIST2013", "LIST2018", "VMA3", "LIST2004",
               "LIST2001", "LIST2011", "VMA6", "R1.169", "LIST2017", "Rufa19")

# For 2018
loci_2018 <- c("LIST2007", "LIST2003", "LIST2013", "LIST2018", "VMA3", "LIST2001", 
               "LIST2011", "VMA6", "R1.169", "LIST2017", "Rufa19", "D3.15", "VMA4")

# Create table
# 2017
HWE_2017 <- data.frame(loci_2017, HW_all_p_values_vec_2017, HW_all_adj_p_values_vec_2017)
colnames(HWE_2017) <- c("Locus", "p-value", "adjusted p-value (Bonferroni)")

# 2018
HWE_2018 <- data.frame(loci_2018, HW_all_p_values_vec_2018, HW_all_adj_p_values_vec_2018)
colnames(HWE_2018) <- c("Locus", "p-value", "adjusted p-value (Bonferroni)")

# Save these files 
write.csv(HWE_2017, "Results/GenePop/05a_BWS_Data_50_NoSibs_2017_HardyWeinberg_Bonferroni.csv", row.names =  FALSE)
write.csv(HWE_2018, "Results/GenePop/05b_BWS_Data_50_NoSibs_2018_HardyWeinberg_Bonferroni.csv", row.names = FALSE)

