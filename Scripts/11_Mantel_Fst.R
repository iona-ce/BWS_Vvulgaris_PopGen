# Mantel test and Fst 

#Libraries 

library(hierfstat)
library(ade4)

## 1. Mantel test ####

# Create a list of objects to analyze
obj_list <- list(bws17_genind, bws18_regional_genind, bws18_genind_cl17, bws18_genind_cl22, bws18_genind_cl29,bws18_genind_cl57,
                 bws18_genind_cl62, bws18_genind_cl122, bws18_genind_cl193, bws18_genind_cl205, bws18_genind_cl286,
                 bws18_genind_cl328)

names(obj_list) <- c("National Scale", "Regional Scale", "Walkford Cluster", "Hastings Cluster", 
                     "Poole Cluster","Edinburgh Cluster", "Norwood Cluster", "Crawley Cluster", 
                     "Ferndown Cluster", "Shawford Cluster", "Wrecclesham Cluster", "Walton-on-Thames Cluster")

# Loop over objects and perform analysis
for(i in seq_along(obj_list)){
  
  # 1.b. Calculate genetic distance 
  Dgen <- dist(obj_list[[i]], method ="euclidian")
  
  # 1.c. Calculate geographic distances
  Dgeo <- dist(obj_list[[i]]@other$xy)
  Dgeo <- log10(Dgeo + 1)
  
  # 1.d. Perform Mantel test between two distance matrices
  ibd <- mantel.randtest(Dgen, Dgeo, nrepet = 10000)
  print(ibd)
  
  # 1.2. Plot Geographic and Genetic distances
  plot(Dgeo, Dgen,
       xlab = "Log10 (1 + Geographic distance)",
       ylab = "Genetic distance (Euclidian)")
  title(main = paste0(names(obj_list)[i]), cex.main = 3)
  abline(lm(Dgen ~ Dgeo), col = "red", lty = 2)
}

# Repeating Mantel without NI samples 
Dgen <- dist(bws17_genind_NI, method ="euclidian")

# 1.c. Calculate geographic distances
Dgeo <- dist(bws17_genind_NI@other$xy)
Dgeo <- log10(Dgeo + 1)

# 1.d. Perform Mantel test between two distance matrices
ibd <- mantel.randtest(Dgen, Dgeo, nrepet = 10000)
print(ibd)

# 2. F statistics using Nei's statistic ####

nei_fst_list <- list()

for(i in seq_along(obj_list)){
  
  # 2.a. Fst value 
  nei_fst <- pairwise.neifst(obj_list[[i]])
  
  #Transform negative values to 0 
  nei_fst[nei_fst < 0]  <- 0
  
  #Visualise table 
  nei_fst <- round(nei_fst, digits = 3)
  nei_fst_list[[i]] <- nei_fst
  names(nei_fst_list)[i] <- names(obj_list)[i]
  
}

# See results
nei_fst_list

# Optional step 
  
#Get maximum, minimum and mean Fst values 
# max(nei_fst, na.rm = TRUE)
# mean(nei_fst, na.rm = TRUE)
# median(nei_fst, na.rm = TRUE)

# 2.b. Bootstrap Fst 

# Long example; this seems to be only possible to do manually 

boot_regional <- boot.ppfst(bws18_regional_genind,
                nboot = 100,
                quant=c(0.025,0.975),
                diploid = TRUE)

#Create table of lower limit 
bootfst_regional_ll <- boot_regional$ll

#Transform negative values to 0
bootfst_regional_ll[bootfst_regional_ll < 0] <- 0

#Round to 3 digits 
bootfst_regional_ll <- round(bootfst_regional_ll, digits = 3)

#Create table of upper limit 
bootfst_regional_ul <- boot_regional$ul

#Transform negative values to 0
bootfst_regional_ul[bootfst_regional_ul < 0] <- 0

#Round to 3 digits 
bootfst_regional_ul <- round(bootfst_regional_ul, digits = 3)

# National scale
boot_national <- boot.ppfst(bws17_genind,
                            nboot = 100,
                            quant=c(0.025,0.975),
                            diploid = TRUE)

bootfst_national_ll <- boot_national$ll
bootfst_national_ll[bootfst_national_ll < 0] <- 0
bootfst_national_ll <- round(bootfst_national_ll, digits = 3)
bootfst_national_ul <- boot_national$ul
bootfst_national_ul[bootfst_national_ul < 0] <- 0
bootfst_national_ul <- round(bootfst_national_ul, digits = 3)

# Cluster 17
boot_cl17 <- boot.ppfst(bws18_genind_cl17,
                        nboot = 100,
                        quant=c(0.025,0.975),
                        diploid = TRUE)

bootfst_cl17_ll <- boot_cl17$ll
bootfst_cl17_ll[bootfst_cl17_ll < 0] <- 0
bootfst_cl17_ll <- round(bootfst_cl17_ll, digits = 3)
bootfst_cl17_ul <- boot_cl17$ul
bootfst_cl17_ul[bootfst_cl17_ul < 0] <- 0
bootfst_cl17_ul <- round(bootfst_cl17_ul, digits = 3)

# Cluster 22
boot_cl22 <- boot.ppfst(bws18_genind_cl22,
                        nboot = 100,
                        quant=c(0.025,0.975),
                        diploid = TRUE)

bootfst_cl17_ll <- boot_cl22$ll
bootfst_cl17_ll[bootfst_cl17_ll < 0] <- 0
bootfst_cl17_ll <- round(bootfst_cl17_ll, digits = 3)
bootfst_cl17_ul <- boot_cl22$ul
bootfst_cl17_ul[bootfst_cl17_ul < 0] <- 0
bootfst_cl17_ul <- round(bootfst_cl17_ul, digits = 3)

# Plot figures for manuscript 
layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), nrow = 2, byrow = TRUE))

Dgennat <- dist(bws17_genind, method ="euclidian")
Dgeonat <- dist(bws17_genind@other$xy)
Dgeonat <- log10(Dgeonat + 1)
plot(Dgeonat, Dgennat,
     xlab = "Log10 (1 + Geographic distance)",
     ylab = "Genetic distance (Euclidian)", cex.lab = 1.5)
title(main = "National", cex.main = 2)
abline(lm(Dgennat ~ Dgeonat), col = "red", lty = 2)

Dgenreg <- dist(bws18_regional_genind, method ="euclidian")
Dgeoreg <- dist(bws18_regional_genind@other$xy)
Dgeoreg <- log10(Dgeoreg + 1)
plot(Dgeoreg, Dgenreg,
     xlab = "Log10 (1 + Geographic distance)",
     ylab = "Genetic distance (Euclidian)", cex.lab = 1.5)
title(main = "Regional", cex.main = 2)
abline(lm(Dgenreg ~ Dgeoreg), col = "red", lty = 2)

Dgen1 <- dist(bws18_genind_cl22, method ="euclidian")
Dgeo1 <- dist(bws18_genind_cl22@other$xy)
Dgeo1 <- log10(Dgeo1 + 1)
plot(Dgeo1, Dgen1,
     xlab = "Log10 (1 + Geographic distance)",
     ylab = "Genetic distance (Euclidian)",
     cex.lab = 1.5)
title(main = "Hastings", cex.main = 2)
abline(lm(Dgen1 ~ Dgeo1), col = "red", lty = 2)

Dgen2 <- dist(bws18_genind_cl62, method ="euclidian")
Dgeo2 <- dist(bws18_genind_cl62@other$xy)
Dgeo2 <- log10(Dgeo2 + 1)
plot(Dgeo2, Dgen2,
     xlab = "Log10 (1 + Geographic distance)",
     ylab = "Genetic distance (Euclidian)", cex.lab = 1.5)
title(main = "Norwood", cex.main = 2)
abline(lm(Dgen2 ~ Dgeo2), col = "red", lty = 2)

Dgen3 <- dist(bws18_genind_cl205, method ="euclidian")
Dgeo3 <- dist(bws18_genind_cl205@other$xy)
Dgeo3 <- log10(Dgeo3 + 1)
plot(Dgeo3, Dgen3,
     xlab = "Log10 (1 + Geographic distance)",
     ylab = "Genetic distance (Euclidian)", cex.lab = 1.5)
title(main = "Shawford", cex.main = 2)
abline(lm(Dgen3 ~ Dgeo3), col = "red", lty = 2)




