# Revisions 
# METHOD COMPARISON ----

# Aim: Assess effect of liquid/number of days that the trap was left out for/extraction method on DNA extraction success 


# Load libraries
library(stringr)

# Read necessary datasets 
trapinfo_2017 <- read.csv("Data/14_Revisions/14a_TrapInfo2017.csv")
bws_data <- read.csv("Data/01_BWS_DATA_COMPLETE.csv")
trapinfo_2018 <- read.csv("Data/14_Revisions/14b_TrapInfo2018.csv")

# Add NA rows 
bws_data$Region[is.na(bws_data$Region)] <- 1
bws_data$cluster[is.na(bws_data$cluster)] <- 1

bws_data$count_na_row <- (rowSums(is.na(bws_data)))/2 

# Subset bws_data to 2017 and 2018
bws_data_2017 <- subset(bws_data, bws_data$Year == 2017)
bws_data_2018 <- subset(bws_data, bws_data$Year == 2018)

#new_date <- as.Date(date, format = "%d/%m/%Y")
trapinfo_2017$inDate <- as.Date(trapinfo_2017$inDate, format = "%d/%m/%Y")
trapinfo_2017$outDate <- as.Date(trapinfo_2017$outDate, format = "%d/%m/%Y")

trapinfo_2018$inDate <- as.Date(trapinfo_2018$inDate, format = "%d/%m/%Y")
trapinfo_2018$outDate <- as.Date(trapinfo_2018$outDate, format = "%d/%m/%Y")

# Get difference in time that traps were set out for 
trapinfo_2017$days_total <- difftime(time1 = trapinfo_2017$inDate, time2 = trapinfo_2017$outDate, , units = "days")
trapinfo_2018$days_total <- difftime(time1 = trapinfo_2018$inDate, time2 = trapinfo_2018$outDate, , units = "days")

# Remove word "days"
trapinfo_2017$days_total <- gsub(" days", "", trapinfo_2017$days_total)
trapinfo_2018$days_total <- gsub(" days", "", trapinfo_2018$days_total)

# Change to as.numeric 
trapinfo_2017$days_total <- as.numeric(trapinfo_2017$days_total)
trapinfo_2018$days_total <- as.numeric(trapinfo_2018$days_total)

trapinfo_2017 <- trapinfo_2017[ , names(trapinfo_2017) %in% c( "Trap_ID", "days_total", "liquid")]
trapinfo_2018 <- trapinfo_2018[ , names(trapinfo_2018) %in% c( "ID", "days_total", "type")]

trapinfo_2018 <- trapinfo_2018[, c("ID", "days_total", "type")]

# Change colnames 
colnames(trapinfo_2018) <- colnames(trapinfo_2017)

# Match the two together 
bws_data_2017$Liquid <- trapinfo_2017[(match(bws_data_2017$Trap_ID, trapinfo_2017$Trap_ID)),]$liquid
bws_data_2018$Liquid <- trapinfo_2018[(match(bws_data_2018$Trap_ID, trapinfo_2018$ID)),]$liquid

bws_data_2017$Days <- trapinfo_2017[(match(bws_data_2017$Trap_ID, trapinfo_2017$Trap_ID)),]$days_total
bws_data_2018$Days <- trapinfo_2018[(match(bws_data_2018$Trap_ID, trapinfo_2018$ID)),]$days_total

# Merge these two dataframes 
bws_data <- rbind(bws_data_2017, bws_data_2018)

# Perform linear regression on the data
model <- lm(count_na_row ~ Method + Trap_ID + Liquid, data = bws_data)
summary(model)
model1 <- lm(count_na_row ~ Method + Days + Trap_ID + Liquid , data = bws_data)
summary(model1)

# Based on the results, it appears that the main impacts on PCR success depends on the method,
# and a small amount by the number of days that the sample was contained in the trap 

# Looking at different methods 

bws_dneasy <- subset(bws_data, bws_data$Method == "Dneasy")
bws_chelex <- subset(bws_data, bws_data$Method == "Chelex")

mean(bws_dneasy$count_na_row)
mean(bws_chelex$count_na_row)
t.test(bws_dneasy$count_na_row, bws_chelex$count_na_row) # p = 3.537e-10; 
# DNeasy is significantly more efficient than Chelex on BWS samples 

# CLARIFICATION OF SAMPLE SIZES ----

# Total number of samples selected for the study

nrow(bws_data) #383. 
# This does not include an additional 10 DNeasy 2017 samples that we failed to extract DNA from. 
# The total number of samples therefore is 383 + 10 = 393.

# Total number of samples depending on year 

# 2017
nrow(subset(bws_data, bws_data$Year == 2017)) #127 
nrow(bws_data_2017)
# 127 samples from 2017, + the 10 mentioned previously = 137

# 2018
nrow(subset(bws_data, bws_data$Year == 2018)) #256 
nrow(bws_data_2018)
# 256 samples from 2018 

# Total number of samples depending on method 
# DNeasy
nrow(subset(bws_data, bws_data$Method == "Dneasy")) #81
# + the 10 samples mentioned previously = 91 

# Chelex 
nrow(subset(bws_data, bws_data$Method == "Chelex")) #302

# Total number of randomly selected samples from 2017
# This is equivalent to the DNeasy samples, therefore = 91 

# Total number of samples that were selected from under-sampled areas 
# These are the 2017 Chelex samples 
nrow(subset(bws_data, bws_data$Year == 2017 & bws_data$Method == "Chelex")) #46

# Checking the number of Chelex 2018 samples 
nrow(subset(bws_data, bws_data$Year == 2018 & bws_data$Method == "Chelex")) #256

nrow(subset(bws_data, bws_data$Year == 2017 & bws_data$Method == "Dneasy")) #81

# Understanding success rates 
nrow(subset(bws_data_50, bws_data_50$Method == "Chelex")) #222
nrow(bws_chelex) - nrow(subset(bws_data_50, bws_data_50$Method == "Chelex")) #80 
nrow(subset(bws_data_50, bws_data_50$Method == "Dneasy")) #79

nrow(subset(bws_data_50, bws_data_50$Method == "Chelex" & bws_data_50$Year == 2017)) #26
nrow(subset(bws_data_50, bws_data_50$Method == "Chelex" & bws_data_50$Year == 2018)) #196

# Number of samples removed 
nrow(subset(bws_data, bws_data$Method == "Chelex" & bws_data$Year == 2017)) - nrow(subset(bws_data_50, bws_data_50$Method == "Chelex" & bws_data_50$Year == 2017))
# 20 samples removed from 2017 with Chelex method 
nrow(subset(bws_data, bws_data$Method == "Chelex" & bws_data$Year == 2018)) - nrow(subset(bws_data_50, bws_data_50$Method == "Chelex" & bws_data_50$Year == 2018))
# 60 samples removed from 2018 with Chelex method 
nrow(subset(bws_data, bws_data$Method == "Dneasy" & bws_data$Year == 2017)) - nrow(subset(bws_data_50, bws_data_50$Method == "Dneasy" & bws_data_50$Year == 2017))
# 2 samples removed from 2017 with DNeasy method 

# Final sample sizes with different methods 
nrow(subset(bws_data_50, bws_data_50$Method == "Chelex")) #222
nrow(subset(bws_data_50, bws_data_50$Method == "Dneasy")) #79

# Calculating percentages after removing samples with low amplification success
bws_50_chelex <- subset(bws_data_50, bws_data_50$Method == "Chelex")
bws_50_dneasy <- subset(bws_data_50, bws_data_50$Method == "Dneasy")

# Check amplification success across Chelex 50 samples 
1 - colSums(bws_50_chelex != 0)/nrow(bws_50_chelex)

# Check amplification success across all samples 
colSums(!is.na(bws17))/nrow(bws17)







# Delving deeper into the 2018 samples 

# Change cluster numbers if they are not equal to 17, 22, 29, 57, 62, 122, 193, 205, 286, 328
library(dplyr)
bws_data <- bws_data %>%
  mutate(cluster = ifelse(cluster %in% c(17, 22, 29, 57, 62, 122, 193, 205, 286, 328), cluster, 1))




