# Revisions 

# Aim: Assess effect of liquid/number of days that the trap was left out for/extraction method on DNA extraction success 

# Adding trap liquid data for manuscript revisions 

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

