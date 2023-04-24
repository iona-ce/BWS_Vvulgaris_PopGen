#### ANALYSING STRUCTURE OUTPUTS ####

# 1. Libraries ####

#Using pophelper http://www.royfrancis.com/pophelper/articles/index.html
#How to instgall this package pophelper
#install_github('royfrancis/pophelper')

library(devtools)
library(pophelper)
library(ggplot2)
library(gridExtra)

# Create slist 
slist <- readQ(c("Data/11_Structure/bws17_structurek1_results_f", "Data/11_Structure/bws17_structurek2_results_f", 
                 "Data/11_Structure/bws17_structurek3_results_f","Data/11_Structure/bws17_structurek4_results_f", 
                 "Data/11_Structure/bws17_structurek5_results_f", "Data/11_Structure/bws17_structurek6_results_f",
                 "Data/11_Structure/bws17_structurek1_2_results_f", "Data/11_Structure/bws17_structurek1_3_results_f",
                 "Data/11_Structure/bws17_structurek2_2_results_f", "Data/11_Structure/bws17_structurek2_3_results_f",
                 "Data/11_Structure/bws17_structurek3_2_results_f", "Data/11_Structure/bws17_structurek3_3_results_f",
                 "Data/11_Structure/bws17_structurek4_2_results_f", "Data/11_Structure/bws17_structurek4_3_results_f",
                 "Data/11_Structure/bws17_structurek5_2_results_f", "Data/11_Structure/bws17_structurek5_3_results_f", 
                 "Data/11_Structure/bws17_structurek6_2_results_f", "Data/11_Structure/bws17_structurek6_3_results_f"),
               filetype = "structure", indlabfromfile = T)

#Check that it has worked correctly 
head(slist[[1]])

## 2. TabulateQ ####

#Tabulate q to produce table of runs from qlist
tr1 <- tabulateQ(qlist=slist)
sr1 <- summariseQ(tr1)
summariseQ(tr1, writetable=TRUE, exportpath=getwd())

## 3. Evanno's Method to determine best ∆K ####

#Do Evanno's method with all the other runs 
evannoMethodStructure(data=sr1)

#Export plot 
evannoMethodStructure(data=sr1,exportplot=T,exportpath=getwd())

#Create plots
p <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p) #Take highest ∆K (here 28.8 for K = 2) 

## 4. Plotting the data ####

#Create pretty plots
slist <- alignK(slist)

# #Checking that the data are loaded correctly 
# p1 <- plotQ(slist,imgoutput="join",returnplot=T,exportplot=F,basesize=11)
# grid.arrange(p1$plot[[1]])

#Read in the data with the population info - the pop info needs to be in a 
#data frame format to be able to incorporate it into the graph
bws17_structure <- read.table("Data/11_Structure/bws17_structure.txt", stringsAsFactors=F)
bws17_structure$V2 <- as.character(bws17_structure$V2)

bws17_structure_pop <- as.data.frame(cbind(bws17_structure$V1, bws17_structure$V2))

#Group labels with one group label set
labk2<- data.frame(lab1=bws17_structure$V2)
labk2$lab1 <- as.character(labk2$lab1)
rownames(labk2) <- bws17_structure$V1
colnames(labk2) <- c("Pop")

#Rename 

labk2[labk2 == 1] <- "N. England"
labk2[labk2 == 2] <- "S. England"
labk2[labk2 == 3] <- "Scotland"
labk2[labk2 == 4] <- "E. England"
labk2[labk2 == 5] <- "W. Eng/Wales"
labk2[labk2 == 6] <- "N. Ireland"

regions <- c("S. England", "E. England","N. England", "W. Eng/Wales","Scotland", "N. Ireland")

#Test with only one K
# p1 <- plotQ(slist[13],returnplot=T,exportplot=F,basesize=11,
#             grplab=labk2,grplabsize=5,linesize=0.5,pointsize=6,
#             showdiv = TRUE, divsize = 1, ordergrp = TRUE,
#             sortind="all", subsetgrp = regions,
#             grplabangle = 90, grplabjust = 0.4)
# 
# grid.arrange(p1$plot[[1]])

#Choose which Ks to display (only need one per K of interest) - here for K=2, K=5, K=6

p1 <- plotQ(slist[c(4, 17)],imgoutput="join",sharedindlab = F,returnplot=T,exportplot=F,basesize=11,
            grplab=labk2, grplabsize=5, linesize=0.5, pointsize=6,
            showdiv = TRUE, divsize = 1, ordergrp = TRUE,
            splab = c("K=2", "K=6"), splabsize = 17,
            sortind="all", subsetgrp = regions,
            grplabangle = -90, grplabjust = 0,
            grplabpos = 0.8,
            panelratio = c(1.5,1), linepos = 0.9,
            titlelab = "National", showtitle = TRUE, titlehjust = 0.5, 
            titlesize = 14, titlecol = "black", titlespacer = 10, titleface = "bold",
            splabcol = "black")

grid.arrange(p1$plot[[1]])

# So the same for 2018 

slist2 <- readQ(c("Data/11_Structure/bws18_structurek1_results_f", "Data/11_Structure/bws18_structurek2_results_f",
                 "Data/11_Structure/bws18_structurek3_results_f", "Data/11_Structure/bws18_structurek4_results_f",
                 "Data/11_Structure/bws18_structurek5_results_f", "Data/11_Structure/bws18_structurek6_results_f",
                 "Data/11_Structure/bws18_structurek7_results_f", "Data/11_Structure/bws18_structurek8_results_f",
                 "Data/11_Structure/bws18_structurek9_results_f", "Data/11_Structure/bws18_structurek1_2_results_f",
                 "Data/11_Structure/bws18_structurek1_3_results_f", "Data/11_Structure/bws18_structurek2_2_results_f",
                 "Data/11_Structure/bws18_structurek2_3_results_f", "Data/11_Structure/bws18_structurek3_2_results_f",
                 "Data/11_Structure/bws18_structurek3_3_results_f", "Data/11_Structure/bws18_structurek4_2_results_f",
                 "Data/11_Structure/bws18_structurek4_3_results_f", "Data/11_Structure/bws18_structurek5_2_results_f",
                 "Data/11_Structure/bws18_structurek5_3_results_f", "Data/11_Structure/bws18_structurek6_2_results_f",
                 "Data/11_Structure/bws18_structurek6_3_results_f", "Data/11_Structure/bws18_structurek7_2_results_f",
                 "Data/11_Structure/bws18_structurek7_3_results_f", "Data/11_Structure/bws18_structurek8_2_results_f",
                 "Data/11_Structure/bws18_structurek8_3_results_f", "Data/11_Structure/bws18_structurek9_2_results_f",
                 "Data/11_Structure/bws18_structurek9_3_results_f"),
               filetype = "structure", indlabfromfile = T)

## 2. TabulateQ ####

#Tabulate q to produce table of runs from qlist
tr2 <- tabulateQ(qlist=slist2)
sr2 <- summariseQ(tr2)
summariseQ(tr2, writetable=TRUE, exportpath=getwd())

## 3. Evanno's Method to determine best ∆K ####

#Do Evanno's method with all the other runs 
evannoMethodStructure(data=sr2)

#Export plot 
evannoMethodStructure(data=sr2,exportplot=T,exportpath=getwd())

#Create plots
p_evanno2 <- evannoMethodStructure(data=sr2,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p_evanno2) #Take highest 119.46 for K=2

## 4. Plotting the data ####

#Create pretty plots
slist2 <- alignK(slist2)

# Make labels for plot 
bws18_pop<- data.frame(Pop=bws18_structure$cluster)
bws18_pop$Pop <- as.character(bws18_pop$Pop)
rownames(bws18_pop) <- bws18_structure$Sample

# Replace numbers with names
bws18_pop[bws18_pop == 17] <- "Walkford"
bws18_pop[bws18_pop == 22] <- "Hastings"
bws18_pop[bws18_pop == 29] <- "Poole"
bws18_pop[bws18_pop == 62] <- "Norwood"
bws18_pop[bws18_pop == 122] <- "Crawley"
bws18_pop[bws18_pop == 193] <- "Ferndown"
bws18_pop[bws18_pop == 205] <- "Shawford"
bws18_pop[bws18_pop == 286] <- "Wrecclesham"
bws18_pop[bws18_pop == 328] <- "Walton"

order_cluster <- c("Poole", "Ferndown", "Walkford", "Shawford", "Wrecclesham", "Walton", "Norwood", "Crawley", "Hastings")

# Visualise

p2 <- plotQ(slist2[c(4, 27)], imgoutput = "join", returnplot = TRUE, exportplot = FALSE, basesize = 11,
            sharedindlab = FALSE, sortind = "all",
            grplab = bws18_pop, grplabsize = 5, 
            linesize=0.5, pointsize=6, 
            ordergrp = TRUE,
            subsetgrp = order_cluster, 
            showdiv = TRUE, divsize = 1,
            splab = c("K=2", "K=9"), splabsize = 17,
            grplabangle = -90, grplabjust = 0,
            grplabpos = 0.8,
            panelratio = c(1.5,1), linepos = 0.9,
            titlelab = "Regional", showtitle = TRUE, titlehjust = 0.5, 
            titlesize = 14, titlecol = "black", titlespacer = 10, titleface = "bold",
            splabcol = "black")

grid.arrange(p2$plot[[1]])

grid.arrange(p1$plot[[1]], p2$plot[[1]], ncol = 2)

