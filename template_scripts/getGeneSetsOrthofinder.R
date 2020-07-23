#This script takes the output of Orthofinder of the species you wish to run. Both the model training set and the targets.
#Finds Orthofinder/VERSION_RUN/SingleCopyOrthogroups.txt Orthofinder/VERSION_RUN/Orthogroups.csv Orthofinder/VERSION_RUN/Orthogroups.GeneCount.csv
#You need to set the lines below

# Remove all objects.
rm(list = ls())

# Species all.
HEREspecies.all
# Training species
HEREspecies.train
# Testing species
HEREspecies.test
# Set the working directory.
setwd(paste(".", sep = ""))

library(tximport) # For importing resultsanscript estimates.
library(tximportData)
library(edgeR) # For differential expression analysis. 
library(seqinr) # Import seqinr for editing fasta files.
library(stringr)


# Number of test species.
num.species.test <- length(species.test)
# Number of train species.
num.species.train <- length(species.train)

# Get the single copy orthologues.
file.single.copy <- paste("DATA/Orthofinder/VERSION_RUN/SingleCopyOrthogroups.txt", sep = "")
data.single.copy <- read.table(file.single.copy)

# Get the orthologue matrix.
file.orthologues <- paste("DATA/Orthofinder/VERSION_RUN/Orthogroups.csv", sep = "")
data.orthologues <- read.csv2(file.orthologues, sep = "\t", colClasses = "character")
rownames(data.orthologues) <- data.orthologues$X
data.orthologues <- data.orthologues[ , -1]
#colnames(data.orthologues) <- gsub(".contig", "", colnames(data.orthologues))

# Get the orthologue matrix counts.
file.orthologues.counts <- paste("DATA/Orthofinder/VERSION_RUN/Orthogroups.GeneCount.csv", sep = "")
data.orthologues.counts <- read.csv2(file.orthologues.counts, sep = "\t", colClasses = "character")
rownames(data.orthologues.counts) <- data.orthologues.counts$X
data.orthologues.counts <- data.orthologues.counts[ , -1]
#data.orthologues.counts <- data.orthologues.counts[ , -length(colnames(data.orthologues.counts))]
#colnames(data.orthologues.counts) <- gsub(".contig", "", colnames(data.orthologues.counts))

# Get the orthologue matrix counts for train and test
data.orthologues.counts.train <- data.orthologues.counts[ , species.train]
data.orthologues.counts.test <- data.orthologues.counts[ , species.test]

# Get the single copy orthologues from the training data.
data.orthologues.single.copy <- data.orthologues[rowSums(as.matrix(data.orthologues.counts.train == 1)) == dim(data.orthologues.counts.train)[2], ]
data.orthologues.counts.single.copy <- rownames(data.orthologues.single.copy)

# Get the gene ids.
data.orthologues.single.copy.genes <- data.orthologues.single.copy

for(i in 1:dim(data.orthologues)[1]){
  for(j in 1:dim(data.orthologues)[2]){
    data.orthologues.single.copy.genes[i, j] <- gsub('__.*', '', data.orthologues.single.copy[i, j])
    #data.orthologues.single.copy.genes[i, j] <- gsub('_i.*', '', data.orthologues.single.copy[i, j])
  }
}

##### Delete big files, that are not used later:
data.orthologues<- NULL
data.orthologues.single.copy <- NULL

# Output the data.

# Save the workspace so it can be loaded from start without having to repeat all these steps.
file.rdata <- paste("DATA/Orthofinder/VERSION_RUN/Orthogroups.RData")
#savehistory ("DATA/Orthofinder/VERSION_RUN/History_step1.txt")
save(list = ls(), file = file.rdata)


