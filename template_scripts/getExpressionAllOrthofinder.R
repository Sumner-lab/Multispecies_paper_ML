# For importing resultsanscript estimates.
library(tximport)
library(tximportData)
# For differential expression analysis. 
library(edgeR)
# Import seqinr for editing fasta files.
library(seqinr)
# Load pheatmap.
library(pheatmap)
# For importing resultsanscript estimates.
library(tximport)
library(tximportData)
# For differential expression analysis. 
library(edgeR)
# Import seqinr for editing fasta files.
library(seqinr)
library(stringr)

# Remove all objects.
rm(list = ls())

# Set the working directory.
#setwd(paste("./DATA", sep = ""))

# Load previous data.
file.rdata <- paste("DATA/Orthofinder/VERSION_RUN/Orthogroups.RData")
load(file = file.rdata)

##### 

# Set parameters

# Use a hard-coded dispersion since we don't have replicates from which to perform 'accurate' estimates.
bcv <- 0.1 # square-root-dispersion - 0.1 for human data, 0.1 for genetically identical model organisms, 0.01 for technical replicates
# Set a cutoff p-value after multiple matrix.counts.log2.quantile.species.scaled.orthing correction.
p.value.cutoff = 0.05
# Set a log-fold cut-off.
log.fc.cutoff <- 1.5
## Load the orthologues
data.oma.orthologous.matrix <- data.orthologues.single.copy.genes
# Change the colnames.
#colnames(data.oma.orthologous.matrix) <- gsub("_transcriptome", "", colnames(data.oma.orthologous.matrix))

data.oma.orthologous.matrix <- data.oma.orthologous.matrix
species.names <- colnames(data.oma.orthologous.matrix)
species.names.short EATMYSHORTS

# data.oma.orthologous.matrix <- data.oma.orthologous.matrix[ , 3, drop = FALSE]
# species.names <- colnames(data.oma.orthologous.matrix)
# species.names.short <- c("BM")

# Number of species.
num.species = length(species.names)

#####

# Loop over each species, load the data, perform a basic DGE analyses, and store outputs for downstream analysis in lists.

# Store outputs in lists.
species.data.counts <- list()
species.data.tpm <- list()
species.data.counts.orth <- list()
species.data.tpm.orth <- list()
species.data.counts.species.scaled <- list()
species.data.tpm.species.scaled <- list()
species.data.counts.species.scaled.orth <- list()
species.data.tpm.species.scaled.orth <- list()

species.data.counts.log2.quantile <- list()
species.data.tpm.log2.quantile <- list()
species.data.counts.log2.quantile.orth <- list()
species.data.tpm.log2.quantile.orth <- list()
species.data.counts.log2.quantile.species.scaled <- list()
species.data.tpm.log2.quantile.species.scaled <- list()
species.data.counts.log2.quantile.species.scaled.orth <- list()
species.data.tpm.log2.quantile.species.scaled.orth <- list()

species.data.edgeR.results <- list()
species.genes.ge <- list()
species.genes.orth <- list()

#Chris addition
species.data.counts.orth.tpm <- list()
species.data.counts.orth.tpm.log2.quantile <- list()
species.data.counts.orth.tpm.log2.quantile.species.scaled <- list()
species.data.counts.orth.tpm.log2.species.scaled <- list()
species.data.counts.orth.tpm.log2 <- list()
species.data.counts.orth.log2.quantile.species.scaled <-  list()

for(i in 1:num.species){
  
  cat(species.names[i])
  
  # Load the trinity fasta.
  file.fasta.trinity.fnn <- paste("DATA/FIX_FOLDER/", species.names[i], "/Pool/Pool_trinity.fnn", sep = "")
  fasta.trinity.fnn <- read.fasta(file = file.fasta.trinity.fnn, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
  
  # Get the transcript example for each gene.
  transcript.names <- getName(fasta.trinity.fnn)
  gene.names <- gsub("_i.*", "", transcript.names)
  #Delete:
  fasta.trinity.fnn <- NULL
  # Transcript to gene map.
  tx2gene <- data.frame(transcripts = transcript.names, 
                        genes = gene.names, 
                        stringsAsFactors = FALSE)
  tx2gene <- tx2gene[order(tx2gene$transcripts), ]
  rownames(tx2gene) <- tx2gene$transcripts
  
  # Get the queen and worker sample names.
  queen.samples <- list.dirs(path = paste("DATA/FIX_FOLDER/", species.names[i] , "/", sep = ""), full.names = FALSE)
  queen.samples <- queen.samples[grepl("Queen|^Reproductive", queen.samples)]
  worker.samples <- list.dirs(path = paste("DATA/FIX_FOLDER/", species.names[i] , "/", sep = ""), full.names = FALSE)
  worker.samples <- worker.samples[grepl("Worker|^NonReproductive", worker.samples)]
  
  # Get the numbers of queens and workers.
  num.queens <- length(queen.samples)
  num.workers <- length(worker.samples)
  
  # Get the number of samples.
  num.samples <- num.queens + num.workers
  
  # Get sample names.
  sample.names <- c(queen.samples, worker.samples)
  sample.names.short <- paste(species.names.short[i], sample.names, sep = " ")
  
  # Get the groups.
  sample.groups <- c(rep("Queen", times = num.queens), rep("Worker", times = num.workers))
  
  # Get the rsem file list.
  queen.files <- list.files(path = paste("DATA/FIX_FOLDER/", species.names[i], "/", queen.samples, sep = ""), pattern = ".isoforms.results$", full.names = TRUE)
  worker.files <- list.files(path = paste("DATA/FIX_FOLDER/", species.names[i], "/", worker.samples, sep = ""), pattern = ".isoforms.results$", full.names = TRUE)
  
  all.files <- c(queen.files, worker.files)
  names(all.files) <- sample.names.short
  
  # Load the rsem files.
  # Note: txi.rsem$abundance = TPM scores.
  # Note: txi.rsem$counts = Raw counts.
  # Note: txi.rsem$length = Effective length.
  txi.rsem <- tximport(all.files, tx2gene = tx2gene, type = "rsem", txIn = TRUE, countsFromAbundance = "no")
  


  ## Get differential expression.

  # Get count data.
  data.counts <- txi.rsem$counts
  data.tpm <- txi.rsem$abundance

  #chris remove unused files:
  tx2gene <- NULL
  txi.rsem <- NULL
  
  # Get the non-NA orthologues.
  orths.present <- data.oma.orthologous.matrix[ , species.names[i], drop = FALSE]
  orths.present <- orths.present[(!is.na(orths.present[ , 1]) & (!orths.present[ , 1] == "")), , drop = FALSE]
  print (nrow(orths.present))
    
  # Create a data.frame to store the one-to-one orthologue expression values.
  data.counts.orth <- matrix(data = NA, nrow = dim(data.oma.orthologous.matrix)[1], ncol = num.samples)
  rownames(data.counts.orth) <- rownames(data.oma.orthologous.matrix)
  colnames(data.counts.orth) <- sample.names.short

  # Set the expression values for the non-NA orths.
  data.counts.orth[rownames(orths.present), ] <- data.counts[orths.present[ , 1], ]

  # Create a data.frame to store the one-to-one orthologue expression values.
  data.tpm.orth <- matrix(data = NA, nrow = dim(data.oma.orthologous.matrix)[1], ncol = num.samples)
  rownames(data.tpm.orth) <- rownames(data.oma.orthologous.matrix)
  colnames(data.tpm.orth) <- sample.names.short
  
  # Set the expression values for the non-NA orths.
  data.tpm.orth[rownames(orths.present), ] <- data.tpm[orths.present[ , 1], ]
  
  # Get within-species scaled data.
  data.counts.species.scaled <- t(apply(data.counts, 1, function(x){ x <- (x - mean(x)) / mean(x) }))
  data.tpm.species.scaled <- t(apply(data.tpm, 1, function(x){ x <- (x - mean(x)) / mean(x) }))
  
  # Create a data.frame to store the one-to-one orthologue expression values.
  data.counts.species.scaled.orth <- matrix(data = NA, nrow = dim(data.oma.orthologous.matrix)[1], ncol = num.samples)
  rownames(data.counts.species.scaled.orth) <- rownames(data.oma.orthologous.matrix)
  colnames(data.counts.species.scaled.orth) <- sample.names.short
  
  # Set the expression values for the non-NA orths.
  data.counts.species.scaled.orth[rownames(orths.present), ] <- data.counts.species.scaled[orths.present[ , 1], ]

  # Create a data.frame to store the one-to-one orthologue expression values.
  data.tpm.species.scaled.orth <- matrix(data = NA, nrow = dim(data.oma.orthologous.matrix)[1], ncol = num.samples)
  rownames(data.tpm.species.scaled.orth) <- rownames(data.oma.orthologous.matrix)
  colnames(data.tpm.species.scaled.orth) <- sample.names.short
  
  # Set the expression values for the non-NA orths.
  data.tpm.species.scaled.orth[rownames(orths.present), ] <- data.tpm.species.scaled[orths.present[ , 1], ]
  
  # Form a DGEList, modelling data on counts.
  data.edgeR <- DGEList(counts = data.counts, group = sample.groups)
  
  # Optionally filter out genes with low counts per million (CPM) except those that are known orthologues. Here we keep genes that are represented at least 1cpm reads in at least 1 sample.
  idx.keep <- rowSums(cpm(data.edgeR) > CPM_filter) >= floor((num.queens + num.workers) / 2)
  paste(sum(idx.keep))
  
  # We record the names of the transcripts and orthologues that fail this criterion so that we can set them to NA later.
  names.remove <- names(idx.keep[idx.keep == FALSE])
  names.remove.orth <- as.character(rownames(orths.present[match(names.remove[names.remove %in% orths.present[ , 1]], orths.present[ , 1]), , drop = FALSE]))
  
  # New DGEList filtered.
  data.edgeR <- DGEList(counts = data.counts[idx.keep, ], group = sample.groups)
  #write.table( data.edgeR, file=paste("DGE_result" , species.names[i], sep="_"), sep="\t")

  # Normalise for RNA composition effecounts using results mean of M-values (TMM). Note: can't do lfc shrink like in DESeq2.
  data.edgeR <- calcNormFactors(data.edgeR)
  
  # Estimate or hard code dispersion.
  # data.edgeR.genes <- estimateGLMCommonDisp(data.edgeR.genes, method = "deviance", robust = TRUE, subset = NULL)
  data.edgeR$common.dispersion <- bcv^2
  
  # Fit a negative binomial generalized log-linear model to the read counts for each gene. 
  design <- model.matrix(~sample.groups)
  data.edgeR.fit <- glmFit(data.edgeR, design)
  
  # Conduct genewise statistical matrix.counts.log2.quantile.species.scaled.orths for a given coefficient or contrast relative to a specified fold-change threshold.
  data.edgeR.results <- glmTreat(data.edgeR.fit, coef = 2, lfc = log.fc.cutoff)
  
  # Correct for multiple matrix.counts.log2.quantile.species.scaled.orthing.
  data.edgeR.results$table$fdr <- p.adjust(data.edgeR.results$table$PValue, method = "fdr")
  
  # Identify which genes are significantly differentially expressed from an edgeR fit object containing p-values and matrix.counts.log2.quantile.species.scaled.orth statistics.
  data.edgeR.results$de <- decideTestsDGE(data.edgeR.results)
  
  dir.create("DEGS/VERSION_RUN/")
  path_here = paste("DEGS/VERSION_RUN/", species.names[i], "deg.res", sep = "")
  write.table(data.edgeR.results$table, path_here)
  
  # Summarise the results.
  summary.absolute <- summary(data.edgeR.results$de)
  print(summary.absolute)
  summary.relative <- summary.absolute / sum(summary.absolute)
  print(summary.relative)

  #Print a volcano
  path_here = paste("DEGS/VERSION_RUN/Smear_", species.names[i], "deg.res.pdf", sep = "")
  pdf (path_here, height=6, width=6)
  plotSmear(data.edgeR.results,de.tags = rownames(data.edgeR.results$table)[which(data.edgeR.results$table$PValue<0.05)], smearWidth=0.1, pch=19, cex=0.4)
  dev.off()
  volcanoData2 <- cbind(data.edgeR.results$table$logFC, -log10(data.edgeR.results$table$fdr))
  colnames(volcanoData2) <- c("logFC", "fdr")
  path_hereN = paste("DEGS/VERSION_RUN/Volc_", species.names[i], "deg.res.pdf", sep = "")
  pdf (path_hereN, height=6, width=6)
  plot(volcanoData2, pch=19)
  dev.off()

  #Chris addition (to try to get orth before TPM:   ############
  scaling.factor <- colSums(data.counts.orth, , na.rm=T) # Often divided by  10^6 to get transcripts per million.
  data.counts.orth.tpm <- sweep(data.counts.orth, 2, scaling.factor, `/`)*1000000
  #data.counts.orth.tpm.log2.quantile <- normalizeBetweenArrays(log2(data.counts.orth.tpm+1), method = "quantile") Old quantile does not work, makes some genes reversed!!!
  data.counts.orth.tpm.log2.quantile <- data.counts.orth.tpm
  data.counts.orth.tpm.log2.quantile.species.scaled <- t(apply(data.counts.orth.tpm.log2.quantile, 1, function(x){ x <- (x - mean(x)) / mean(x) }))

  #Chris without quantile normalisation step.
  data.counts.orth.tpm.log2 <- log2(data.counts.orth.tpm+1)
  data.counts.orth.tpm.log2.species.scaled <- t(apply(data.counts.orth.tpm.log2, 1, function(x){ x <- (x - mean(x)) / mean(x) }))

  #Chris Addition (without tpm, chris version)
  #data.counts.orth.log2.quantile <- normalizeBetweenArrays(log2(data.counts.orth+1), method = "quantile") 
  data.counts.orth.log2.quantile <- data.counts.orth
  data.counts.orth.log2.quantile.species.scaled<-t(apply(data.counts.orth.log2.quantile, 1, function(x){ x <- (x - mean(x)) / mean(x) }))

  #Repeat for 1 to 1 orthogroups present.
  #First remove NA rows. Else DESeq complains.
  data.counts.orth <- data.counts.orth[rowSums(is.na(data.counts.orth)) == 0,]

  #Then repeat exactly the same as before.
  data.edgeR2 <- DGEList(counts = data.counts.orth, group = sample.groups) 
   # Optionally filter out genes with low counts per million (CPM) except those that are known orthologues. Here we keep genes that are represented at least 1cpm reads in at least 1 sample.
  idx.keep2 <- rowSums(cpm(data.edgeR2) > 1) >= floor((num.queens + num.workers) / 2)
  paste(sum(idx.keep2))
  
  # We record the names of the transcripts and orthologues that fail this criterion so that we can set them to NA later.
  names.remove2 <- names(idx.keep2[idx.keep2 == FALSE])
  names.remove2.orth <- as.character(rownames(orths.present[match(names.remove2[names.remove2 %in% orths.present[ , 1]], orths.present[ , 1]), , drop = FALSE]))
  
  # New DGEList filtered.
  data.edgeR2 <- DGEList(counts = data.counts.orth[idx.keep2, ], group = sample.groups)
  #write.table( data.edgeR2, file=paste("DGE_result" , species.names[i], sep="_"), sep="\t")

  # Normalise for RNA composition effecounts using results mean of M-values (TMM). Note: can't do lfc shrink like in DESeq2.
  data.edgeR2 <- calcNormFactors(data.edgeR2)
  
  # Estimate or hard code dispersion.
  # data.edgeR2.genes <- estimateGLMCommonDisp(data.edgeR2.genes, method = "deviance", robust = TRUE, subset = NULL)
  data.edgeR2$common.dispersion <- bcv^2
  
  # Fit a negative binomial generalized log-linear model to the read counts for each gene. 
  design <- model.matrix(~sample.groups)
  data.edgeR2.fit <- glmFit(data.edgeR2, design)
  
  # Conduct genewise statistical matrix.counts.log2.quantile.species.scaled.orths for a given coefficient or contrast relative to a specified fold-change threshold.
  data.edgeR2.results <- glmTreat(data.edgeR2.fit, coef = 2, lfc = log.fc.cutoff)
  
  # Correct for multiple matrix.counts.log2.quantile.species.scaled.orthing.
  data.edgeR2.results$table$fdr <- p.adjust(data.edgeR2.results$table$PValue, method = "fdr")
  
  # Identify which genes are significantly differentially expressed from an edgeR fit object containing p-values and matrix.counts.log2.quantile.species.scaled.orth statistics.
  data.edgeR2.results$de <- decideTestsDGE(data.edgeR2.results)

  path_here = paste("DEGS/VERSION_RUN/", species.names[i], "deg.res.orth", sep = "")
  write.table(data.edgeR2.results$table, path_here)

  # Summarise the results.
  summary.absolute <- summary(data.edgeR2.results$de)
  print(summary.absolute)
  summary.relative <- summary.absolute / sum(summary.absolute)
  print(summary.relative)

  #Print a volcano
  path_here = paste("DEGS/VERSION_RUN/Smear_", species.names[i], "deg.res.orth.pdf", sep = "")
  pdf (path_here, height=6, width=6)
  plotSmear(data.edgeR2.results,de.tags = rownames(data.edgeR2.results$table)[which(data.edgeR2.results$table$PValue<0.05)], smearWidth=0.1, pch=19, cex=0.4)
  dev.off()
  volcanoData2 <- cbind(data.edgeR2.results$table$logFC, -log10(data.edgeR2.results$table$fdr))
  colnames(volcanoData2) <- c("logFC", "fdr")
  path_hereN = paste("DEGS/VERSION_RUN/Volc_", species.names[i], "deg.res.orth.pdf", sep = "")
  pdf (path_hereN, height=6, width=6)
  plot(volcanoData2, pch=19)
  dev.off()
  
  ## Get the abundance data for plotting.
  
  # Quantile normalise.
  data.counts.log2.quantile <- normalizeBetweenArrays(log2(data.counts+1), method = "quantile")
  data.tpm.log2.quantile <- normalizeBetweenArrays(log2(data.tpm+1), method = "quantile")
  
  # Create a data.frame to store the one-to-one orthologue expression values.
  data.counts.log2.quantile.orth <- matrix(data = NA, nrow = dim(data.oma.orthologous.matrix)[1], ncol = num.samples)
  rownames(data.counts.log2.quantile.orth) <- rownames(data.oma.orthologous.matrix)
  colnames(data.counts.log2.quantile.orth) <- sample.names.short
  
  # Set the expression values for the non-NA orths.
  data.counts.log2.quantile.orth[rownames(orths.present), ] <- data.counts.log2.quantile[orths.present[ , 1], ]
  
  # Create a data.frame to store the one-to-one orthologue expression values.
  data.tpm.log2.quantile.orth <- matrix(data = NA, nrow = dim(data.oma.orthologous.matrix)[1], ncol = num.samples)
  rownames(data.tpm.log2.quantile.orth) <- rownames(data.oma.orthologous.matrix)
  colnames(data.tpm.log2.quantile.orth) <- sample.names.short
  
  # Set the expression values for the non-NA orths.
  data.tpm.log2.quantile.orth[rownames(orths.present), ] <- data.tpm.log2.quantile[orths.present[ , 1], ]
  
  ## Get within-species scaled data.
  data.counts.log2.quantile.species.scaled <- t(apply(data.counts.log2.quantile, 1, function(x){ x <- (x - mean(x)) / mean(x) }))
  data.tpm.log2.quantile.species.scaled <- t(apply(data.tpm.log2.quantile, 1, function(x){ x <- (x - mean(x)) / mean(x) }))
  
  # Create a data.frame to store the one-to-one orthologue expression values.
  data.counts.log2.quantile.species.scaled.orth <- matrix(data = NA, nrow = dim(data.oma.orthologous.matrix)[1], ncol = num.samples)
  rownames(data.counts.log2.quantile.species.scaled.orth) <- rownames(data.oma.orthologous.matrix)
  colnames(data.counts.log2.quantile.species.scaled.orth) <- sample.names.short

  # Set the expression values for the non-NA orths.
  data.counts.log2.quantile.species.scaled.orth[rownames(orths.present), ] <- data.counts.log2.quantile.species.scaled[orths.present[ , 1], ]

  # Create a data.frame to store the one-to-one orthologue expression values.
  data.tpm.log2.quantile.species.scaled.orth <- matrix(data = NA, nrow = dim(data.oma.orthologous.matrix)[1], ncol = num.samples)
  rownames(data.tpm.log2.quantile.species.scaled.orth) <- rownames(data.oma.orthologous.matrix)
  colnames(data.tpm.log2.quantile.species.scaled.orth) <- sample.names.short
  
  # Set the expression values for the non-NA orths.
  data.tpm.log2.quantile.species.scaled.orth[rownames(orths.present), ] <- data.tpm.log2.quantile.species.scaled[orths.present[ , 1], ]
  
  # Add all data to lists.
  species.data.counts[[i]] <- data.counts
  species.data.tpm[[i]] <- data.tpm
  species.data.counts.orth[[i]] <- data.counts.orth
  species.data.tpm.orth[[i]] <- data.tpm.orth
  species.data.counts.species.scaled[[i]] <- data.counts.species.scaled
  species.data.tpm.species.scaled[[i]] <- data.tpm.species.scaled
  species.data.counts.species.scaled.orth[[i]] <- data.counts.species.scaled.orth
  species.data.tpm.species.scaled.orth[[i]] <- data.tpm.species.scaled.orth
  
  species.data.counts.log2.quantile[[i]] <- data.counts.log2.quantile
  species.data.tpm.log2.quantile[[i]] <- data.tpm.log2.quantile
  species.data.counts.log2.quantile.orth[[i]] <- data.counts.log2.quantile.orth
  species.data.tpm.log2.quantile.orth[[i]] <- data.tpm.log2.quantile.orth
  species.data.counts.log2.quantile.species.scaled[[i]] <- data.counts.log2.quantile.species.scaled
  species.data.tpm.log2.quantile.species.scaled[[i]] <- data.tpm.log2.quantile.species.scaled
  species.data.counts.log2.quantile.species.scaled.orth[[i]] <- data.counts.log2.quantile.species.scaled.orth
  species.data.tpm.log2.quantile.species.scaled.orth[[i]] <- data.tpm.log2.quantile.species.scaled.orth
  
  species.data.edgeR.results[[i]] <- data.edgeR.results$table
  species.genes.ge[[i]] <- rownames(data.edgeR.results$table[data.edgeR.results$table$fdr < 0.05, ])
  species.genes.orth[[i]] <- data.oma.orthologous.matrix[!is.na(data.oma.orthologous.matrix[ , species.names[i]]), species.names[i], drop = FALSE]
  
  #Chris added:
  species.data.counts.orth.tpm.log2.quantile.species.scaled[[i]] <- data.counts.orth.tpm.log2.quantile.species.scaled
  species.data.counts.orth.tpm[[i]] <- data.counts.orth.tpm
  species.data.counts.orth.tpm.log2.species.scaled[[i]] <- data.counts.orth.tpm.log2.species.scaled
  species.data.counts.orth.tpm.log2[[i]] <- data.counts.orth.tpm.log2

  species.data.counts.orth.log2.quantile.species.scaled[[i]] <- data.counts.orth.log2.quantile.species.scaled

  # Set CPM fails to NA.
  species.data.counts[[i]][names.remove, ] <- NA
  species.data.tpm[[i]][names.remove, ] <- NA
  species.data.counts.orth[[i]][names.remove.orth, ] <- NA
  species.data.tpm.orth[[i]][names.remove.orth, ] <- NA
  species.data.counts.species.scaled[[i]][names.remove, ] <- NA
  species.data.tpm.species.scaled[[i]][names.remove, ] <- NA
  species.data.counts.species.scaled.orth[[i]][names.remove.orth, ] <- NA
  species.data.tpm.species.scaled.orth[[i]][names.remove.orth, ] <- NA
  
  species.data.counts.log2.quantile[[i]][names.remove, ] <- NA
  species.data.tpm.log2.quantile[[i]][names.remove, ] <- NA
  species.data.counts.log2.quantile.orth[[i]][names.remove.orth, ] <- NA
  species.data.tpm.log2.quantile.orth[[i]][names.remove.orth, ] <- NA
  species.data.counts.log2.quantile.species.scaled[[i]][names.remove, ] <- NA
  species.data.tpm.log2.quantile.species.scaled[[i]][names.remove, ] <- NA
  species.data.counts.log2.quantile.species.scaled.orth[[i]][names.remove.orth, ] <- NA
  species.data.tpm.log2.quantile.species.scaled.orth[[i]][names.remove.orth, ] <- NA
  
  #Chris added:Set CPM fails to NA.:
  species.data.counts.orth.tpm.log2.quantile.species.scaled[[i]][names.remove.orth, ] <- NA
  species.data.counts.orth.tpm[[i]][names.remove.orth, ] <- NA
  species.data.counts.orth.tpm.log2.species.scaled[[i]][names.remove.orth, ] <- NA
  species.data.counts.orth.tpm.log2[[i]][names.remove.orth, ] <- NA
  species.data.counts.orth.log2.quantile.species.scaled[[i]][names.remove.orth, ] <- NA

}

# Add names to the list.
names(species.data.counts) <- species.names
names(species.data.tpm) <- species.names
names(species.data.counts.orth) <- species.names
names(species.data.tpm.orth) <- species.names
names(species.data.counts.species.scaled) <- species.names
names(species.data.tpm.species.scaled) <- species.names
names(species.data.counts.species.scaled.orth) <- species.names
names(species.data.tpm.species.scaled.orth) <- species.names

names(species.data.counts.log2.quantile) <- species.names
names(species.data.tpm.log2.quantile) <- species.names
names(species.data.counts.log2.quantile.orth) <- species.names
names(species.data.tpm.log2.quantile.orth) <- species.names
names(species.data.counts.log2.quantile.species.scaled) <- species.names
names(species.data.tpm.log2.quantile.species.scaled) <- species.names
names(species.data.counts.log2.quantile.species.scaled.orth) <- species.names
names(species.data.tpm.log2.quantile.species.scaled.orth) <- species.names

names(species.data.edgeR.results) <- species.names
names(species.genes.ge) <- species.names
names(species.genes.orth) <- species.names

#Chris added:
names(species.data.counts.orth.tpm.log2.quantile.species.scaled) <- species.names
names(species.data.counts.orth.tpm) <- species.names
names(species.data.counts.orth.tpm.log2.species.scaled) <- species.names 
names(species.data.counts.orth.tpm.log2) <- species.names 
names(species.data.counts.orth.log2.quantile.species.scaled) <- species.names 

#chris remove unused files (save memory):
species.data.tpm.log2.quantile.species.scaled <- NULL
species.data.tpm.species.scaled <- NULL
species.data.tpm.log2.quantile <- NULL
species.data.tpm <- NULL
species.data.counts.species.scaled <- NULL
species.data.counts.log2.quantile.species.scaled <- NULL
species.data.counts.log2.quantile <- NULL
species.data.counts  <- NULL

#####

# Output the data.

# Save the workspace so it can be loaded from start without having to repeat all these steps.
file.rdata <- paste("DATA/Orthofinder/VERSION_RUN/OrthogroupsExpression.RData")
save(list = ls(), file = file.rdata)
