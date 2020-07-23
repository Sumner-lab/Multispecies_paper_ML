# For importing resultsanscript estimates.
library(tximport)
library(tximportData)

# For differential expression analysis. 
library(edgeR)

# Import seqinr for editing fasta files.
library(seqinr)

# Load pheatmap.
library(pheatmap)

library(stringr)

# Interface to the libsvm library in C.
library(e1071)
library(probsvm)

# Remove all objects.
rm(list = ls())
#dev.off()

# Set the working directory.
setwd(paste(".", sep = ""))

# Load previous data.
file.rdata <- paste("DATA/Orthofinder/VERSION_RUN/OrthogroupsExpression.RData")
load(file = file.rdata)

#Creating raw tpm (Michael style) output matrix(Chris addition)
CREATECOUNTAGAIN2
matrix.data.tpm.orth.filter <- matrix.data.tpm.orth[rowSums(!is.na(matrix.data.tpm.orth)) == dim(matrix.data.tpm.orth)[2], ]
matrix.data.tpm.orth.filter.scale <- scale(t(scale(t(matrix.data.tpm.orth.filter))))
matrix.data.tpm.orth.filter.scale.train <- matrix.data.tpm.orth.filter.scale[ , -GREPFITHERE , colnames(matrix.data.tpm.orth.filter.scale))]

#Creating normalised tpm (Michael style) matrix(Chris addition)
CREATETPMHERE
matrix.tpm.log2.quantile.species.scaled.orth.filter <- matrix.tpm.log2.quantile.species.scaled.orth[rowSums(!is.na(matrix.tpm.log2.quantile.species.scaled.orth)) == dim(matrix.tpm.log2.quantile.species.scaled.orth)[2], ]
matrix.tpm.log2.quantile.species.scaled.orth.filter.scale <- scale(t(scale(t(matrix.tpm.log2.quantile.species.scaled.orth.filter))))
matrix.tpm.log2.quantile.species.scaled.orth.filter.scale.train <- matrix.tpm.log2.quantile.species.scaled.orth.filter.scale[ , -GREPFITHERE , colnames(matrix.tpm.log2.quantile.species.scaled.orth.filter.scale))]

#creating file for Chris real TPMs.
CREATE_REAL_TPM_HERE
matrix_real.tpm.filter <- matrix_real.tpm[rowSums(!is.na(matrix_real.tpm)) == dim(matrix_real.tpm)[2], ]
matrix_real.tpm.filter.scale <- scale(t(scale(t(matrix_real.tpm.filter))))
matrix_real.tpm.filter.scale.train <- matrix_real.tpm.filter.scale[ , -GREPFITHERE , colnames(matrix_real.tpm.filter.scale))]

#creating file for Chris real TPM without quatile nos, full normalisation.
CREATE_UNQUANTILED
matrix_real.UNQ.filter <- matrix_real.UNQ[rowSums(!is.na(matrix_real.UNQ)) == dim(matrix_real.UNQ)[2], ]
matrix_real.UNQ.filter.scale <- scale(t(scale(t(matrix_real.UNQ.filter))))
matrix_real.UNQ.filter.scale.train <- matrix_real.UNQ.filter.scale[ , -GREPFITHERE , colnames(matrix_real.UNQ.filter.scale))]

#creating file for Chris real TPM without quatile nos log.
CREATE_LOG2_TPM_UNQUANT
matrix_real.TPM_UNQUANT.filter <- matrix_real.TPM_UNQUANT[rowSums(!is.na(matrix_real.TPM_UNQUANT)) == dim(matrix_real.TPM_UNQUANT)[2], ]
matrix_real.TPM_UNQUANT.filter.scale <- scale(t(scale(t(matrix_real.TPM_UNQUANT.filter))))
matrix_real.TPM_UNQUANT.filter.scale.train <- matrix_real.TPM_UNQUANT.filter.scale[ , -GREPFITHERE , colnames(matrix_real.TPM_UNQUANT.filter.scale))]

#Creating orth first, tpm adjusted values (Chris addition)
CREATEorthTPMHERE
matrix.data.counts.orth.tpm.log2.quantile.species.scaled.filter <- matrix.data.counts.orth.tpm.log2.quantile.species.scaled[rowSums(!is.na(matrix.data.counts.orth.tpm.log2.quantile.species.scaled)) == dim(matrix.data.counts.orth.tpm.log2.quantile.species.scaled)[2], ]
matrix.data.counts.orth.tpm.log2.quantile.species.scaled.filter.scale <- scale(t(scale(t(matrix.data.counts.orth.tpm.log2.quantile.species.scaled.filter))))
matrix.data.counts.orth.tpm.log2.quantile.species.scaled.filter.scale.train <- matrix.data.counts.orth.tpm.log2.quantile.species.scaled.filter.scale[ , -GREPFITHERE , colnames(matrix.data.counts.orth.tpm.log2.quantile.species.scaled.filter.scale))]


# Creating the output. Michael style.
CREATECOUNTHERE
matrix.counts.log2.quantile.species.scaled.orth.filter <- matrix.counts.log2.quantile.species.scaled.orth[rowSums(!is.na(matrix.counts.log2.quantile.species.scaled.orth)) == dim(matrix.counts.log2.quantile.species.scaled.orth)[2], ]

##### Choose whether to normalise between arrays at this point.. ####

# matrix.counts.log2.quantile.species.scaled.orth.quantile <- normalizeBetweenArrays(matrix.counts.log2.quantile.species.scaled.orth.filter, method = "quantile")
matrix.counts.log2.quantile.species.scaled.orth.quantile <- matrix.counts.log2.quantile.species.scaled.orth.filter
matrix.counts.log2.quantile.species.scaled.orth.quantile.scale <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile))))

matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train <- matrix.counts.log2.quantile.species.scaled.orth.quantile.scale[ , -GREPFITHERE , colnames(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale))]
matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test <- matrix.counts.log2.quantile.species.scaled.orth.quantile.scale[ , GREPFITHERE , colnames(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale))]
matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.together <- cbind(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train, matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test)

dir.create("FIGURES/Figure_of_Classifiers/VERSION_RUN/", showWarnings = FALSE)

caste <- grepl("Worker|NonReproductive", colnames(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train))
caste <- unlist(lapply(caste, function(x) { 
  if(x){ 
    x <- 0
  } else{ 
    x <- 1
  }
}))
caste <- as.numeric(caste)
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Boxplot_filter.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
boxplot(matrix.counts.log2.quantile.species.scaled.orth.filter)
dev.off()
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Boxplot_quantile.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
boxplot(matrix.counts.log2.quantile.species.scaled.orth.quantile)
dev.off()
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Boxplot_scale.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
boxplot(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale)
dev.off()
b.coeffs <- c()
p.coeffs <- c()

for(i in 1:dim(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[ , ])[1]){
  
  data <- data.frame(caste = caste, expr = as.numeric(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[i, ]))
  regression <- lm(caste ~ expr, data)
  s <- summary(regression)
  b.coeffs[i] <- regression$coefficients[2]
  p.coeffs[i] <- s$coefficients[8]
  # boxplot(expr ~ caste, data)
  
}

top.n <- 50 # dim(data.orthologues.tpm.exprs.log2.quantile.scale)[1] # max is dim(data.orthologues.tpm.exprs.log2.quantile.scale)[1]
cutoff <- abs(b.coeffs)[order(abs(b.coeffs), decreasing = TRUE)][top.n]
b.coeffs.sig <- abs(b.coeffs) >= cutoff
print(sum(b.coeffs.sig))
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/B_P_coefficients.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
plot(b.coeffs, p.coeffs)
p.coeffs.sig <- p.coeffs < 0.05
q.coeffs <- p.adjust(p.coeffs)
q.coeffs.sig <- q.coeffs < 0.05
print(length(p.coeffs.sig))
print(sum(p.coeffs.sig))
print(sum(q.coeffs.sig))
dev.off()

#Set up output tables:
play<-as.data.frame(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train, h=T)
play$p.coeffs<-p.coeffs
play$q.coeffs<-q.coeffs
withExpr<-merge(matrix.data.tpm.orth.filter, play, by="row.names")
write.table(withExpr, "FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_table_with_significance.tsv", sep="\t", quote=F)

# Train data
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_Heatmap.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[])
dev.off()

data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[]))


pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_PCA.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
plot(data.pca$x[ , 1], data.pca$x[ , 2])
text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)
dev.off()

pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_cnt_PCA.pdf"), width = 1.5 * 3.54, height = 1.5 * 3.54)
plot(data.pca$x[ , 1], data.pca$x[ , 2])
text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)
dev.off()

data.pca_raw <- prcomp(t(matrix.data.tpm.orth.filter.scale.train[]))
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_raw_PCA.pdf"), width = 1.5 * 3.54, height = 1.5 * 3.54)
plot(data.pca_raw$x[ , 1], data.pca_raw$x[ , 2])
text(data.pca_raw$x[ , 1], data.pca_raw$x[ , 2], labels=rownames(data.pca_raw$x), cex= 0.7, pos=3)
dev.off()
#Set up output tables:
play_raw<-as.data.frame(matrix.data.tpm.orth.filter.scale.train, h=T)
write.table(play_raw, "FIGURES/Figure_of_Classifiers/Primitive_eusocial_merged_2comb/Train_table_raw.tsv", sep="\t", quote=F)

data.pca_tpm <- prcomp(t(matrix.tpm.log2.quantile.species.scaled.orth.filter.scale.train[]))
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_tpm_PCA.pdf"), width = 1.5 * 3.54, height = 1.5 * 3.54)
plot(data.pca_tpm$x[ , 1], data.pca_tpm$x[ , 2])
text(data.pca_tpm$x[ , 1], data.pca_tpm$x[ , 2], labels=rownames(data.pca_tpm$x), cex= 0.7, pos=3)
dev.off()
#Set up output tables:
play_tpm<-as.data.frame(matrix.tpm.log2.quantile.species.scaled.orth.filter.scale.train, h=T)
write.table(play_tpm, "FIGURES/Figure_of_Classifiers/Primitive_eusocial_merged_2comb/Train_table_tpm.tsv", sep="\t", quote=F)


#plot(data.pca$x[ , 1], data.pca$x[ , 3])
#text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

#plot(data.pca$x[ , 2], data.pca$x[ , 3])
#text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)
#dev.off()

pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_Heatmap.p_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[p.coeffs.sig, ])
dev.off()

data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[p.coeffs.sig, ]))


pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_PCA.p_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
plot(data.pca$x[ , 1], data.pca$x[ , 2])
text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)
dev.off()

#plot(data.pca$x[ , 1], data.pca$x[ , 3])
#text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

#plot(data.pca$x[ , 2], data.pca$x[ , 3])
#text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_Heatmap.q_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[q.coeffs.sig, ])
dev.off()

data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[q.coeffs.sig, ]))

pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Train_PCA.q_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
plot(data.pca$x[ , 1], data.pca$x[ , 2])
text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)
dev.off()

#plot(data.pca$x[ , 1], data.pca$x[ , 3])
#text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

#plot(data.pca$x[ , 2], data.pca$x[ , 3])
#text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

# Test data
if (length(species.test) >= 2){
  pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Test_Heatmap.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
  pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[])
  data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[]))
  dev.off()

  pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Test_PCA.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
  plot(data.pca$x[ , 1], data.pca$x[ , 2])
  text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)
  dev.off()

  #plot(data.pca$x[ , 1], data.pca$x[ , 3])
  #text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

  #plot(data.pca$x[ , 2], data.pca$x[ , 3])
  #text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

  pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Test_Heatmap.p_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
  pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[p.coeffs.sig, ])
  data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[p.coeffs.sig, ]))
  dev.off()

  pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Test_PCA.p_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
  plot(data.pca$x[ , 1], data.pca$x[ , 2])
  text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)
  dev.off()

  #plot(data.pca$x[ , 1], data.pca$x[ , 3])
  #text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

  #plot(data.pca$x[ , 2], data.pca$x[ , 3])
  #text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

  pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Test_Heatmap.q_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
  pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[q.coeffs.sig, ])
  data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[q.coeffs.sig, ]))
  dev.off()

  pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Test_PCA.q_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
  plot(data.pca$x[ , 1], data.pca$x[ , 2])
  text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)
  dev.off()

  #plot(data.pca$x[ , 1], data.pca$x[ , 3])
  #text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

  #plot(data.pca$x[ , 2], data.pca$x[ , 3])
  #text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)
}


# Together
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Together_Heatmap.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.together[])
dev.off()

pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Together_PCA.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.together[]))
plot(data.pca$x[ , 1], data.pca$x[ , 2])
text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)
dev.off()

#plot(data.pca$x[ , 1], data.pca$x[ , 3])
#text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

#plot(data.pca$x[ , 2], data.pca$x[ , 3])
#text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)


pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Together_Heatmap.p_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.together[p.coeffs.sig, ])
dev.off()
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Together_PCA.p_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.together[p.coeffs.sig, ]))
plot(data.pca$x[ , 1], data.pca$x[ , 2])
text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)

#plot(data.pca$x[ , 1], data.pca$x[ , 3])
#text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

#plot(data.pca$x[ , 2], data.pca$x[ , 3])
#text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)
dev.off()

pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Together_Heatmap.q_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
pheatmap(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.together[q.coeffs.sig, ])
dev.off()
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Together_PCA.q_coeff.pdf"), width = 3 * 1.5 * 3.54, height = 1.5 * 3.54)
data.pca <- prcomp(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.together[q.coeffs.sig, ]))
plot(data.pca$x[ , 1], data.pca$x[ , 2])
text(data.pca$x[ , 1], data.pca$x[ , 2], labels=rownames(data.pca$x), cex= 0.7, pos=3)

#plot(data.pca$x[ , 1], data.pca$x[ , 3])
#text(data.pca$x[ , 1], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)

#plot(data.pca$x[ , 2], data.pca$x[ , 3])
#text(data.pca$x[ , 2], data.pca$x[ , 3], labels=rownames(data.pca$x), cex= 0.7, pos=3)
dev.off()



##### SVM #####

##### Predictions from leave-one-out SVM on the training data - all orths.

data <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[, ]))))
predictionAll <- data.frame(matrix(nrow = 1, ncol = 2 * num.species.train))
colnames(predictionAll) <- colnames(data)
set <- 1:(2 * num.species.train)

for(j in 1:num.species.train){
  
  idx.test <- c(2 * j - 1, 2 * j)
  idx.train <- set[set != idx.test]
  
  data.train <- data[, idx.train]
  data.test <- data[, idx.test]
  
  # Perform a grid search to optimise SVM parameters.
  tuneResult <- tune("svm", train.x = t(data.train), train.y = caste[idx.train], probability = TRUE, kernel = "radial", tunecontrol = tune.control(sampling = "cross", cross = num.species.train - 2))
  
  # Final classifier.
  classifier <- tuneResult$best.model
  
  # Make predictions for the test data.
  pred <- predict(classifier, t(data.test), type = "class", probability = TRUE)
  
  # Add to data frame.
  predictionAll[1, idx.test] <- pred
  
}


## Predictions from leave-one-out SVM on the training data - p-value significant orths.

data <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[p.coeffs.sig, ]))))
predictionP <- data.frame(matrix(nrow = 1, ncol = 2 * num.species.train))
colnames(predictionP) <- colnames(data)
set <- 1:(2 * num.species.train)

for(j in 1:num.species.train){
  
  idx.test <- c(2 * j - 1, 2 * j)
  idx.train <- set[set != idx.test]
  
  data.train <- data[, idx.train]
  data.test <- data[, idx.test]
  
  # Perform a grid search to optimise SVM parameters.
  tuneResult <- tune("svm", train.x = t(data.train), train.y = caste[idx.train], probability = TRUE, kernel = "radial", tunecontrol = tune.control(sampling = "cross", cross = num.species.train - 2))
  
  # Final classifier.
  classifier <- tuneResult$best.model
  
  # Make predictions for the test data.
  pred <- predict(classifier, t(data.test), type = "class", probability = TRUE)
  
  # Add to data frame.
  predictionP[1, idx.test] <- pred
  
}

# Predictions from leave-one-out SVM on the training data - q-value significant orths.

data <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[q.coeffs.sig, ]))))
predictionQ <- data.frame(matrix(nrow = 1, ncol = 2 * num.species.train))
colnames(predictionQ) <- colnames(data)
set <- 1:(2 * num.species.train)

for(j in 1:num.species.train){
  
  idx.test <- c(2 * j - 1, 2 * j)
  idx.train <- set[set != idx.test]
  
  data.train <- data[, idx.train]
  data.test <- data[, idx.test]
  
  # Perform a grid search to optimise SVM parameters.
  tuneResult <- tune("svm", train.x = t(data.train), train.y = caste[idx.train], probability = TRUE, kernel = "radial", tunecontrol = tune.control(sampling = "cross", cross = num.species.train - 2))
  
  # Final classifier.
  classifier <- tuneResult$best.model
  
  # Make predictions for the test data.
  pred <- predict(classifier, t(data.test), type = "class", probability = TRUE)
  
  # Add to data frame.
  predictionQ[1, idx.test] <- pred
  
}

# Create pdf for plot.
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Classification_leave-one-out-species-scaled.pdf"),
    width = 3 * 1.5 * 3.54,
    height = 1.5 * 3.54)

# Plot
par(las=2)
par(mfrow=c(1,3))
par(mar=c(8,8,1,1)) # adjust as needed
barplot(as.matrix(predictionAll[1,  ]), ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "All")
abline(h=0.5)
barplot(as.matrix(predictionP[1, ]), ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "P-val sig")
abline(h=0.5)
barplot(as.matrix(predictionQ[1, ]), ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "Q-val sig")
abline(h=0.5)

dev.off()







##### Predictions on the test data - all orths.

data.train <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[, ]))))
data.test <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[ ]))))

# Perform a grid search to optimise SVM parameters.
tuneResult <- tune("svm", train.x = t(data.train), train.y = caste, probability = TRUE, kernel = "linear", tunecontrol = tune.control(sampling = "cross", cross = num.species.train))

# Final classifier.
classifier <- tuneResult$best.model

# Make predictions for the test data.
predictionAll <- predict(classifier, t(data.test), type = "class", probability = TRUE)



## Predictions on the test data - p-value significant orths.

data.train <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[p.coeffs.sig, ]))))
data.test <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[p.coeffs.sig, ]))))

# Perform a grid search to optimise SVM parameters.
tuneResult <- tune("svm", train.x = t(data.train), train.y = caste, probability = TRUE, kernel = "linear", tunecontrol = tune.control(sampling = "cross", cross = num.species.train))

# Final classifier.
classifier <- tuneResult$best.model

# Make predictions for the test data.
predictionP <- predict(classifier, t(data.test), type = "class", probability = TRUE)



## Predictions on the test data - q value significant orths. 

data.train <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[q.coeffs.sig, ]))))
data.test <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.test[q.coeffs.sig, ]))))

# Perform a grid search to optimise SVM parameters.
tuneResult <- tune("svm", train.x = t(data.train), train.y = caste, probability = TRUE, kernel = "linear", tunecontrol = tune.control(sampling = "cross", cross = num.species.train))

# Final classifier.
classifier <- tuneResult$best.model

# Make predictions for the test data.
predictionQ <- predict(classifier, t(data.test), type = "class", probability = TRUE)


# Create pdf for plot.
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Classification_Pcan-species-scaled.pdf"),
    width = 3 * 1.5 * 3.54,
    height = 1.5 * 3.54)

par(las=2)
par(mfrow=c(1,3))
par(mar=c(8,8,1,1)) # adjust as needed

barplot(as.matrix(predictionAll)[,1], ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "All")
abline(h=0.5)
barplot(as.matrix(predictionP)[,1], ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "P-val sig")
abline(h=0.5)
barplot(as.matrix(predictionQ)[,1], ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "Q-val sig")
abline(h=0.5)

dev.off()


# Create pdf for plot.
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Classification_Pcan-species-scaled.wide.pdf"),
    width = 3 * 1.5 * 3.54 * 2,
    height = 1.5 * 3.54)

par(las=2)
par(mfrow=c(1,3))
par(mar=c(8,8,1,1)) # adjust as needed

barplot(as.matrix(predictionAll)[,1], ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "All")
abline(h=0.5)
barplot(as.matrix(predictionP)[,1], ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "P-val sig")
abline(h=0.5)
barplot(as.matrix(predictionQ)[,1], ylim = c(0, 1), ylab = "Likelihood of being Queen", main = "Q-val sig")
abline(h=0.5)

dev.off()


##### Predictions from leave-one-out ordering genes by beta-coefficient and progressively filtering.

data <- scale(t(scale(t(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train[, ]))))
num.genes <- dim(matrix.counts.log2.quantile.species.scaled.orth.quantile.scale.train)[1]
predictionAll <- data.frame(matrix(nrow = num.genes, ncol = 2 * num.species.train))
colnames(predictionAll) <- colnames(data)
set <- 1:(2 * num.species.train)

filter.pc <- c(seq(0.01, 1, 0.01))
filter.num <- floor(filter.pc * num.genes)

for(i in 1:length(filter.num)){

  top.n <- filter.num[i] 
  cutoff <- abs(b.coeffs)[order(abs(b.coeffs), decreasing = TRUE)][top.n]
  b.coeffs.sig <- abs(b.coeffs) >= cutoff  
  
  for(j in 1:num.species.train){
    
    idx.test <- c(2 * j - 1, 2 * j)
    idx.train <- set[set != idx.test]
    
    data.train <- data[b.coeffs.sig, idx.train]
    data.test <- data[b.coeffs.sig, idx.test]
    
    # Perform a grid search to optimise SVM parameters.
    tuneResult <- tune("svm", train.x = t(data.train), train.y = caste[idx.train], probability = TRUE, kernel = "radial", tunecontrol = tune.control(sampling = "cross", cross = num.species.train - 2))
    
    # Final classifier.
    classifier <- tuneResult$best.model
    
    # Make predictions for the test data.
    pred <- predict(classifier, t(data.test), type = "class", probability = TRUE)
    
    # Add to data frame.
    predictionAll[i, idx.test] <- pred
    
  }
  
  cat(paste("Processe step (out of 100): ", i, "\n"))
  
}

# Fix out-of-bounds errors.
predictionAll[predictionAll < 0] <- 0
predictionAll[predictionAll > 1] <- 1

species <- colnames(predictionAll)
colours <- c("red", "blue", "green", "orange", "purple","yellow","lightgreen","black", "grey")

# Create pdf for plot.
pdf(file = paste("FIGURES/Figure_of_Classifiers/VERSION_RUN/Classification_confidence_leave-one-out-filterting-species-scaled.pdf"),
    width = 1.5 * 3.54,
    height = 1.5 * 3.54)

plot(rev(predictionAll[1:100,1]), type = "l", ylim = c(0, 1), col=colours[1], xlab="Percentage kept after filtering", ylab="Classification confidence values", xaxt = "n", main = "Change in classification confidence as % of genes decreased - species-scaled")
lines(rev(predictionAll[1:100,3]), col=colours[2])
lines(rev(predictionAll[1:100,5]), col=colours[3])
lines(rev(predictionAll[1:100,7]), col=colours[4])
lines(rev(predictionAll[1:100,9]), col=colours[5])
#lines(rev(predictionAll[1:100,11]), col=colours[6])
#lines(rev(predictionAll[1:100,13]), col=colours[7])
#lines(rev(predictionAll[1:100,15]), col=colours[8])
#lines(rev(predictionAll[1:100,17]), col=colours[9])
legend(0,1,legend = species[c(1,3,5,7,9)], col = colours, lty=1)
axis(1, at=1:100,labels=c(100:1))

dev.off()
