
#D<-read.table("data.txt",header=TRUE )
install.packages("tidyverse")
install.packages("limma")

library("magrittr")
library("ggplot2")
library("ggrepel")
library("tidyverse")
library("dplyr")
library("limma")
library(vegan)

data <- read.table("D:\\مشروع_التخرج\\data\\Human__TCGA_LAML__BDGSC__miRNASeq__GA_miR__01_28_2016__BI__Gene__Firehose_RPM_log2.cct", header = TRUE, row.names = 1)
dim(data)
summary(data)
head(data[1:5,1:5])
data =t(data )
data
class(data)
head(data[1:5,1:5])

data=as.data.frame(data)
data <- type.convert(data, as.is = TRUE)
head(data[1:5,1:5])

data=data %>% mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))

head(data[1:5,1:5])

num_zeros <- rowSums(data == 0.0000)
print(ncol(data))
# Calculate the proportion of zeros in each row
prop_zeros <- num_zeros / ncol(data)
# Plot the distribution of proportion of zeros
hist(prop_zeros)

prop_zeros <- apply(data, 2, function(x) sum(x == 0.000) / length(x))

threshold <- 0.1
cols_to_remove <- which(prop_zeros > threshold) #there is som columns to remove
data=data[,-cols_to_remove]
dim(data)

plot(density(apply(data, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "after_remove_zeros", xlab = "", 
     ylab = "")

boxplot(data[, 1:10], names = colnames(data[, 1:10]), las = 2, col = "lightgreen", horizontal = T , main="After Normalization by remove zeros")

colSums(is.na(data))[1:5]
Log2Transformed <- log2((abs(data))^(sign(data)))
head(Log2Transformed[1:5,1:5])
density(apply(Log2Transformed, 2, mean, na.rm = TRUE))
plot(density(apply(Log2Transformed, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
boxplot(Log2Transformed[, 1:10], names = colnames(Log2Transformed[, 1:10]), las = 2, col = "lightgreen", horizontal = T)
scaled_df=scale(Log2Transformed,center = TRUE, scale = TRUE)
head(scaled_df[1:5,1:5])

plot(density(apply(scaled_df, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")

boxplot(scaled_df[, 1:10], names = colnames(scaled_df[, 1:10]), las = 2, col = "lightgreen", horizontal = T)
dim(data)
write.csv(scaled_df,"Normalized mirna_data.csv")

##################################################################################################
data <- read.table("D:\\مشروع_التخرج\\data\\Human__TCGA_LAML__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct", header = TRUE, row.names = 1)
dim(data)
summary(data)

head(data[1:5,1:5])
data =t(data )
class(data)
data=as.data.frame(data)
data <- type.convert(data, as.is = TRUE)
head(data[1:5,1:5])
data=data %>% mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))

num_zeros <- rowSums(data == 0.0000)

# Calculate the proportion of zeros in each row
prop_zeros <- num_zeros / ncol(data)

# Plot the distribution of proportion of zeros
hist(prop_zeros)
# Calculate proportion of 0.00 values in each column
prop_zeros <- apply(data, 2, function(x) sum(x == 0.000) / length(x))

# Set threshold for proportion of 0.00 values
threshold <- 0.1   #most common
# Identify columns with proportion of 0.00 values above threshold
cols_to_remove <- which(prop_zeros > threshold) #there is no columns to remove
dim(data)
plot(density(apply(data, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
boxplot(data[, 1:10], names = colnames(data[, 1:10]), las = 2, col = "lightgreen", horizontal = T)
Log2Transformed <- log2((abs(data))^(sign(data)))
head(Log2Transformed[1:5,1:5])
plot(density(apply(Log2Transformed, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
boxplot(Log2Transformed[, 1:10], names = colnames(Log2Transformed[, 1:10]), las = 2, col = "lightgreen", horizontal = T)
scaled_df=scale(Log2Transformed,center = TRUE, scale = TRUE)
plot(density(apply(scaled_df, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
boxplot(scaled_df[, 1:10], names = colnames(scaled_df[, 1:10]), las = 2, col = "lightgreen", horizontal = T)
colSums(is.na(scaled_df))[1:5]
write.csv(scaled_df,"Normalized rnaseq_data.csv")
data <- read.table("D:\\مشروع_التخرج\\data\\Human__TCGA_LAML__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header = TRUE, row.names = 1)
dim(data)
data =t(data )
data=as.data.frame(data)
data <- type.convert(data, as.is = TRUE)
data<- na.omit(data)
dim(data)
write.csv(scaled_df,"Normalized_clinical_data.csv")


data <- read.table("D:\\مشروع_التخرج\\data\\Human__TCGA_LAML__JHU_USC__Methylation__Meth450__01_28_2016__BI__Gene__Firehose_Methylation_Prepocessor.cct", header=TRUE, na.strings=c("NA","NaN", ""),row.names=1, sep="\t")
dim(data)
head(data[1:5,1:5])
data =t(data )
data=as.data.frame(data)
data <- type.convert(data, as.is = TRUE)
head(data[1:5,1:5])
data=data %>% mutate_if(is.numeric, ~replace_na(.,mean(., na.rm = TRUE)))
# Count the number of zeros in each row
num_zeros <- rowSums(data == 0.0000)

# Calculate the proportion of zeros in each row
prop_zeros <- num_zeros / ncol(data)

# Plot the distribution of proportion of zeros
hist(prop_zeros)
# Calculate proportion of 0.00 values in each column
prop_zeros <- apply(data, 2, function(x) sum(x == 0.000) / length(x))

# Set threshold for proportion of 0.00 values
threshold <- 0.1   #most common
# Identify columns with proportion of 0.00 values above threshold
cols_to_remove <- which(prop_zeros > threshold) #there is no columns to remove
dim(data)
plot(density(apply(data, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
boxplot(data[, 1:10], names = colnames(data[, 1:10]), las = 2, col = "lightgreen", horizontal = T)
colSums(is.na(data))[1:5]
Log2Transformed <- log2((abs(data))^(sign(data)))
head(Log2Transformed[1:5,1:5])
plot(density(apply(Log2Transformed, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
boxplot(Log2Transformed[, 1:10], names = colnames(Log2Transformed[, 1:10]), las = 2, col = "lightgreen", horizontal = T)
scaled_df=scale(Log2Transformed,center = TRUE, scale = TRUE)
head(scaled_df[1:5,1:5])
plot(density(apply(scaled_df, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
boxplot(scaled_df[, 1:10], names = colnames(scaled_df[, 1:10]), las = 2, col = "lightgreen", horizontal = T)
colSums(is.na(scaled_df))[1:5]
write.csv(scaled_df,"Normalized methylation_data.csv")

#####################################################################################################

