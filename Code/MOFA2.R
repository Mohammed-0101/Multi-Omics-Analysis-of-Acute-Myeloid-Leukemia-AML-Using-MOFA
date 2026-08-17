if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

install.packages("dplyr")

library("limma")
library("vegan")
library("cluster")
library("factoextra")
library("gridExtra")
library("PerformanceAnalytics")
library("corrplot")
library("Hmisc")
library("RColorBrewer")
library("impute")
library("pathview")
library("glmnet")
library("reshape2")
library("readr")
library("magrittr")
library("ggrepel")
library("tidyverse")
library("ggfortify")
library("ggplot2")
library("ggrepel")
library("EnhancedVolcano")
library("TCGAbiolinks")
library("MOFA2")
library("MOFAdata")
library("data.table")
library("magrittr") 
library("dplyr") 
library("GGally")
library("reticulate")
reticulate::py_config()
mirna <- read.csv("D:\\مشروع_التخرج\\data\\Normalized mirna_data.csv")
row.names(mirna) = mirna$X
mirna <- mirna[-1]
mirna <- as.data.frame(t(mirna))

rna <- read.csv("D:\\مشروع_التخرج\\data\\Normalized rnaseq_data.csv")
row.names(rna) = rna$X
rna <- rna[-1]
rna <- as.data.frame(t(rna))


methyl <- read.csv("D:\\مشروع_التخرج\\data\\Normalized methylation_data.csv")
row.names(methyl) = methyl$X
methyl <- methyl[-1]
methyl <- as.data.frame(t(methyl))

#bring all samples names(columns) common between all omics types
common_cols=intersect(colnames(mirna),colnames(rna))
common_cols=intersect(common_cols,colnames(methyl))

mirna_upd=mirna[,common_cols]
rna_upd=rna[,common_cols]
methyl_upd=methyl[,common_cols]

mirna_mat=as.matrix(mirna_upd)
rna_mat=as.matrix(rna_upd)
methyl_mat=as.matrix(methyl_upd)




clinical_data=read.table("D:\\مشروع_التخرج\\data\\Human__TCGA_LAML__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header=TRUE, na.strings=c("NA","NaN", ""),row.names=1, sep="\t")
class(clinical_data)
#Filter it on common samples (columns)
dim(clinical_data)
clinical_data=t(clinical_data)
clinical_data=data.frame(clinical_data)
clinical_data$race <- replace(clinical_data$race, is.na(clinical_data$race), "white")
clinical_data <- subset(clinical_data, select = -overallsurvival)
#write.csv(clinical_data,"fff.csv")

clinical_data1=clinical_data[,common_cols]
clinical_data1 =t(clinical_data1 )
clinical_data1=as.data.frame(clinical_data1)
clinical_data1 <- type.convert(clinical_data1, as.is = TRUE)
clinical_data1 =t(clinical_data1 )

dim(clinical_data1)

data=list(mirna_mat,rna_mat,methyl_mat)
lapply(data, dim) 
lapply(data, class) 
MOFAobject<-create_mofa(data)

views_names(MOFAobject)=c('mirna','rna_seq','methylation')
MOFAobject
plot_data_overview(MOFAobject)


data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
head(data_opts,1)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 7
head(model_opts,1)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts,1)

Nsamples = sum(MOFAobject@dimensions$N)
clinical_data1_t=data.frame(t(clinical_data1))
clinical_data1_t$sample=rownames(clinical_data1_t)
# samples_metadata(MOFAobject) <- clinical_data1_t
head(MOFAobject@samples_metadata, n=3)
#------------------------------------------------------------------------
if (!requireNamespace("reticulate", quietly = TRUE)) {
    install.packages("reticulate")
  }

# Load the reticulate package

library("devtools")


# Specify the name of the Conda environment
conda_env_name <- "mofa_env2"

# Create a new Conda environment with Python 3.7
conda_create( envname = conda_env_name)

# Activate the Conda environment
use_condaenv(conda_env_name, required = TRUE)

# Install mofapy2 using pip
devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))

py_install("mofapy2")

# # Deactivate the Conda environment
# use_condaenv()

# Use the created environment in run_mofa
MOFAobject.trained <- run_mofa(MOFAobject, use_basilisk = T)  
  
#--------------------------------  
  


MOFAobject <- prepare_mofa(object = MOFAobject ,
                           data_options = data_opts , 
                           model_options = model_opts , 
                           training_options = train_opts)


outfile = paste0("D:\\model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject,outfile, use_basilisk=FALSE)
model <- load_model(outfile)
plot_data_overview(model)

head(model@samples_metadata, n=3)

head(model@cache$variance_explained$r2_total[[1]]) 
head(model@cache$variance_explained$r2_per_factor[[1]]) 
plot_variance_explained(model, x="view", y="factor")
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
model <- load_model(outfile)
Nsamples = sum(model@dimensions$N)
#We set sample metadata to the clinical data
sample_metadata=clinical_data1_t
samples_metadata(model) <- sample_metadata
column_to_rownames(samples_metadata)
head(model@samples_metadata, n=3)
head(model@cache$variance_explained$r2_total[[1]])
print(model)

# gender ethnicity  status   race  

#mirna     rna_seq     methylation


plot_factor(model, 
            factor = 1:3,
            color_by  = "race")
plot_factors(model, 
            factor = 1:3,
            color_by  = "race")

plot_factor(model, 
                 factors = c(1,2,3),
                 color_by  = "race",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)
# The output of plot_factor is a ggplot2 object that we can edit
p <- p + 
  scale_color_manual(values=c(t1="black", t2="red" , t3="green", t4="blue")) +
  scale_fill_manual(values=c(t1="black", t2="red",t3="green",t4="blue"))

print(p)

# gender ethnicity  status   race  

plot_data_heatmap(model,
                  view = "mirna",         # view of interest
                  factor = 1,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  cluster_rows = TRUE, cluster_cols = TRUE,
show_rownames = TRUE, show_colnames = TRUE,annotation_samples = "gender")

plot_data_scatter(model,
                  view = "mirna",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "gender"
)





plot_weights(model,
             view = "mirna",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(model,
             view = "rna_seq",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)
plot_weights(model,
             view = "methylation",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

###############################################
plot_top_weights(model,
                 view = "mirna",
                 factor = 1,
                 nfeatures = 10
)
