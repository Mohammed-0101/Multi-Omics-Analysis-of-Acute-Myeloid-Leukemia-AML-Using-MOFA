if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",version =3.18)
BiocManager::install("cyPlot", force = TRUE)
BiocManager::install("RCytoscape")
library(BioNERO)



path = "D:\\مشروع_التخرج\\data\\row_data"


data<- dfs2one(path, pattern =".cct$")

exp_filt <- replace_na(data)
sum(is.na(exp_filt ))
exp_filt<-exp_filt[0:188]

exp_filt <- remove_nonexp(exp_filt, method = "median", min_exp = 1)
dim(exp_filt)
exp_filt <- ZKfiltering(exp_filt)
dim(exp_filt)
exp_filt <- PC_correction(exp_filt)
final_exp<-exp_filt
#final_exp<-t(final_exp)
p <- plot_heatmap(final_exp, type = "samplecor", show_rownames = FALSE)
p



sft <- SFT_fit(final_exp, net_type = "signed hybrid")
sft
power <- sft$power
power
net <- exp2gcn(
  final_exp, SFTpower = power
)
names(net)
plot_dendro_and_colors(net)
plot_eigengene_network(net)
plot_ngenes_per_module(net)

"""plot_expression_profile(
  exp = final_exp, 
  metadata = final_exp[0:300,(2:4)],
  net = net, 
  plot_module = TRUE, 
  modulename = "yellow"
)"""
hubs <- get_hubs_gcn(final_exp, net)
head(hubs)


edges <- get_edge_list(net)
head(edges)


# Remove edges based on optimal scale-free topology fit
edges_filtered <- get_edge_list(net, filter = TRUE)
## The correlation threshold that best fits the scale-free topology is 0.7
dim(edges_filtered)
## [1] 588   3

# Remove edges based on p-value
edges_filtered <- get_edge_list(
  net, 
  filter = TRUE, method = "pvalue", 
  nSamples = ncol(final_exp)
)
dim(edges_filtered)
## [1] 921   3

# Remove edges based on minimum correlation
edges_filtered <- get_edge_list(
  net, 
  filter = TRUE, method = "min_cor", rcutoff = 0.7
)
dim(edges_filtered)
## [1] 588   3



ii<-plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net, 
  color_by = "module", 
  hubs = hubs
)

class(ii)
//"gg"     "ggplot"

igg<-plot_gcn(
  edgelist_gcn = edges_filtered, 
  net = net,
  color_by = "module",
  hubs = hubs,
  interactive = TRUE,
  dim_interactive = c(500, 500)
)





library(htmlwidgets)

saveWidget(igg, "network_plot.html", selfcontained = TRUE)



interpro_enrichment <- module_enrichment(
  net = net, 
  background_genes = rownames(final_exp)
)
interpro_enrichment[, -6]

####################################################################

library(RCy3)
library(igraph)

cytoscapePing()



createNetworkFromDataFrames(nodes =  nn ,edges =  e ,title="my first network", collection="DataFrame Example")

Edges<-edges_filtered
class(edges_filtered[1])

colnames(hubs)[1]="id"
colnames(hubs)[2]="group"
colnames(hubs)[3]="score"
h=hubs[1]
colnames(Edges)[1]="source"
colnames(Edges)[2]="target"

library(dplyr)

# Filter rows where target value exists in hubs$id and drop duplicates
nodesssss <- Edges %>%
  filter(Edges$target %in% Edges$source & Edges$source %in% Edges$target) %>%
  distinct()

n=nodesssss[!duplicated(nodesssss$source),1]
nn=data.frame(id=n)
nn

target<-Edges$target[c(Edges$target %in% hubs$id )]
cc=source[c(source%in%target)]
edddd<-Edges[c(Edges$source %in% nn$id ),]
e=edddd[!duplicated(edddd$source),]
write.csv(e,"edges_data.csv")
write.csv(nn,"nodes_data.csv")
nn=rbind(nn,c("hsa-mir-532"))
,"hsa-mir-493","hsa-mir-503","hsa-mir-532"
x=cc[1:47]
x=data.frame(source=x)
target=data.frame(target=target)
Edges<=data.frame(target,x)
eee=merge(target,x)
head(E)
head(Edges)
E=Edges[,1:2]
Edges1<-Edges[1:10,]
Edges1=cbind(Edges1,interaction=c("inhibits","interacts","activates","interacts","inhibits","interacts","inhibits","interacts","activates","interacts"))
Edges1


# Load necessary library
library(dplyr)

# Filter rows where target value exists in hubs$id and drop duplicates
edgesss <- Edges %>%
  filter(Edges$target %in% hubs$id & Edges$source %in% hubs$id) %>%
  distinct()

ee=edgesss[0:8,]

ee=cbind(ee,interaction=c("inhibits","interacts","activates","interacts","inhibits","interacts","inhibits","interacts"))




install.packages("RCy3")
install.packages("htmlwidgets")
install.packages("h2o")
library(RCy3)
library(htmlwidgets)
library(networkD3)
library(h2o)
cyPlot(igg)

####

install.packages("cyPlot")
url="C:/Users/DELL/Documents/network_plot.html"

loadNetworkFromURL(url, format = "cyjs", collection = "My Network")


nodes <- data.frame(id=c("node 0","node 1","node 2","node 3"),
                    group=c("A","A","B","B"), # categorical strings
                    score=as.double(c(20.26466,10.65468,15.5465,5.5468)), # integers
                    stringsAsFactors=FALSE)
edges <- data.frame(source=c("node 0","node 0","node 0","node 2"),
                    target=c("node 1","node 2","node 3","node 3"),
                    interaction=c("inhibits","interacts","activates","interacts"),  # optional
                    weight=c(5.1,3.0,5.2,9.9), # numeric
                    stringsAsFactors=FALSE)

createNetworkFromDataFrames(nodes,edges, title="my first network", collection="DataFrame Example")
colnames(nodes)[2]="dwdc"
n=nodes[1]
colnames(edges)[3]="Gene"
e=edges[,1:2]
str(e)
rownames(nodes)=NULL
rownames(Edges)=NULL

####

nodes <- data.frame(id=c("hsa-mir-134","hsa-mir-379","hsa-mir-127","hsa-mir-654","hsa-mir-655"),
                    group=c("green","green","blue","blue","blue") # categorical strings
                    )


edges <- data.frame(source=c("hsa-mir-134","hsa-mir-134","hsa-mir-134","hsa-mir-127"),
                    target=c("hsa-mir-379","hsa-mir-127","hsa-mir-654","hsa-mir-654"),
                    weight=c(.154,.540,.256,.9546))

createNetworkFromDataFrames(nodes =  nodes ,edges =  edges ,title="my first network", collection="DataFrame Example")
