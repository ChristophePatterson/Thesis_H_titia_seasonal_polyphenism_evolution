## Data extraction

#install.packages("mvMORPH")
library(mvMORPH)
library(ggplot2)
library(patchwork)
library(ggtree)
library(lubridate)
library(tidyverse)

# Raw phenotype data
raw.pheno <- read.csv("data/MASTER_titia_phenotypes_v1.0.csv")
peak.trough <- read.csv("data/peak_trough_v2.0_logit.csv")
peak.trough.season <- read.csv("data/peak_trough_season_dates.csv")

# read in tree
#tree <- treeio::read.mega(file = "data/SNAPP/Hetaerina_titia_ddRAD_titia_dg-SNAPP-hydro_5_max_cov.trees_old.Anon")
# Get in ape phylo format
#tree.phylo <- tree@phylo
# Remove tips that arn't retained in cluster
#tree.phylo <- drop.tip(tree.phylo, tip = which(!tree.phylo$tip.label%in%peak.trough$treename))

tree.phylo <-ape::read.tree(file="data/titia_tree_trimmed.tree")

# Rename tips based on cluster names
tree.phylo$tip.label <- peak.trough$CLUSTER[match(tree.phylo$tip.label, peak.trough$treename)]

# Convert to ggtree
tree.plot <- ggtree(tree.phylo, size = 1.2) + geom_tiplab(size = 6)
# Get order in which trees are plotted
plot.tip.order <- tree.plot$data$label[order(tree.plot$data$y[tree.plot$data$isTip], decreasing = T)]

peak.trough$CLUSTER <- factor(peak.trough$CLUSTER, levels = plot.tip.order)
raw.pheno$CLUSTER <- factor(raw.pheno$CLUSTER, levels = plot.tip.order)

# Get wether each raw phenotype measurement was collected in the peak or trough
pheno_season <- apply(raw.pheno ,MARGIN = 1, function(x){
  day <- as.numeric(x["Julian.date"])
  cluster <- x["CLUSTER"]
  if(is.na(cluster)){return(NA)}
  cluster_season <- peak.trough.season[which(peak.trough.season$cluster==cluster),]
  if(day>=cluster_season$peak.start&day<=cluster_season$peak.end){return("peak")}
  if(day<91|day>287){return("trough")}
  if(!(day<91|day>287)&!(day>=cluster_season$peak.start&day<=cluster_season$peak.end)){return("off_season")}
  
}
)
# add to dataset
raw.pheno$cluster_season <- pheno_season

names(raw.pheno)

# Basic phenotype data
raw.pheno$Date.R <- as.Date(raw.pheno$DATE, format = "%d/%m/%Y")

range(raw.pheno[!is.na(raw.pheno$Standard_image),]$Date.R)
raw.pheno$Standard_image

table(raw.pheno$best.estimate.source)
sum(table(raw.pheno$best.estimate.source))
raw.pheno %>%
  group_by(by = best.estimate.source) %>%
  summarise(min.date = min(Date.R), max.date = max(Date.R), num.records = length(Date.R))

## What were the total number of phenotype measurements from each source
sum(!is.na(raw.pheno$iNat.AI.proportion))
sum(!is.na(raw.pheno$iNat.manual.proportion))
sum(!is.na(raw.pheno$Standard.manual.hw.proportion))
sum(!is.na(raw.pheno$Standard.AI.hw.proportion))

pheno.summary <- c(table(paste(!is.na(raw.pheno$iNat.AI.proportion),!is.na(raw.pheno$iNat.manual.proportion),
            !is.na(raw.pheno$Standard.manual.hw.proportion), !is.na(raw.pheno$Standard.AI.hw.proportion), 
            sep = "-")),c("TRUE-TRUE-TRUE-TRUE" = 0))
  
pheno_comb <-  apply(combn(c(T,F,T,F,T,F,T,F), 4), MARGIN = 2, function(x) paste(x, collapse = "-"))
pheno_comb <- pheno_comb[!duplicated(pheno_comb)]

combn(letters[1:4], 2)

### How many data points got assigned to each tip
raw.pheno$CLUSTER

tip.measurement.type <- raw.pheno %>%
  group_by(CLUSTER, best.estimate.source) %>%
  summarise(num.records = length(Date.R))

tip.labs <- stringr::str_split_i(tree@phylo$tip.label, "-", 2)

samples <- read.csv("data/All samples held in Durham_v17.csv", header = T)
samples <- samples[match(tip.labs, samples$Unique.ID),c("Unique.ID", "Long","Lat")]




