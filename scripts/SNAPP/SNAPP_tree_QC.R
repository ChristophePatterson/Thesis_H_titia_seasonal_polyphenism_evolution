library(ggtree)
library(tidyverse)
library(ape)
library(patchwork)

# Read in args
## Get command line argument
args <- commandArgs(trailingOnly = TRUE)
run.name <- args[1]

plot.dir <- "plots/SNAPP"
dir.create(plot.dir)

# Print name
print(run.name)

#  Get trace file
trace.df <- read_table(paste0(run.name, ".log"))
trace.df

# Get trees
# trees <- read_lines(paste0(run.name, ".trees"))
# writeLines(c(trees, "End;"), paste0(run.name, ".trees.tmp"))
trees <- ape::read.nexus(paste0(run.name, ".trees.tmp"))

# Set burn in
burn.in <- 0.2
min.tree <- round(length(trees)*burn.in)
max.tree <- length(trees)

#Extracting sample information
# Read in sample data
sample_map <- read_csv("data/All samples held in Durham_v17.csv") %>%
    select(Unique.ID, species,Year, Site.ID, Lat, Long, Country, Province, Ocean.drainage)
sites <- tibble(data.frame(samples = str_split_i(trees$STATE_0$tip.label, "-", 2)))
sites
# Remove samples not included
sample_map <- sample_map[sample_map$Unique.ID%in%sites$samples,]
names(sample_map)
sites <- merge(sites, sample_map, by.x = "samples", by.y ="Unique.ID")
sites$Lat <- as.numeric(sites$Lat)
sites$Long <- as.numeric(sites$Long)

# Join tree and sample
time.trees <- lapply(1:length(trees), function(i) {
  tree <- trees[[i]]
  tree <- fortify(tree) %>% mutate(tree=factor(i, levels=as.numeric(1:length(trees))),
                           samples = str_split_i(label, "-", 2))
  left_join(tree, sites)

})

time.trees[[1]]

cbPalette <- c("#E69F00", "#009E73","#D55E00","#0072B2","#999999", "#F0E442", "#56B4E9", "#CC79A7", "black")
# Plot multiple trees with aligned tips from muliple time points ## aes(colour=as.numeric(tree))
p.tree <- ggdensitree(time.trees[min.tree:max.tree], colour='black', alpha = .3) + geom_tiplab() +
  geom_tippoint(aes(color = Country, fill=Ocean.drainage), shape = 21, size = 2, stroke = 1.2) +
  scale_x_continuous(expand = c(0.1,1)) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = rev(cbPalette)) +
  theme_bw()

# Create trace plot without burn-in
p.trace.tmp <- ggplot(trace.df[min.tree:max.tree,]) +
    geom_line(aes(Sample, treeHeightLogger))

p.trace.inset <- ggplot(trace.df) +
    geom_line(aes(Sample, treeHeightLogger))


# Add in inset for all trees
p.trace <- p.trace.tmp # + inset_element(p.trace.inset, left = 0.6, right = 0.9, top = 0.6, bottom = 0.1)

## p.trace <- ggplot(trace.df[min.tree:max.tree,]) +
##     geom_line(aes(Sample, treeHeightLogger))

ggsave(paste0(plot.dir, "/", run.name,"_SNAPP.png"), p.tree / p.trace, width = 12, height = 8)
ggsave(paste0(plot.dir, "/", run.name,"_SNAPP.pdf"), p.tree / p.trace, width = 12, height = 8)