## Code extract position of samples based on river basin
## Extract random samples fromthose samples and input into SNAPP
#install.packages("elevatr")
library(sf)
library(ggplot2)
library(elevatr)
library(vcfR)
library(adegenet)
library(ape)
library(adegenet)
library(poppr)

# Read in world map
worldmap <- read_sf("data/world_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")

# read in level 3 and 4 of hydrosheds
hydrobasins_3 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev03_v1c.shp")
hydrobasins_4 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev04_v1c.shp")
hydrobasins_5 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev05_v1c.shp")

# Read in VCF
# Output file location
SNP.library.name <- "titia.mysnps"
species <- "titia"
dir.path <- paste0("data/SNPs/")
filter_para <- ".thin0_10000.missfilter"
# Plot output file location
plot.dir <- paste0("plots/SNPs")
dir.create(plot.dir)

vcf <- read.vcfR(paste0(dir.path,SNP.library.name,filter_para,".vcf.gz"))

# Get sample names
samples <- colnames(vcf@gt)
sample_map <- read.csv("data/Sequenced_samples140_v4.csv")

#Extracting sample information
sites <- data.frame(samples)
sites$samples
sites <- merge(sites, sample_map, by.x = "samples", by.y = "Novogene.id")
sites$Lat <- as.numeric(sites$Lat)
sites$Long <- as.numeric(sites$Long)

# Rename vcf with lab codes
sites$Unique.ID[which(sites$Unique.ID=="TM01")] <- c("TM01c1", "TM01c2")
colnames(vcf@gt) <- sites$Unique.ID[match(colnames(vcf@gt), sites$samples)]
colnames(vcf@gt)[1] <- "FORMAT"
colnames(vcf@gt)[duplicated(colnames(vcf@gt))]


# Convert to sf
dim(sites)
sites_sf <- st_as_sf(sites, coords = c("Long", "Lat"))
dim(sites_sf)
sf_use_s2(FALSE)
st_crs(sites_sf) <- st_crs(hydrobasins_4)

plot(sites_sf)

# Extract basin info for each sample point
samples_basin_3 <- st_intersection(hydrobasins_3, sites_sf)
samples_basin_4 <- st_intersection(hydrobasins_4, sites_sf)
samples_basin_5 <- st_intersection(hydrobasins_5, sites_sf)

sites_sf[!sites_sf$samples%in%samples_basin_3$samples,]

# ggplot() +
#  geom_sf(data = hydrobasins_5, aes(fill = as.factor(HYBAS_ID)), show.legend = F) +
#  geom_sf(data = samples_basin_5, aes(fill = as.factor(HYBAS_ID)), size = 3, shape = 21, show.legend = F)

#number of river basins assigned to samples
length(unique(samples_basin_3$HYBAS_ID))
length(unique(samples_basin_4$HYBAS_ID))
length(unique(samples_basin_5$HYBAS_ID))

# Get highest coverage sample from each region
vcf.bi <- vcf[is.biallelic(vcf),]
vcf.bi@fix[,1] <- gsub("\\.", "_", vcf.bi@fix[,1])

my_genind_ti <- vcfR2genind(vcf.bi, sep = "/", return.alleles = TRUE)

# Get proporation of missing SNPS for each sample
sample.miss <- propTyped(my_genind_ti, by = "ind")
hist(sample.miss)

# Combine hydrobasin and sample coverage
sites$samples==samples_basin_3$samples
sites$basin_3 <- samples_basin_3$HYBAS_ID
sites$basin_4 <- samples_basin_4$HYBAS_ID
sites$basin_5 <- samples_basin_5$HYBAS_ID
sites$covarage <- sample.miss

# Get the sample which the highest coverage for each river basin
max.coverage.basin_4 <- mapply(unique(sites$basin_4), FUN = function(x) sites$Unique.ID[sites$basin_4==x][which.max(sites$covarage[sites$basin_4==x])])
max.coverage.basin_5 <- mapply(unique(sites$basin_5), FUN = function(x) sites$Unique.ID[sites$basin_5==x][which.max(sites$covarage[sites$basin_5==x])])


# Get random N samples from each hydrobasins
run.num <- 5
N <- 1
rand.basin_5 <- list()
for(run in 1:run.num){
  rand.basin_5[[run]] <- mapply(unique(sites$basin_5), FUN = function(x) {
    if(length(sites$samples[sites$basin_5==x])>=N) {
      return(sample(sites$samples[sites$basin_5==x], size = N))
    }
    else {
      return(sample(sites$samples[sites$basin_5==x], size = length(sites$samples[sites$basin_5==x])))
    }
  }
  )
}
rand.basin_5
# max.coverage.basin_5 <- do.call("c", max.coverage.basin_5)


ggplot(sites[sites$Unique.ID%in%max.coverage.basin_5,]) +
  geom_sf(data = hydrobasins_5, aes(fill = as.factor(HYBAS_ID)), show.legend = F) +
  geom_point(aes(Long, Lat, size = covarage))

####################
## HYDROBASINS 4 ###
####################


#Subset geneid to the top coverage samples
my_genind_ti_SNPs.short <- my_genind_ti[max.coverage.basin_4,]
sites.SNPS.short <- sites[sites$Unique.ID%in%max.coverage.basin_4,]

#Remove non-varient sites
SNP.allele.num <- vector()
for(i in 1:dim(my_genind_ti_SNPs.short@tab)[2]){
  temp.gt <- unique(my_genind_ti_SNPs.short@tab[,i])
  temp.gt <- temp.gt[!is.na(temp.gt)]
  SNP.allele.num[i] <- length(temp.gt)
}
table(SNP.allele.num)

my_genind_ti_SNPs.short <- my_genind_ti_SNPs.short[,SNP.allele.num!=1]

# Get smaples from vcf file
vcf.bi.short <- vcf.bi[,c("FORMAT",sites.SNPS.short$Unique.ID)]
poly.snps <- substr(colnames(my_genind_ti_SNPs.short@tab), start = 1, stop =  nchar(colnames(my_genind_ti_SNPs.short@tab))-2)
all_snps <- paste(vcf.bi.short@fix[,1],vcf.bi.short@fix[,2], sep = "_")

vcf.bi.short <- vcf.bi.short[all_snps%in%poly.snps,]
vcf.bi.short <- vcf.bi.short[is.biallelic(vcf.bi.short),]

gt <- extract.gt(vcf.bi.short, return.alleles = TRUE, convertNA = TRUE)

gt[gt=="A/A"] <- "A"
gt[gt=="T/T"] <- "T"
gt[gt=="G/G"] <- "G"
gt[gt=="C/C"] <- "C"
gt[gt=="A/G"] <- "R"
gt[gt=="G/A"] <- "R"
gt[gt=="C/T"] <- "Y"
gt[gt=="T/C"] <- "Y"
gt[gt=="A/C"] <- "M"
gt[gt=="C/A"] <- "M"
gt[gt=="G/T"] <- "K"
gt[gt=="T/G"] <- "K"
gt[gt=="C/G"] <- "S"
gt[gt=="G/C"] <- "S"
gt[gt=="A/T"] <- "W"
gt[gt=="T/A"] <- "W"
gt[gt=="."] <- "?"

# Rename file with population name
colnames(gt)==sites.SNPS.short$Unique.ID
colnames(gt) <- paste(sites.SNPS.short$basin_4, sites.SNPS.short$samples, sep = "_")

# Create in put

N <- 1
SNAPP_run_name <- "max_cov"
#Write nexus
write.nexus.data(t(gt), file = paste0(dir.path, "SNAPP/", SNP.library.name,"-SNAPP-hydro_4_max_cov.nex"), 
                 format = "DNA", missing = "?", interleaved = F)

## Set up files following https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md
## (1) Phylip file
colnames(gt) <- sites.SNPS.short$samples

ape::write.dna(t(gt), file = paste0(dir.path, "SNAPP/", SNP.library.name,"-SNAPP-hydro_4_max_cov.phy"),
               format ="interleaved", nbcol = -1, colsep = "")
## (2) Species table
species.df <- data.frame(species = sites.SNPS.short$basin_4, sample = sites.SNPS.short$Unique.ID)
write.table(species.df, file = paste0(dir.path, "SNAPP/", SNP.library.name,"-SNAPP-hydro_4_max_cov.txt"),
            row.names = F,quote = F, sep = "\t")
# (3) Theta priors from Stranding et al 2022
# All Hetaerina samples
# normal(offset,mean,sigma)
# Titia divergence from G-Phocs run without migration
hist(rnorm(1000, 3.4, 0.1))
contrant.df <- data.frame(x = "normal(0,3.4,0.1)", y = "crown", 
                          z = paste0(unique(sites.SNPS.short$basin_4), sep = ",",collapse = ""))
# Americana/calverti
#remove trailing comma
contrant.df$z <- substr(contrant.df$z, start = 1 , stop = nchar(contrant.df$z)-1)

write.table(contrant.df, file = paste0(dir.path, "SNAPP/", SNP.library.name,"-SNAPP-hydro_4_max_cov.con.txt"),
            row.names = F, col.names = F, quote = F, sep = "\t")


####################
## HYDROBASINS 5 ###
####################


#Subset geneid to the top coverage samples
max.coverage.basin_5 <- mapply(unique(sites$basin_5), FUN = function(x) sites$Unique.ID[sites$basin_5==x][which.max(sites$covarage[sites$basin_5==x])])


my_genind_ti_SNPs.short <- my_genind_ti[max.coverage.basin_5,]
sites.SNPS.short <- sites[sites$Unique.ID%in%max.coverage.basin_5,]

#Remove non-varient sites
SNP.allele.num <- vector()
for(i in 1:dim(my_genind_ti_SNPs.short@tab)[2]){
  temp.gt <- unique(my_genind_ti_SNPs.short@tab[,i])
  temp.gt <- temp.gt[!is.na(temp.gt)]
  SNP.allele.num[i] <- length(temp.gt)
}
table(SNP.allele.num)

my_genind_ti_SNPs.short <- my_genind_ti_SNPs.short[,SNP.allele.num!=1]

# Get smaples from vcf file
vcf.bi.short <- vcf.bi[,c("FORMAT",sites.SNPS.short$Unique.ID)]
poly.snps <- substr(colnames(my_genind_ti_SNPs.short@tab), start = 1, stop =  nchar(colnames(my_genind_ti_SNPs.short@tab))-2)
all_snps <- paste(vcf.bi.short@fix[,1],vcf.bi.short@fix[,2], sep = "_")

vcf.bi.short <- vcf.bi.short[all_snps%in%poly.snps,]
vcf.bi.short <- vcf.bi.short[is.biallelic(vcf.bi.short),]

gt <- extract.gt(vcf.bi.short, return.alleles = TRUE, convertNA = TRUE)

gt[gt=="A/A"] <- "A"
gt[gt=="T/T"] <- "T"
gt[gt=="G/G"] <- "G"
gt[gt=="C/C"] <- "C"
gt[gt=="A/G"] <- "R"
gt[gt=="G/A"] <- "R"
gt[gt=="C/T"] <- "Y"
gt[gt=="T/C"] <- "Y"
gt[gt=="A/C"] <- "M"
gt[gt=="C/A"] <- "M"
gt[gt=="G/T"] <- "K"
gt[gt=="T/G"] <- "K"
gt[gt=="C/G"] <- "S"
gt[gt=="G/C"] <- "S"
gt[gt=="A/T"] <- "W"
gt[gt=="T/A"] <- "W"
gt[gt=="."] <- "?"

# Rename file with population name
colnames(gt)==sites.SNPS.short$Unique.ID
colnames(gt) <- paste(sites.SNPS.short$basin_5, sites.SNPS.short$samples, sep = "_")


# Stats on how many samples and SNPs there are
colnames(gt)
dim(gt)

ggplot(sites_sf[sites_sf$samples%in%max.coverage.basin_5,]) +
  geom_sf(data = hydrobasins_5, aes(fill = as.factor(HYBAS_ID)), show.legend = F) +
  geom_sf(size = 3)

# Create in put

N <- 1
SNAPP_run_name <- "max_cov"
#Write nexus
write.nexus.data(t(gt), file = paste0(dir.path, "SNAPP/", SNP.library.name,"-SNAPP-hydro_5_max_cov.nex"), 
                 format = "DNA", missing = "?", interleaved = F)

## Set up files following https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md
## (1) Phylip file
colnames(gt) <- sites.SNPS.short$Unique.ID

ape::write.dna(t(gt), file = paste0(dir.path, "SNAPP/", SNP.library.name,"-SNAPP-hydro_5_max_cov.phy"),
               format ="interleaved", nbcol = -1, colsep = "")
## (2) Species table
species.df <- data.frame(species = paste0(sites.SNPS.short$basin_5,"-",sites.SNPS.short$Unique.ID), sample = sites.SNPS.short$samples)
write.table(species.df, file = paste0(dir.path, "SNAPP/", SNP.library.name,"-SNAPP-hydro_5_max_cov.txt"),
            row.names = F,quote = F, sep = "\t")
# (3) Theta priors from Stranding et al 2022
# All Hetaerina samples
# normal(offset,mean,sigma)
# Titia divergence from SNAPP run without migration
hist(rnorm(1000, 3.4, 0.1), breaks = 30)
contrant.df <- data.frame(x = "normal(0,3.5,0.5)", y = "crown", 
                          z = paste0(sites.SNPS.short$basin_5,"-",sites.SNPS.short$Unique.ID, sep = ",",collapse = ""))
# Americana/calverti
#remove trailing comma
contrant.df$z <- substr(contrant.df$z, start = 1 , stop = nchar(contrant.df$z)-1)

write.table(contrant.df, file = paste0(dir.path, "SNAPP/", SNP.library.name,"-SNAPP-hydro_5_max_cov.con.txt"),
            row.names = F, col.names = F, quote = F, sep = "\t")

####################
## HYDROBASINS 5 - 1 random sample ###
####################
run <- 1
for(run in 1:run.num){
  #Subset geneid to the top coverage samples
  run.i.coverage.basin_5 <- rand.basin_5[[run]]
  
  my_genind_ti_SNPs.short <- my_genind_ti[run.i.coverage.basin_5,]
  sites.SNPS.short <- sites[sites$samples%in%run.i.coverage.basin_5,]
  
  #Remove non-varient sites
  SNP.allele.num <- vector()
  for(i in 1:dim(my_genind_ti_SNPs.short@tab)[2]){
    temp.gt <- unique(my_genind_ti_SNPs.short@tab[,i])
    temp.gt <- temp.gt[!is.na(temp.gt)]
    SNP.allele.num[i] <- length(temp.gt)
  }
  table(SNP.allele.num)
  
  my_genind_ti_SNPs.short <- my_genind_ti_SNPs.short[,SNP.allele.num!=1]
  
  # Get smaples from vcf file
  vcf.bi.short <- vcf.bi[,c("FORMAT",sites.SNPS.short$samples)]
  poly.snps <- substr(colnames(my_genind_ti_SNPs.short@tab), start = 1, stop =  nchar(colnames(my_genind_ti_SNPs.short@tab))-2)
  all_snps <- paste(vcf.bi.short@fix[,1],vcf.bi.short@fix[,2], sep = "_")
  
  vcf.bi.short <- vcf.bi.short[all_snps%in%poly.snps,]
  vcf.bi.short <- vcf.bi.short[is.biallelic(vcf.bi.short),]
  
  gt <- extract.gt(vcf.bi.short, return.alleles = TRUE, convertNA = TRUE)
  
  gt[gt=="A/A"] <- "A"
  gt[gt=="T/T"] <- "T"
  gt[gt=="G/G"] <- "G"
  gt[gt=="C/C"] <- "C"
  gt[gt=="A/G"] <- "R"
  gt[gt=="G/A"] <- "R"
  gt[gt=="C/T"] <- "Y"
  gt[gt=="T/C"] <- "Y"
  gt[gt=="A/C"] <- "M"
  gt[gt=="C/A"] <- "M"
  gt[gt=="G/T"] <- "K"
  gt[gt=="T/G"] <- "K"
  gt[gt=="C/G"] <- "S"
  gt[gt=="G/C"] <- "S"
  gt[gt=="A/T"] <- "W"
  gt[gt=="T/A"] <- "W"
  gt[gt=="."] <- "?"
  
  # Rename file with population name
  colnames(gt)==sites.SNPS.short$samples
  colnames(gt) <- paste(sites.SNPS.short$basin_5, sites.SNPS.short$samples, sep = "_")
  
  
  # Stats on how many samples and SNPs there are
  colnames(gt)
  dim(gt)
  
  ggplot(sites_sf[sites_sf$samples%in%run.i.coverage.basin_5,]) +
    geom_sf(data = hydrobasins_5, aes(fill = as.factor(HYBAS_ID)), show.legend = F) +
    geom_sf(size = 3)
  
  # Create in put
  
  
  SNAPP_run_name <- paste0("-SNAPP-hydro_5_run", run)
  #Write nexus
  write.nexus.data(t(gt), file = paste0(dir.path, "SNAPP/", SNP.library.name, SNAPP_run_name,".nex"), 
                   format = "DNA", missing = "?", interleaved = F)
  
  ## Set up files following https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md
  ## (1) Phylip file
  colnames(gt) <- sites.SNPS.short$samples
  
  ape::write.dna(t(gt), file = paste0(dir.path, "SNAPP/", SNP.library.name,SNAPP_run_name,".phy"),
                 format ="interleaved", nbcol = -1, colsep = "")
  ## (2) Species table
  species.df <- data.frame(species = paste0(sites.SNPS.short$basin_5,"-",sites.SNPS.short$samples), sample = sites.SNPS.short$samples)
  write.table(species.df, file = paste0(dir.path, "SNAPP/", SNP.library.name,SNAPP_run_name,".txt"),
              row.names = F,quote = F, sep = "\t")
  # (3) Theta priors from Stranding et al 2022
  # All Hetaerina samples
  # normal(offset,mean,sigma)
  # Titia divergence from SNAPP run without migration
  ## hist(rnorm(10000, 3.4, 0.1), breaks = 30)
  contrant.df <- data.frame(x = "normal(0,3.5,0.5)", y = "crown", 
                            z = paste0(sites.SNPS.short$basin_5,"-",sites.SNPS.short$samples, sep = ",",collapse = ""))
  # Americana/calverti
  #remove trailing comma
  contrant.df$z <- substr(contrant.df$z, start = 1 , stop = nchar(contrant.df$z)-1)
  
  write.table(contrant.df, file = paste0(dir.path, "SNAPP/", SNP.library.name,SNAPP_run_name,".con.txt"),
              row.names = F, col.names = F, quote = F, sep = "\t")
  
}

