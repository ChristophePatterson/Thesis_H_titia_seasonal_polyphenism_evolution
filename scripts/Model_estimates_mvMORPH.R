## Modelling BM vs OU

#### NOTE
#  Code for simulating NULL datasets has been hashed out  using ###
#  Simulated datasets should be read in using readRDS

#### DATA ####
# Read in data and load packages
source("scripts/Data_extraction.R")
library(patchwork)
# Convert into format of mvMORPH
# Trait values with row names corrisponding to tips
var.mean <- cbind.data.frame(peak = peak.trough$mean[peak.trough$time=="peak"&!is.na(peak.trough$CLUSTER)], 
                             trough = peak.trough$mean[peak.trough$time=="trough"&!is.na(peak.trough$CLUSTER)])
row.names(var.mean) <- peak.trough$CLUSTER[peak.trough$time=="peak"&!is.na(peak.trough$CLUSTER)]

# Errors on trait values
var.se <- cbind.data.frame(peak = peak.trough$se[peak.trough$time=="peak"&!is.na(peak.trough$CLUSTER)], 
                           trough = peak.trough$se[peak.trough$time=="trough"&!is.na(peak.trough$CLUSTER)])
row.names(var.se) <- peak.trough$CLUSTER[peak.trough$time=="peak"&!is.na(peak.trough$CLUSTER)]
# Convert standar error into squared standard errors
var.sq.se <- var.se^2


# plot tree for determine different regimes
plot(tree.phylo)
nodelabels()
tiplabels()

#### TREES ####
# create tree with different regemes for Pacfic and Atlantic
# Separtate Pacific and Atlantic into different regemes
tree.phylo.r1 <- paintSubTree(tree.phylo, node = 9, state = "All")
tree.phylo.r2.P <- paintSubTree(tree.phylo, node=10, state="Pacific", anc.state="group_1",stem=TRUE)
# Separate most northern lineage 
tree.phylo.r2.HL <- paintSubTree(tree.phylo, node=6, state="High_Lat", anc.state="group_1",stem=TRUE)
# Separate both Pacific and most northern
tree.phylo.r3.HL.P <- paintSubTree(tree.phylo.r2.HL, node=10, state="Pacific", anc.state="group_1",stem=TRUE)

cols <- setNames(c("black", "orange", "skyblue"),colnames(tree.phylo.r3.HL.P$mapped.edge))
png("plots/Regime_trees.png", width = 2000, height = 1000)
par(mfrow = (c(1,4)))
plot(tree.phylo.r1, mar=c(0.2,0.1,1.1,0.2), xlim = c(0,7.5), lwd = 5, fsize = 8) ; text(1,8,"(a)", cex = 6)
plot(tree.phylo.r2.P, mar=c(0.2,0.1,1.1,0.2), xlim = c(0,7.5),lwd = 5, fsize = 8, colors = cols) ; text(1,8,"(b)", cex = 6)
plot(tree.phylo.r2.HL, mar=c(0.2,0.1,1.1,0.2), xlim = c(4.8,7),lwd = 5, fsize = 8, colors = cols) ; text(5,8,"(c)", cex = 6)
plot(tree.phylo.r3.HL.P, mar=c(0.2,0.1,1.1,0.2), xlim = c(4.8,7),lwd = 5, fsize = 8, colors = cols) ; text(5,8,"(d)", cex = 6)
par(mfrow= (c(1,1)))
dev.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#### MODELS ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

##### Co-evolution with BM ####
# Off-peak and peak melansistion are not correlated
## Brownian model with co-evolution
modelBM.ce <- mvBM(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r1)
## Brownian model without co-evolution
# NOTE CHECK WITH JONATHAN THAT THIS IS CORRECT
no.ce.matrix <- matrix(c(1,NA,NA,1), nrow = 2)
modelBM.no.ce <- mvBM(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r1, 
                      param = list(constraint = no.ce.matrix))
# Model comparison using corrected Akaike's Information Criterion
aicw(list(modelBM.ce, modelBM.no.ce), aicc=TRUE)
# Model comparison using non-parametric Likelihood Ratio Test
LRT(modelBM.ce, modelBM.no.ce)

## Hypothesis 2
##### Selective Regimes ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Test for different selective regimes
# Model brownian motion
modelBM1 <- mvBM(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r1)
# Model OU using varying different regimes
modelOUr1 <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo)
modelOUr2.P <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r2.P)
modelOUr2.HL <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r2.HL)
modelOUr3.HL.P <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r3.HL.P)
# Models without co-evolve
modelBM1.no.ce <- mvBM(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r1, param = list(constraint = no.ce.matrix))
modelOUr1.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo, param = list(sigma="constraint", alpha="constraint"))
modelOUr2.P.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r2.P, param = list(sigma="constraint", alpha="constraint"))
modelOUr2.HL.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r2.HL, param = list(sigma="constraint", alpha="constraint"))
modelOUr3.HL.P.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r3.HL.P, param = list(sigma="constraint", alpha="constraint"))

modelOUr3.HL.P.no.ce.RR <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r3.HL.P,
                                param = list(sigma="constraint", alpha="constraint",vcv="randomRoot"))


# Caculate AIC and compare
results <- list(modelBM1, modelOUr1, modelOUr2.P, modelOUr2.HL, modelOUr3.HL.P,
                modelBM1.no.ce, modelOUr1.no.ce, modelOUr2.P.no.ce, modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce)
weights <- aicw(results, aicc = TRUE)
weights

# Which has lowest AIC
all.AIC.df <- do.call("cbind.data.frame",weights[1:5])
all.AIC.df$regimes <- c("r1","r1","r2.Pac","r2.HL", "r3.HL.Pac")
all.AIC.df$co.evolve <- c(T,T,T,T,T,F,F,F,F,F)
all.AIC.df[order(all.AIC.df$AIC, decreasing = F),]

write.table(all.AIC.df[order(all.AIC.df$AIC, decreasing = F),], file = "data/models/AICc_model_regime_comparison.txt")

# Model OU with and without co-evolution
modelOUr1.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r1)
modelOUr1.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r1,
                        param = list(sigma="constraint", alpha="constraint"))


# table results
all.results <- list(modelBM.ce, modelBM.no.ce, modelOUr1.ce, modelOUr1.no.ce)
all.AIC <- aicw(all.results, aicc=T)
all.AIC 
all.AIC.df <- do.call("cbind.data.frame",all.AIC[1:5])
all.AIC.df$co_evove <- c(T,F,T,F)
all.AIC.df$regimes <- "r1"
write.table(all.AIC.df[order(all.AIC.df$AIC, decreasing = F),], file = "data/models/AICc_model_co-evolve_comparison.txt")

best.model <- all.results[[which.min(all.AIC$AIC)]]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#### LRT model comparison ####
# Model comparison using non-parametric LRT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Takes -2*(LL1-LL2) where LL1 & LL2 are the likelihood of two models you want to compare. 
# Then simulate data under the simplest model of the two
# fit both models on the simulated data, and compute the statistic from the model fit on all these simulated data. 
# Compare the statistic obtained on the observed to the null obtained by simulations.
## Number of simulations 
nsims <- 1000
asin.limits <- c(asin(sqrt(c(0,1))))

### Brownian motion simualtion (not valid)
modelBMr1.no.ce.sims <- mvSIM(tree = tree.phylo.r1, nsim = nsims, model = c("BM1"), param = modelBM.no.ce)


# Number of simulations that exceed max extent of trait range
sim.exceed.limits <- lapply(modelBMr1.no.ce.sims, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
unlist((table(unlist(sim.exceed.limits))/nsims)*100)
sim.exceed.limits.ntips <- lapply(modelBMr1.no.ce.sims, function(x) sum(x<asin.limits[1]|x>asin.limits[2]))
plot((table(unlist(sim.exceed.limits.ntips))/nsims)*100)

plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = c(-1,2), ylim=asin.limits)
lapply(modelBMr1.no.ce.sims, function(x) points(x, col = as.factor(row.names(x))))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

##### Co-evolution ####
# Simluate data using mvSIM to test the power of our data to compare between co-evolution
# for LRT the similar model should be simulated
modelOUr1.no.ce.sims <- mvSIM(tree = tree.phylo, nsim = nsims, model = c("OU1"), param = list(theta = modelOUr1.no.ce$theta,
                                                                                              alpha = modelOUr1.no.ce$alpha,
                                                                                              sigma = modelOUr1.no.ce$sigma))
# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = c(-1,2), ylim=asin.limits)
lapply(modelOUr1.no.ce.sims, function(x) points(x, col = as.factor(row.names(x))))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# Does simulations include values outside of valid phenotype arc sin of the square root of 0 to 1)
sim.exceed.limits <- lapply(modelOUr1.no.ce.sims, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
(table(unlist(sim.exceed.limits))/nsims)*100
sim.exceed.limits.ntips <- lapply(modelOUr1.no.ce.sims, function(x) sum(x<asin.limits[1]|x>asin.limits[2]))
(table(unlist(sim.exceed.limits.ntips))/nsims)*100

# Test for co-evolution
i <- 1
cov2cor(stationary(modelOUr1))

### modelOUr1.no.ce.sims.LRT <- lapply(modelOUr1.no.ce.sims, function(x){
###   
###   # If there are any values within the simulation that are outside of valid range make them equal to the limit
###   x.rm.lim <- x
###   x.rm.lim[x.rm.lim<(asin.limits[1])] <- 0
###   x.rm.lim[x.rm.lim>(asin.limits[2])] <- asin(sqrt(1))
###   
###   # Creates temp model with co-evolution
###  temp.modelOUr1.ce <- mvOU(tree = tree.phylo, data = x, echo = F, diagnostic = F)
###  # Creates temp model with no co-evolution
###  temp.modelOUr1.no.ce <- mvOU(tree = tree.phylo, data = x, param = list(sigma="constraint", alpha="constraint"), echo = F,diagnostic = F)
###  
###  temp.modelOUr1.ce.lim <- mvOU(tree = tree.phylo, data = x.rm.lim, echo = F, diagnostic = F)
###  # Creates temp model with no co-evolution
###  temp.modelOUr1.no.ce.lim <- mvOU(tree = tree.phylo, data = x.rm.lim, param = list(sigma="constraint", alpha="constraint"), echo = F,diagnostic = F)
###  
###  # Compares models to determin if "correct" model is returned
###  temp.LWD <- LRT(temp.modelOUr1.ce, temp.modelOUr1.no.ce, echo = F)
###  temp.LWD.lim <- LRT(temp.modelOUr1.ce.lim, temp.modelOUr1.no.ce.lim, echo = F)
###  return(data.frame(LRT = temp.LWD$ratio, LRT.lim = temp.LWD.lim$ratio,
###                    corr = cov2cor(stationary(temp.modelOUr1.ce))[1,2], corr.lim = cov2cor(stationary(temp.modelOUr1.ce.lim))[1,2]))
### })

### modelOUr1.no.ce.sims.LRT <- do.call("rbind", modelOUr1.no.ce.sims.LRT)
### modelOUr1.no.ce.sims.LRT$exceed.limits <- sim.exceed.limits

# Save output to avoid rerunning
### saveRDS(modelOUr1.no.ce.sims.LRT, file = "data/models/modelOUr1.no.ce.sims.LRT.cov2cor.rds")
modelOUr1.no.ce.sims.LRT <- readRDS("data/models/modelOUr1.no.ce.sims.LRT.cov2cor.rds")

hist(modelOUr1.no.ce.sims.LRT$LRT, freq = T, ylim = c(0,400), breaks = 20)
abline(v = LRT(modelOUr1.no.ce, modelOUr1.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr1.ce, modelOUr1.no.ce, echo = F)$ratio, 350, cex = 4, pch = 19)

# Number of simulated LRT that were greater than observed LRT
LRT.no.ce <- table(modelOUr1.no.ce.sims.LRT$LRT.lim>LRT(modelOUr1.ce, modelOUr1.no.ce, echo = F)$ratio)/nsims*100

max.y <- 250
text.size = 20

p <- ggplot(modelOUr1.no.ce.sims.LRT[modelOUr1.no.ce.sims.LRT$LRT.lim>=0,]) +
  geom_histogram(aes(LRT.lim, 
                     fill = as.character(exceed.limits)), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr1.no.ce, modelOUr1.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr1.no.ce, modelOUr1.ce, echo = F)$ratio, y = max.y*0.8), size = 8, ) +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  theme_classic() +
  ylim(c(0, max.y)) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("Co-evolution") +
  labs(subtitle = paste0(LRT.no.ce[2], "% greater than the observed models")) +
  theme(text = element_text(size = text.size), title = element_text(size = text.size*0.75))
p
pp <- ggplot(modelOUr1.no.ce.sims.LRT[modelOUr1.no.ce.sims.LRT$LRT.lim>=0,]) +
  geom_histogram(aes(corr.lim, fill = as.character(exceed.limits)), binwidth = 0.05, show.legend = F) + 
  geom_vline(xintercept = cov2cor(stationary(modelOUr1.ce))[1,2]) +
  geom_point(aes(x = cov2cor(stationary(modelOUr1.ce))[1,2], y = 70*0.8), size = 6, shape = 17) +
  scale_fill_manual(values = c("orange", "orangered")) +
  theme_classic() +
  ylim(c(0, 70)) +
  xlim(-1,1) +
  xlab("Standardised stationary covariance") +
  ylab("Frequency") +
  theme(text = element_text(size = text.size*0.5), title = element_text(size = text.size*0.5))
pp
#variance
p <- p + inset_element(pp, align_to = "plot", left = 0.4, right = 0.95, bottom = 0.4, top = 0.95, ignore_tag = F)
p

##### Pac Atl regimes ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simluate data using mvSIM to test the power of our data to compare between different pacific and atlantic selective regemes
# Model data without different regimes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

model.sims.OUr0.no.ce <- mvSIM(tree = tree.phylo, nsim = nsims, model = c("OU1"), param = list(theta = modelOUr1.no.ce$theta,
                                                                                               alpha = modelOUr1.no.ce$alpha,
                                                                                               sigma = modelOUr1.no.ce$sigma))

# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lapply(model.sims.OUr0.no.ce, function(x) points(x, col = as.factor(row.names(x))))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)

points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# Does simulations include values outside of valid phenotype arc sin of the square root of 0 to 1)
sim.exceed.limits <- lapply(model.sims.OUr0.no.ce, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
table(unlist(sim.exceed.limits))

# Test if model is true model is rejected
### sim.OU.r1.r2P.LRT <- lapply(model.sims.OUr0.no.ce, function(x, i = 1) {
###   
###   # If there are any values within the simulation that are outside of valid range make them equal to the limit
###   x.rm.lim <- x
###   x.rm.lim[x.rm.lim<(asin.limits[1])] <- 0
###   x.rm.lim[x.rm.lim>(asin.limits[2])] <- asin(sqrt(1))
###   
###   # Creates temp model
###   temp.model.mvOUM <- mvOU(tree = tree.phylo.r2.P, model = "OUM", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOU1 <- mvOU(tree = tree.phylo, model = "OU1", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   
###   # Creates temp model
###   temp.model.mvOUM.lim <- mvOU(tree = tree.phylo.r2.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOU1.lim <- mvOU(tree = tree.phylo, model = "OU1", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   # Compares models to determin if "correct" model is returned
###   temp.LWD <- LRT(temp.model.mvOUM, temp.model.mvOU1, echo = F)
###   temp.LWD.lim <- LRT(temp.model.mvOUM.lim, temp.model.mvOU1.lim, echo = F)
###   return(data.frame(LRT = temp.LWD$ratio, LRT.lim = temp.LWD.lim$ratio))
### })
### 
### sim.OU.r1.r2P.LRT <- do.call("rbind", sim.OU.r1.r2P.LRT)
### sim.OU.r1.r2P.LRT$exceed.limits <- sim.exceed.limits
### saveRDS(sim.OU.r1.r2P.LRT, file = "data/models/sim.OU.r1.r2P.LRT.rds")
sim.OU.r1.r2P.LRT <- readRDS("data/models/sim.OU.r1.r2P.LRT.rds")
hist(sim.OU.r1.r2P.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr1.no.ce, modelOUr2.P.no.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr1.no.ce, modelOUr2.P.no.ce, echo = F)$ratio, 350, cex = 4, pch = 19)

# Number of simulated LRT that were greater than observed LRT
LRT.r1.no.ce <- table(sim.OU.r1.r2P.LRT$LRT.lim>LRT(modelOUr1.no.ce, modelOUr2.P.no.ce, echo = F)$ratio)/nsims*100

q <- ggplot(sim.OU.r1.r2P.LRT[sim.OU.r1.r2P.LRT$LRT.lim>0,]) +
  geom_histogram(aes(LRT.lim, 
                     fill = as.character(sim.OU.r1.r2P.LRT$exceed.limits[sim.OU.r1.r2P.LRT$LRT.lim>0])), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr1.no.ce, modelOUr2.P.no.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr1.no.ce, modelOUr2.P.no.ce, echo = F)$ratio, y = max.y*0.8), size = 8) +
  theme_classic() +
  ylim(c(0, max.y)) +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("Pacific and Atlantic regimes") +
  labs(subtitle = paste0(LRT.r1.no.ce[2], "% greater than the observed models")) +
  theme(text = element_text(size = text.size), title = element_text(size = text.size*0.75))

pq <- p + q + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
pq
## ggsave("plots/LDR_coevolve_and_regeimes.png", plot = pq, width = 12, height = 5)


##### HL regime ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simluate data using mvSIM to test the power of our data to compare between different HL and all other regemes
# Model data without different regimes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

set.seed(52398)
model.sims.OUr0.no.ce <- mvSIM(tree = tree.phylo, nsim = nsims, model = c("OU1"), param = list(theta = modelOUr1.no.ce$theta,
                                                                                               alpha = modelOUr1.no.ce$alpha,
                                                                                               sigma = modelOUr1.no.ce$sigma))

# Plot simulated data and true data

plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lapply(model.sims.OUr0.no.ce, function(x) points(x, col = as.factor(row.names(x))))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)

points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# Does simulations include values outside of valid phenotype arc sin of the square root of 0 to 1)
sim.exceed.limits <- lapply(model.sims.OUr0.no.ce, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
table(unlist(sim.exceed.limits))

# Test if model is true model is rejected
### sim.OU.r1.r2HL.LRT <- lapply(model.sims.OUr0.no.ce, function(x, i = 1) {
###   # If there are any values within the simulation that are outside of valid range make them equal to the limit
###   x.rm.lim <- x
###   x.rm.lim[x.rm.lim<(asin.limits[1])] <- 0
###   x.rm.lim[x.rm.lim>(asin.limits[2])] <- asin(sqrt(1))
###   
###   # Creates temp model
###   temp.model.mvOUM <- mvOU(tree = tree.phylo.r2.HL, model = "OUM", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOU1 <- mvOU(tree = tree.phylo, model = "OU1", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   # Limited sim data
###   temp.model.mvOUM.lim <- mvOU(tree = tree.phylo.r2.HL, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOU1.lim <- mvOU(tree = tree.phylo, model = "OU1", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   
###   # Compares models to determin if "correct" model is returned
###   temp.LWD <- LRT(temp.model.mvOUM, temp.model.mvOU1, echo = F)
###   temp.LWD.lim <- LRT(temp.model.mvOUM.lim, temp.model.mvOU1.lim, echo = F)
###   return(data.frame(LRT = temp.LWD$ratio, LRT.lim = temp.LWD.lim$ratio))
### })

### sim.OU.r1.r2HL.LRT <- do.call("rbind", sim.OU.r1.r2HL.LRT)
### sim.OU.r1.r2HL.LRT$exceed.limits <- sim.exceed.limits
### saveRDS(sim.OU.r1.r2HL.LRT, file = "data/models/sim.OU.r1.r2HL.LRT.rds")

sim.OU.r1.r2HL.LRT <- readRDS("data/models/sim.OU.r1.r2HL.LRT.rds")

hist(sim.OU.r1.r2HL.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr1.no.ce, modelOUr2.HL.no.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr1.no.ce, modelOUr2.HL.no.ce, echo = F)$ratio, 350, cex = 4, pch = 19)

# Number of simulated LRT that were greater than observed LRT
# any LL lower than 0 
LL.neg <- which(!sim.OU.r1.r2HL.LRT$LRT.lim<0)
LRT.r2.HL.no.ce <- table(sim.OU.r1.r2HL.LRT$LRT.lim[LL.neg] > LRT(modelOUr1.no.ce, modelOUr2.HL.no.ce, echo = F)$ratio)/length(LL.neg)*100

## Are any values of LR lower than zero?
sim.OU.r1.r2HL.LRT$LRT<0|sim.OU.r1.r2HL.LRT$LRT.lim<0

v <- ggplot(sim.OU.r1.r2HL.LRT[sim.OU.r1.r2HL.LRT$LRT.lim>0,]) +
  geom_histogram(aes(LRT.lim, 
                     fill = as.character(exceed.limits)), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr1.no.ce, modelOUr2.HL.no.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr1.no.ce, modelOUr2.HL.no.ce, echo = F)$ratio, y = 200*0.8), size = 8) +
  theme_classic() +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  ylim(c(0, 200)) +
  #xlim(c(-0, max(sim.OU.r1.r2HL.LRT$LRT.lim))) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("High latitude regime") +
  labs(subtitle = paste0(round(LRT.r2.HL.no.ce[2], digits = 1), "% greater than the observed models")) +
  theme(text = element_text(size = text.size), title = element_text(size = text.size*0.75))
v
pqv <- (p+q)/v + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
pqv
### ggsave(plot = (p+q)/v, filename= "plots/LDR_coevolve_and_regeimes.png", width = 12, height = 10)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
##### Pac Atl and HL regimes ####
# Simluate data using mvSIM to test the power of our data to compare between different HL and all other regemes
# Model data without different regimes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

model.sims.OUr0.no.ce <- mvSIM(tree = tree.phylo, nsim = nsims, model = c("OU1"),
                               param = list(theta = modelOUr1.no.ce$theta,
                                            alpha = modelOUr1.no.ce$alpha,
                                            sigma = modelOUr1.no.ce$sigma))

# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lapply(model.sims.OUr0.no.ce, function(x) points(x, col = as.factor(row.names(x))))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# Does simulations include values outside of valid phenotype arc sin of the square root of 0 to 1)
sim.exceed.limits <- lapply(model.sims.OUr0.no.ce, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
table(unlist(sim.exceed.limits))

# Test if model is true model is rejected
### sim.OU.r1.r3.HL.P.LRT <- lapply(model.sims.OUr0.no.ce, function(x, i = 1) {
###   # If there are any values within the simulation that are outside of valid range make them equal to the limit
###   x.rm.lim <- x
###   x.rm.lim[x.rm.lim<(asin.limits[1])] <- 0
###   x.rm.lim[x.rm.lim>(asin.limits[2])] <- asin(sqrt(1))
###   # Creates temp model
###   temp.model.mvOUM <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOU1 <- mvOU(tree = tree.phylo, model = "OU1", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   # Limit models
###   temp.model.mvOUM.lim <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOU1.lim <- mvOU(tree = tree.phylo, model = "OU1", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   
###   # Compares models to determin if "correct" model is returned
###   temp.LWD <- LRT(temp.model.mvOUM, temp.model.mvOU1, echo = F)
###   temp.LWD.lim <- LRT(temp.model.mvOUM.lim, temp.model.mvOU1.lim, echo = F)
###   return(data.frame(LRT = temp.LWD$ratio, LRT.lim = temp.LWD.lim$ratio))
### })
### 
### sim.OU.r1.r3.HL.P.LRT <- do.call("rbind", sim.OU.r1.r3.HL.P.LRT)
### sim.OU.r1.r3.HL.P.LRT$exceeds.lim <- unlist(sim.exceed.limits)
### 
### saveRDS(sim.OU.r1.r3.HL.P.LRT, file = "data/models/sim.OU.r1.r3.HL.P.LRT.rds")
sim.OU.r1.r3.HL.P.LRT <- readRDS("data/models/sim.OU.r1.r3.HL.P.LRT.rds")
hist(sim.OU.r1.r3.HL.P.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, 350, cex = 4, pch = 19)

# Number of simulated LRT that were greater than observed LRT
LRT.OU.r1.r3.HL.P <- table(sim.OU.r1.r3.HL.P.LRT$LRT.lim > LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/nsims*100

sim.OU.r1.r3.HL.P.LRT$exceeds.lim
t <- ggplot(sim.OU.r1.r3.HL.P.LRT[sim.OU.r1.r3.HL.P.LRT$LRT.lim>0,]) +
  geom_histogram(aes(LRT.lim,
                     fill = sim.OU.r1.r3.HL.P.LRT$exceeds.lim[sim.OU.r1.r3.HL.P.LRT$LRT.lim>0]), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, y = 100*0.8), size = 8) +
  theme_classic() +
  ylim(c(0, 100)) +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("High latitude, Pacific, and Atlantic regimes") +
  labs(subtitle = paste0(LRT.OU.r1.r3.HL.P[2], "% greater than the observed models")) +
  theme(text = element_text(size = text.size), title = element_text(size = text.size*0.75))

t
pqvt <- (p+q)/(v+t) + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
pqvt
ggsave(plot = pqvt , filename= "plots/LDR_coevolve_and_regeimes.png", width = 12, height = 12)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
##### Pac Atl or HL regimes ####
# Simluate data using mvSIM to test the power of our data to compare between different HL and and separate Pacific and Atlantic regimes
# Model data without different regimes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Simulate evolution with a different selective regime between HL and all other lineages
model.sims.OUr2.HL.no.ce <- mvSIM(tree = tree.phylo.r2.HL, 
                                  nsim = nsims, model = c("OUM"), 
                                  param = list(theta = modelOUr2.HL.no.ce$theta,
                                               alpha = modelOUr2.HL.no.ce$alpha,
                                               sigma = modelOUr2.HL.no.ce$sigma))

# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lapply(model.sims.OUr2.HL.no.ce, function(x) points(x, col = as.factor(row.names(x))))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# Does simulations include values outside of valid phenotype arc sin of the square root of 0 to 1)
sim.exceed.limits <- lapply(model.sims.OUr2.HL.no.ce, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
table(unlist(sim.exceed.limits))

# Calculate the likelihood ratio models on the simlations
### sim.OU.r2.r3.HL.P.LRT <- lapply(model.sims.OUr2.HL.no.ce, function(x, i = 1) {
###   # If there are any values within the simulation that are outside of valid range make them equal to the limit
###     x.rm.lim <- x
###   x.rm.lim[x.rm.lim<(asin.limits[1])] <- 0
###   x.rm.lim[x.rm.lim>(asin.limits[2])] <- asin(sqrt(1))
###   # Creates temp model
###   temp.model.mvOUMr3 <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOUMr2 <- mvOU(tree = tree.phylo.r2.HL, model = "OUM", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   
###   #Model with data with limits
###   temp.model.mvOUMr3.lim <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOUMr2.lim <- mvOU(tree = tree.phylo.r2.HL, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   
###   # Compares models to determin if "correct" model is returned
###   temp.LWD <- LRT(temp.model.mvOUMr3, temp.model.mvOUMr2, echo = F)
###   temp.LWD.lim <- LRT(temp.model.mvOUMr3.lim, temp.model.mvOUMr2.lim, echo = F)
###   return(data.frame(LRT = temp.LWD$ratio, LRT.lim = temp.LWD.lim$ratio))
### })
### 
### sim.OU.r2.r3.HL.P.LRT <- do.call("rbind", sim.OU.r2.r3.HL.P.LRT)
### sim.OU.r2.r3.HL.P.LRT$exceeds.limits <- unlist(sim.exceed.limits)
### 
### saveRDS(sim.OU.r2.r3.HL.P.LRT, file = "data/models/sim.OU.r2.r3.HL.P.LRT.rds")

sim.OU.r2.r3.HL.P.LRT <- readRDS("data/models/sim.OU.r2.r3.HL.P.LRT.rds")
hist(sim.OU.r2.r3.HL.P.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, 350, cex = 4, pch = 19)

LL.neg.r2.r3.HL.P.LRT <- which(!(sim.OU.r2.r3.HL.P.LRT$LRT<0|sim.OU.r2.r3.HL.P.LRT$LRT.lim<0))
LRT.OU.r2.r3.HL.P <- table(sim.OU.r2.r3.HL.P.LRT$LRT[LL.neg.r2.r3.HL.P.LRT]>LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/length(LL.neg.r2.r3.HL.P.LRT)*100
LRT.OU.r2.r3.HL.P
z <- ggplot(sim.OU.r2.r3.HL.P.LRT[sim.OU.r2.r3.HL.P.LRT$LRT.lim>0,]) +
  geom_histogram(aes(LRT.lim, fill = as.character(exceeds.limits)), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio, y = 120*0.8), size = 8) +
  theme_classic() +
  ylim(c(0, 120)) +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("Pacific and Atlantic regimes with High Latitude") +
  labs(subtitle = paste0(LRT.OU.r2.r3.HL.P[2], "% greater than the observed models")) +
  theme(text = element_text(size = text.size), title = element_text(size = text.size*0.75))
z

z + inset_element(t, align_to = "plot", left = 0.5, right = 0.95, bottom = 0.4, top = 0.95, ignore_tag = T, on_top = T)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
##### Co-evolution with multiple regimes ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simulate evolution with a different selective regime between HL and all other lineages
model.sims.OUr3.HL.P.no.ce <- mvSIM(tree = tree.phylo.r3.HL.P, 
                                    nsim = nsims, model = c("OUM"), 
                                    param = list(theta = modelOUr3.HL.P.no.ce$theta,
                                                 alpha = modelOUr3.HL.P.no.ce$alpha,
                                                 sigma = modelOUr3.HL.P.no.ce$sigma))



# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
lapply(model.sims.OUr3.HL.P.no.ce, function(x) points(x, col = as.factor(row.names(x))))
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# Does simulations include values outside of valid phenotype arc sin of the square root of 0 to 1)
sim.exceed.limits <- lapply(model.sims.OUr3.HL.P.no.ce, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
table(unlist(sim.exceed.limits))

## Actual correlation
cov2cor(stationary(modelOUr3.HL.P))[1,2]
# Calculate the likelihood ratio models on the simlations
### sim.OU.r3.HL.P.ce.LRT <- lapply(model.sims.OUr3.HL.P.no.ce , function(x, i = 1) {
###   # If there are any values within the simulation that are outside of valid range make them equal to the limit
###   x.rm.lim <- x
###   x.rm.lim[x.rm.lim<(asin.limits[1])] <- 0
###   x.rm.lim[x.rm.lim>(asin.limits[2])] <- asin(sqrt(1))
###   # Creates temp model
###   temp.model.mvOUMr3.no.ce <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOUMr3.ce <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x, echo = F, diagnostic = F)
###   # Temp model with data set to limits
###   temp.model.mvOUMr3.no.ce.lim <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, param = list(sigma="constraint", alpha="constraint"))
###   temp.model.mvOUMr3.ce.lim <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F)
###   # Compares models to determine if "correct" model is returned
###   temp.LWD <- LRT(temp.model.mvOUMr3.ce, temp.model.mvOUMr3.no.ce, echo = F)
###   temp.LWD.lim <- LRT(temp.model.mvOUMr3.ce.lim, temp.model.mvOUMr3.no.ce.lim, echo = F)
###   temp.corr <- cov2cor(stationary(temp.model.mvOUMr3.ce))
###   temp.corr.lim <- cov2cor(stationary(temp.model.mvOUMr3.ce.lim))
###   return(data.frame(LRT = temp.LWD$ratio, LRT.lim = temp.LWD.lim$ratio, corr = temp.corr[1,2], corr.lim = temp.corr.lim[1,2]))
### })
### 
### sim.OU.r3.HL.P.ce.LRT <- do.call("rbind", sim.OU.r3.HL.P.ce.LRT)
### sim.OU.r3.HL.P.ce.LRT$exceeds.limits <- unlist(sim.exceed.limits)
### saveRDS(sim.OU.r3.HL.P.ce.LRT, file = "data/models/sim.OU.r3.HL.P.ce.LRT.rds")

sim.OU.r3.HL.P.ce.LRT <- readRDS("data/models/sim.OU.r3.HL.P.ce.LRT.rds")
hist(sim.OU.r3.HL.P.ce.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio, 350, cex = 4, pch = 19)
summary(sim.OU.r3.HL.P.ce.LRT$LRT)

LRT.OU.r3.ce <- table(sim.OU.r3.HL.P.ce.LRT$LRT.lim>LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio)/nsims*100

k <- ggplot() +
  geom_histogram(aes(sim.OU.r3.HL.P.ce.LRT$LRT.lim, 
                     fill = as.character(modelOUr1.no.ce.sims.LRT$exceed.limits)), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio, y = (100)*0.8), size = 8) +
  theme_classic() +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  ylim(c(0, 100)) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("Co-evolution with High Latitude, Pacific, and Atlantic regimes") +
  labs(subtitle = paste0(LRT.OU.r3.ce[2], "% greater than the observed models")) +
  theme(text = element_text(size = text.size), title = element_text(size = text.size*0.75))

table(sim.OU.r3.HL.P.ce.LRT$corr.lim>cov2cor(stationary(modelOUr3.HL.P))[1,2])/nsims*100
table(abs(sim.OU.r3.HL.P.ce.LRT$corr.lim)<abs(cov2cor(stationary(modelOUr3.HL.P))[1,2]))/nsims*100

kk <- ggplot() +
  geom_histogram(aes(sim.OU.r3.HL.P.ce.LRT$corr.lim, fill = as.character(modelOUr1.no.ce.sims.LRT$exceed.limits)), binwidth = 0.05, show.legend = F) + 
  geom_vline(xintercept = cov2cor(stationary(modelOUr3.HL.P))[1,2]) +
  geom_point(aes(x = cov2cor(stationary(modelOUr3.HL.P))[1,2], y = 50*0.8), size = 6, shape = 17) +
  theme_classic() +
  scale_fill_manual(values = c("orange", "orangered")) +
  ylim(c(0, 50)) +
  xlim(-1,1) +
  xlab("Standardised stationary covariance") +
  ylab("Frequency") +
  theme(text = element_text(size = text.size*0.5), title = element_text(size = text.size*0.5))

k <- k + inset_element(kk, align_to = "plot", left = 0.5, right = 0.95, bottom = 0.4, top = 0.95, ignore_tag = F)
pqvta <- ((p+q)/(v+t)/k)
pqvta <- pqvta + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
pqvta

ggsave(plot = pqvta, filename= "plots/LDR_coevolve_and_regeimes_Reg_cov_wlims.png", width = 12, height = 15)

## Model Statistics
stats.df <- matrix(NA, nrow = 5, 3)

# Coevolution
stats.df[1,1] <- "Single regime with co-evolution"
stats.df[1,2] <- (table(modelOUr1.no.ce.sims.LRT$LRT>LRT(modelOUr1.ce, modelOUr1.no.ce, echo = F)$ratio)/nsims*100)[2]
stats.df[1,3] <- (table(modelOUr1.no.ce.sims.LRT$LRT.lim>LRT(modelOUr1.ce, modelOUr1.no.ce, echo = F)$ratio)/nsims*100)[2]

# Pac and Atlantic
stats.df[2,1] <- "Pacfic and Atlantic regime"
stats.df[2,2] <- (table(sim.OU.r1.r2P.LRT$LRT>LRT(modelOUr1.no.ce, modelOUr2.P.no.ce, echo = F)$ratio)/nsims*100)[2]
stats.df[2,3] <- (table(sim.OU.r1.r2P.LRT$LRT.lim>LRT(modelOUr1.no.ce, modelOUr2.P.no.ce, echo = F)$ratio)/nsims*100)[2]

# HL regime
LL.neg..r1.r2HL.LRT <- which(!(sim.OU.r1.r2HL.LRT$LRT.lim<0|sim.OU.r1.r2HL.LRT$LRT<0))
stats.df[3,1] <- "High latitude regime"
stats.df[3,2] <- (table(sim.OU.r1.r2HL.LRT$LRT[LL.neg..r1.r2HL.LRT] > LRT(modelOUr1.no.ce, modelOUr2.HL.no.ce, echo = F)$ratio)/length(LL.neg..r1.r2HL.LRT)*100)[2]
stats.df[3,3] <- (table(sim.OU.r1.r2HL.LRT$LRT.lim[LL.neg..r1.r2HL.LRT] > LRT(modelOUr1.no.ce, modelOUr2.HL.no.ce, echo = F)$ratio)/length(LL.neg..r1.r2HL.LRT)*100)[2]

# HL, Pac, Atl
LL.neg.r1.r3.HL.P.LRT <- which(!(sim.OU.r1.r3.HL.P.LRT$LRT<0|sim.OU.r1.r3.HL.P.LRT$LRT.lim<0))

stats.df[4,1] <- "Pacific Atlantic and High Latitude regime"
stats.df[4,2] <- (table(sim.OU.r1.r3.HL.P.LRT$LRT > LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/length(LL.neg.r1.r3.HL.P.LRT)*100)[2]
stats.df[4,3] <- (table(sim.OU.r1.r3.HL.P.LRT$LRT.lim > LRT(modelOUr1.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/length(LL.neg.r1.r3.HL.P.LRT)*100)[2]

# HL vs Pac and Atl
LL.neg.r2.r3.HL.P.LRT <- which(!(sim.OU.r2.r3.HL.P.LRT$LRT<0|sim.OU.r2.r3.HL.P.LRT$LRT.lim<0))
table(sim.OU.r2.r3.HL.P.LRT$LRT[LL.neg.r2.r3.HL.P.LRT] > LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/length(LL.neg.r2.r3.HL.P.LRT)*100
table(sim.OU.r2.r3.HL.P.LRT$LRT.lim[LL.neg.r2.r3.HL.P.LRT] > LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/length(LL.neg.r2.r3.HL.P.LRT)*100

table(sim.OU.r2.r3.HL.P.LRT$LRT[!sim.OU.r2.r3.HL.P.LRT$exceeds.limits] > LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/length(which(!sim.OU.r2.r3.HL.P.LRT$exceeds.limits))*100
table(sim.OU.r2.r3.HL.P.LRT$LRT.lim[sim.OU.r2.r3.HL.P.LRT$exceeds.limits] > LRT(modelOUr2.HL.no.ce, modelOUr3.HL.P.no.ce, echo = F)$ratio)/length(which(sim.OU.r2.r3.HL.P.LRT$exceeds.limits))*100

# Multiregime coevoltion
stats.df[5,1] <- "Pacific Atlantic and High Latitude regime with co-evolution"
stats.df[5,2] <- (table(sim.OU.r3.HL.P.ce.LRT$LRT>LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio)/nsims*100)[2]
stats.df[5,3] <- (table(sim.OU.r3.HL.P.ce.LRT$LRT.lim>LRT(modelOUr3.HL.P, modelOUr3.HL.P.no.ce, echo = F)$ratio)/nsims*100)[2]

stats.df <- data.frame(stats.df)
colnames(stats.df) <- c("Model comparison", "LRT", "LRT with boundries")
stats.df

write.csv(stats.df, file = "plots/LRT_stats.csv")

# Co-evolution
cov2cor(stationary(modelOUr1.ce))
cov2cor(stationary(modelOUr3.HL.P))

# Stationary variance
stationary(modelOUr3.HL.P.no.ce)
(modelOUr3.HL.P.no.ce$sigma^2)/(2*modelOUr3.HL.P.no.ce$alpha)

# Theta values converted to proporation of wing melansiation
sin(modelOUr3.HL.P.no.ce$theta)^2


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#### Finished ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Bonus analysis from Julien ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

modelOUr3.HL.P.ce.A <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r3.HL.P,
                            param = list(decompSigma="diagonal"))

modelOUr3.HL.P.ce.A.sims <- mvSIM(tree = tree.phylo.r3.HL.P, 
                                  nsim = nsims, model = c("OUM"), 
                                  param = list(theta = modelOUr3.HL.P.ce.A$theta,
                                               alpha = modelOUr3.HL.P.ce.A$alpha,
                                               sigma = modelOUr3.HL.P.ce.A$sigma))
# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
lapply(modelOUr3.HL.P.ce.A.sims, function(x) points(x, col = as.factor(row.names(x))))
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
##### Co-evolution with multiple regimes constrained alpha using decomp ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simulate evolution with a different selective regime between HL and all other lineages

modelOUr3.HL.P.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r3.HL.P, 
                             param = list(decomp = "diagonal", decompSigma="diagonal"))

modelOUr3.HL.P.no.ce.sims <- mvSIM(tree = tree.phylo.r3.HL.P, 
                                   nsim = nsims, model = c("OUM"), 
                                   param = list(theta = modelOUr3.HL.P.no.ce$theta,
                                                alpha = modelOUr3.HL.P.no.ce$alpha,
                                                sigma = modelOUr3.HL.P.no.ce$sigma))

# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
lapply(modelOUr3.HL.P.no.ce.sims, function(x) points(x, col = as.factor(row.names(x))))
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# Does simulations include values outside of valid phenotype arc sin of the square root of 0 to 1)
sim.exceed.limits <- lapply(modelOUr3.HL.P.no.ce.sims, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
table(unlist(sim.exceed.limits))

## Actual correlation
cov2cor(stationary(modelOUr3.HL.P))[1,2]
# Calculate the likelihood ratio models on the simlations
### sims.modelOUr3.HL.P.ce.A.LRT <- lapply(modelOUr3.HL.P.no.ce.sims, function(x, i = 1){
###   # If there are any values within the simulation that are outside of valid range make them equal to the limit
###   # x <- modelOUr3.HL.P.no.ce.sims[[20]]
###   x.rm.lim <- x
###   x.rm.lim[x.rm.lim<(asin.limits[1])] <- 0
###   x.rm.lim[x.rm.lim>(asin.limits[2])] <- asin(sqrt(1))
###   # Creates temp model
###   temp.model.mvOUMr3.no.ce.alp.v <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x, echo = F, diagnostic = F, param = list(decomp = "diagonal", decompSigma="diagonal"))
###   temp.model.mvOUMr3.ce.alp.v <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x, echo = F, diagnostic = F, param = list(decompSigma="diagonal"))
###   # Temp model with data set to limits
###   temp.model.mvOUMr3.no.ce.alp.v.lim <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, param = list(decomp = "diagonal", decompSigma="diagonal"))
###   temp.model.mvOUMr3.ce.alp.v.lim <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, param = list(decompSigma="diagonal"))
###   # Compares models to determine if "correct" model is returned
###   temp.LWD <- LRT(temp.model.mvOUMr3.ce.alp.v, temp.model.mvOUMr3.no.ce.alp.v, echo = F)
###   temp.LWD.lim <- LRT(temp.model.mvOUMr3.ce.alp.v.lim, temp.model.mvOUMr3.no.ce.alp.v.lim, echo = F)
###   temp.corr <- cov2cor(stationary(temp.model.mvOUMr3.ce.alp.v))
###   temp.corr.lim <- cov2cor(stationary(temp.model.mvOUMr3.ce.alp.v.lim))
###   return(data.frame(LRT = temp.LWD$ratio, LRT.lim = temp.LWD.lim$ratio, corr = temp.corr[1,2], corr.lim = temp.corr.lim[1,2]))
### })
### 
### sims.modelOUr3.HL.P.ce.A.LRT <- do.call("rbind", sims.modelOUr3.HL.P.ce.A.LRT)
### sims.modelOUr3.HL.P.ce.A.LRT$exceeds.limits <- unlist(sim.exceed.limits)
### saveRDS(sims.modelOUr3.HL.P.ce.A.LRT, file = "data/models/sims.modelOUr3.HL.P.ce.A.LRT.rds")

sims.modelOUr3.HL.P.ce.A.LRT <- readRDS("data/models/sims.modelOUr3.HL.P.ce.A.LRT.rds")
hist(sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr3.HL.P.no.ce, modelOUr3.HL.P.ce.A, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr3.HL.P.no.ce, modelOUr3.HL.P.ce.A, echo = F)$ratio, 350, cex = 4, pch = 19)
summary(sims.modelOUr3.HL.P.ce.A.LRT$LRT)

table(sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim<0)

#LRT.OU.r3.ce.A <- table(sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim>LRT(modelOUr3.HL.P.ce.A, modelOUr3.HL.P.no.ce, echo = F)$ratio)/nsims*100
LRT.OU.r3.ce.A <- table((sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim>LRT(modelOUr3.HL.P.ce.A, modelOUr3.HL.P.no.ce, echo = F)$ratio)[!sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim<0])
LRT.OU.r3.ce.A <- round(LRT.OU.r3.ce.A/sum(LRT.OU.r3.ce.A)*100, 1)

j <- ggplot() +
  geom_histogram(aes(sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim[sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim>0], 
                     fill = as.character(sims.modelOUr3.HL.P.ce.A.LRT$exceeds.limits[sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim>0])), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr3.HL.P.ce.A, modelOUr3.HL.P.no.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr3.HL.P.ce.A, modelOUr3.HL.P.no.ce, echo = F)$ratio, y = (200)*0.8), size = 8) +
  theme_classic() +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  ylim(c(0, 200)) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("Co-evolution with High Latitude, Pacific, and Atlantic regimes") +
  labs(subtitle = paste0(LRT.OU.r3.ce.A[2], "% greater than the observed models")) +
  theme(text = element_text(size = text.size), title = element_text(size = text.size*0.75))

table(sim.OU.r3.HL.P.ce.LRT$corr.lim>cov2cor(stationary(modelOUr3.HL.P))[1,2])/nsims*100
table(abs(sim.OU.r3.HL.P.ce.LRT$corr.lim)<abs(cov2cor(stationary(modelOUr3.HL.P))[1,2]))/nsims*100

jj <- ggplot() +
  geom_histogram(aes(sims.modelOUr3.HL.P.ce.A.LRT$corr.lim, fill = as.character(sims.modelOUr3.HL.P.ce.A.LRT$exceeds.limits)), binwidth = 0.05, show.legend = F) + 
  geom_vline(xintercept = cov2cor(stationary(modelOUr3.HL.P.ce.A))[1,2]) +
  geom_point(aes(x = cov2cor(stationary(modelOUr3.HL.P.ce.A))[1,2], y = 50*0.8), size = 6, shape = 17) +
  theme_classic() +
  scale_fill_manual(values = c("orange", "orangered")) +
  ylim(c(0, 50)) +
  xlim(-1,1) +
  xlab("Standardised stationary covariance") +
  ylab("Frequency") +
  theme(text = element_text(size = text.size*0.5), title = element_text(size = text.size*0.5))

j <- j + inset_element(jj, align_to = "plot", left = 0.5, right = 0.95, bottom = 0.4, top = 0.95, ignore_tag = F, on_top = T)

pqvtj <- ((p+q)/(v+t)/j)
pqvtj <- pqvtj + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
pqvtj

ggsave(plot = pqvtj, filename= "plots/LDR_coevolve_and_regeimes_Reg_cov_in_alpha_wlims_.png", width = 12, height = 15)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
##### Co-evolution with multiple regimes constrained diagonal alpha using decomp ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simulate evolution with a different selective regime between HL and all other lineages

alpha.set.diagonal <- matrix(c(1,NA,NA,1), nrow = 2)
sigma.set.diagonal <- matrix(c(1,NA,NA,1), nrow = 2)
# Model with same diagonal plots and no co-evolution (NA diagonal alpha)
modelOUr3.HL.P.sda.no.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r3.HL.P, 
                             param = list(decomp = alpha.set.diagonal, decompSigma=sigma.set.diagonal))
# Model with same diagonal alpha values but with co-evolution
modelOUr3.HL.P.sda.ce <- mvOU(data = as.matrix(var.mean), error = as.matrix(var.sq.se), tree = tree.phylo.r3.HL.P, 
                                 param = list(decomp = matrix(c(1,2,2,1), nrow = 2), decompSigma = sigma.set.diagonal))


modelOUr3.HL.P.sda.no.ce.sims <- mvSIM(tree = tree.phylo.r3.HL.P, 
                                   nsim = nsims, model = c("OUM"), 
                                   param = list(theta = modelOUr3.HL.P.sda.no.ce$theta,
                                                alpha = modelOUr3.HL.P.sda.no.ce$alpha,
                                                sigma = modelOUr3.HL.P.sda.no.ce$sigma))

# Plot simulated data and true data
plot(NA, NA, xlab = "peak", ylab = "off-peak", xlim = asin.limits+c(-0.5,0.5), ylim=asin.limits+c(-0.5,0.5))
lines(c(asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1],asin.limits[1]),
      c(asin.limits[1],asin.limits[1],asin.limits[2],asin.limits[2],asin.limits[1]), lwd = 4)
lapply(modelOUr3.HL.P.sda.no.ce.sims, function(x) points(x, col = as.factor(row.names(x))))
points(var.mean, cex = 5, pch = 4, lwd = 10, col = as.factor(row.names(var.mean)))

# Does simulations include values outside of valid phenotype arc sin of the square root of 0 to 1)
sim.exceed.limits <- lapply(modelOUr3.HL.P.sda.no.ce.sims, function(x) any(x<asin.limits[1]|x>asin.limits[2]))
table(unlist(sim.exceed.limits))

## Actual correlation
cov2cor(stationary(modelOUr3.HL.P.sda.ce))[1,2]
# Calculate the likelihood ratio models on the simlations
sims.modelOUr3.HL.P.sda.ce.A.LRT <- lapply(modelOUr3.HL.P.sda.no.ce.sims, function(x, i = 1){
 # If there are any values within the simulation that are outside of valid range make them equal to the limit
 # x <- modelOUr3.HL.P.sda.no.ce.sims[[20]]
 x.rm.lim <- x
 x.rm.lim[x.rm.lim<(asin.limits[1])] <- 0
 x.rm.lim[x.rm.lim>(asin.limits[2])] <- asin(sqrt(1))
 # Creates temp model
 temp.model.mvOUMr3.no.ce.sda <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x, echo = F, diagnostic = F, 
                                        param = list(decomp = alpha.set.diagonal, decompSigma=sigma.set.diagonal))
 temp.model.mvOUMr3.ce.sda <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x, echo = F, diagnostic = F, 
                                     param = list(decomp=matrix(c(1,2,2,1), nrow = 2), decompSigma=sigma.set.diagonal))
 # Temp model with data set to limits
 temp.model.mvOUMr3.no.ce.sda.lim <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, 
                                      param = list(decomp = alpha.set.diagonal, decompSigma=sigma.set.diagonal))
 temp.model.mvOUMr3.ce.sda.lim <- mvOU(tree = tree.phylo.r3.HL.P, model = "OUM", data = x.rm.lim, echo = F, diagnostic = F, 
                                   param = list(decomp=matrix(c(1,2,2,1), nrow = 2), decompSigma=sigma.set.diagonal))
 
 # Compares models to determine if "correct" model is returned
 temp.LWD <- LRT(temp.model.mvOUMr3.no.ce.sda, temp.model.mvOUMr3.ce.sda, echo = F)
 temp.LWD.lim <- LRT(temp.model.mvOUMr3.no.ce.sda.lim, temp.model.mvOUMr3.ce.sda.lim, echo = F)
 temp.corr <- cov2cor(stationary(temp.model.mvOUMr3.ce.sda))
 temp.corr.lim <- cov2cor(stationary(temp.model.mvOUMr3.ce.sda.lim))
 return(data.frame(LRT = temp.LWD$ratio, LRT.lim = temp.LWD.lim$ratio, corr = temp.corr[1,2], corr.lim = temp.corr.lim[1,2]))
})

sims.modelOUr3.HL.P.sda.ce.A.LRT <- do.call("rbind", sims.modelOUr3.HL.P.sda.ce.A.LRT)
sims.modelOUr3.HL.P.sda.ce.A.LRT$exceeds.limits <- unlist(sim.exceed.limits)
saveRDS(sims.modelOUr3.HL.P.sda.ce.A.LRT, file = "data/models/sims.modelOUr3.HL.P.sda.ce.A.LRT.rds")

sims.modelOUr3.HL.P.sda.ce.A.LRT <- readRDS("data/models/sims.modelOUr3.HL.P.sda.ce.A.LRT.rds")
hist(sims.modelOUr3.HL.P.sda.ce.A.LRT$LRT.lim, freq = T, breaks = 20,  ylim = c(0, 400))
abline(v = LRT(modelOUr3.HL.P.sda.no.ce, modelOUr3.HL.P.sda.ce, echo = F)$ratio, lwd  = 3)
points(LRT(modelOUr3.HL.P.sda.no.ce, modelOUr3.HL.P.sda.ce, echo = F)$ratio, 350, cex = 4, pch = 19)
summary(sims.modelOUr3.HL.P.sda.ce.A.LRT$LRT)

table(sims.modelOUr3.HL.P.sda.ce.A.LRT$LRT.lim<0)

#LRT.OU.r3.ce.A <- table(sims.modelOUr3.HL.P.ce.A.LRT$LRT.lim>LRT(modelOUr3.HL.P.ce.A, modelOUr3.HL.P.no.ce, echo = F)$ratio)/nsims*100
LRT.OU.r3.ce.A <- table((sims.modelOUr3.HL.P.sda.ce.A.LRT$LRT.lim>LRT(modelOUr3.HL.P.sda.no.ce, modelOUr3.HL.P.sda.ce, echo = F)$ratio)[!sims.modelOUr3.HL.P.sda.ce.A.LRT$LRT.lim<0])
LRT.OU.r3.ce.A <- round(LRT.OU.r3.ce.A/sum(LRT.OU.r3.ce.A)*100, 1)

f <- ggplot() +
  geom_histogram(aes(sims.modelOUr3.HL.P.sda.ce.A.LRT$LRT.lim[sims.modelOUr3.HL.P.sda.ce.A.LRT$LRT.lim>0], 
                     fill = as.character(sims.modelOUr3.HL.P.sda.ce.A.LRT$exceeds.limits[sims.modelOUr3.HL.P.sda.ce.A.LRT$LRT.lim>0])), binwidth = 0.5, show.legend = F) + 
  geom_vline(xintercept = LRT(modelOUr3.HL.P.sda.no.ce, modelOUr3.HL.P.sda.ce, echo = F)$ratio) +
  geom_point(aes(x = LRT(modelOUr3.HL.P.sda.no.ce, modelOUr3.HL.P.sda.ce, echo = F)$ratio, y = (300)*0.8), size = 8) +
  theme_classic() +
  scale_fill_manual(values = c("skyblue", "deepskyblue")) +
  ylim(c(0, 300)) +
  xlab("likelihood ratio") +
  ylab("Frequency") +
  ggtitle("Co-evolution with High Latitude, Pacific, and Atlantic regimes") +
  labs(subtitle = paste0(LRT.OU.r3.ce.A[2], "% greater than the observed models")) +
  theme(text = element_text(size = text.size), title = element_text(size = text.size*0.75))

table(sim.OU.r3.HL.P.ce.LRT$corr.lim>cov2cor(stationary(modelOUr3.HL.P))[1,2])/nsims*100
table(abs(sim.OU.r3.HL.P.ce.LRT$corr.lim)<abs(cov2cor(stationary(modelOUr3.HL.P))[1,2]))/nsims*100

ff <- ggplot() +
  geom_histogram(aes(sims.modelOUr3.HL.P.sda.ce.A.LRT$corr.lim, fill = as.character(sims.modelOUr3.HL.P.sda.ce.A.LRT$exceeds.limits)), binwidth = 0.05, show.legend = F) + 
  geom_vline(xintercept = cov2cor(stationary(modelOUr3.HL.P.sda.ce))[1,2]) +
  geom_point(aes(x = cov2cor(stationary(modelOUr3.HL.P.sda.ce))[1,2], y = 50*0.8), size = 6, shape = 17) +
  theme_classic() +
  scale_fill_manual(values = c("orange", "orangered")) +
  ylim(c(0, 60)) +
  xlim(-1,1) +
  xlab("Standardised stationary covariance") +
  ylab("Frequency") +
  theme(text = element_text(size = text.size*0.5), title = element_text(size = text.size*0.5))

f <- f + inset_element(ff, align_to = "plot", left = 0.5, right = 0.95, bottom = 0.4, top = 0.95, ignore_tag = F, on_top = T)

pqvtf <- ((p+q)/(v+t)/f)
pqvtf <- pqvtf + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
pqvtf

ggsave(plot = pqvtf, filename= "plots/LDR_coevolve_and_regeimes_Reg_cov_in_alpha_wlims_and_set_alpha_diagonal.png", width = 12, height = 15)


