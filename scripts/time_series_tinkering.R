require("layeranalyzer")
require("data.table")
tit<-read.csv("~/Desktop/titia_data.csv")
tit$hw.logit = 	 log( (tit$proportion.pigmented-0.01) / (1-((tit$proportion.pigmented-0.01))))

tit$two.week.bin = cut(tit$Julian.date,seq(from=1,to=366,by=14)) #making two-week time bins for time series 
tit$two.week.bin.midpoint = as.numeric(sub("\\((.+),.*", "\\1", tit$two.week.bin))+7 
          
          
clusterE<-subset(tit,tip.label=="E")
dt = data.table(clusterE)

clE = dt[,list(mean.prop=mean(hw.logit,na.rm=T),stdev.prop=sd(hw.logit,na.rm=T),n.prop=length(!is.na(proportion.pigmented))),by="two.week.bin.midpoint"]

#create time series object for layeranalyzer
X=layer.data.series(time.points=clE$two.week.bin.midpoint,
  value.points=clE$mean.prop,
  std.dev=clE$stdev.prop, 
  num.meas.per.value=clE$n.prop,name="prop.pigmented")
  
  
#plotting time series
plot(X$time, X$value,type="b",ylim=c(-3.5,5))
for(i in 1:length(X$time))
  lines(c(X$time[i],X$time[i]),
        c(X$value[i]-1.96*X$std.dev[i]/sqrt(X$num.meas.per.value[i]),
          X$value[i]+1.96*X$std.dev[i]/sqrt(X$num.meas.per.value[i])))
          
#explore models with various layers and options to find the best one

models.ml=traverse.standalone.layered(X, max.layers=2, 
   talkative=TRUE, allow.one.feedback.loop=TRUE, 
   just.stationary=FALSE, no.rw=FALSE,    
   time.integrals.possible=FALSE, 
   allow.deterministic.layers=TRUE,
   do.maximum.likelihood = TRUE, maximum.likelihood.numstart = 1000, 
   num.MCMC=1000,spacing=10,burnin=2000, num.temp = 4, prior=p)
          
#best model is model 8, with a random walk layer 2 and an OU tracking layer 1
#fit with MCMC 

X.mod8= layer.series.structure(X, numlayers=2, lin.time=F, prior=p, init.0=TRUE, no.sigma=0,no.pull=TRUE,time.integral=0)       

mod8=layer.analyzer(X.mod8, num.MCMC=1000, burnin=10000,spacing=10,num.temp=6,
   do.model.likelihood = TRUE,maximum.likelihood.numstart=100,
   return.residuals=TRUE)

summary(mod8)

#Description:
#2-layered: Layer 2: RW , Layer 1: OU-like tracking
#
#Coefficients:
#                               Mean    Median Lower 95% Upper 95%
#dt_prop.pigmented_1        7.985321  4.178688  0.583300 40.679426
#sigma_prop.pigmented_1     0.106406  0.075926  0.007596  0.375251
#sigma_prop.pigmented_2     0.234010  0.216432  0.140052  0.424261
#obs_sd_prop.pigmented      0.181102  0.105135  0.007970  0.865200
#init_prop.pigmented_l1_s0 -1.708247 -1.769451 -2.205171 -0.867496
#init_prop.pigmented_l2_s0 -0.558349 -0.666969 -3.011478  2.246714
#
#Model log-likelihood:   -35.681
