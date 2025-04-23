#this script labels the peak and trough mean and standard error for each tip
#current data transformation is a logit, subtracting 0.01 from each value so 1 values don't hit Inf
setwd("~/Thesis_H_titia_seasonal_polyphenism_evolution")
phen<-read.csv("data/MASTER_titia_phenotypes_v1.0.csv")
peak.trough.season <- read.csv("data/peak_trough_season_dates.csv")


#as a result of the normalisation of non standard manual measurements, a few measurements are very small
#consequently, removing all values < 0.05 (the smallest hindwingspot in a standard manually measured wing)
#this removes 35 data points

phen<-subset(phen,draft.estimate.prop.pigment>0.05)

lt<-LETTERS[1:8]
tips.a.to.h<-c("7050000010-NB0103",  "7050002710-ZANAa05","7050054670-BALB05", "7050849600-TXRSa04",  "7050045480-ACSFc04","7050041410-RCJRa01","7050822250-CYa04" , "7050048620-CVa06" ) 

time.vec = vector(length=dim(phen)[1])

for(i in 1:dim(phen)[1]){
	
	if(is.na(phen[i,]$CLUSTER)){
	
	time.vec[i] = "NA"
	
	} else {
	
		if(phen[i,]$Julian.date>=287 || phen[i,]$Julian.date<=91){
	
		time.vec[i] = "trough"
		
		} else {
		
			cl = which(lt==phen[i,]$CLUSTER)
		
			if(phen[i,]$Julian.date>=peak.trough.season[cl,2] && phen[i,]$Julian.date<=peak.trough.season[cl,3]){
	
			time.vec[i] = "peak"
			
			} else {
			
			time.vec[i] = "NA"

			
			}
		
		}
		
	} 

}
	

phen$time.period = time.vec	
	
write.csv(phen,file="~/Desktop/output.csv")