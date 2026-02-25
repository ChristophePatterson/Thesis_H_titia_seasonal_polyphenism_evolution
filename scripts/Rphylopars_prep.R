#this script calculates the peak and trough mean and standard error for each tip
#current data transformation is a logit, subtracting 0.01 from each value so 1 values don't hit Inf
setwd("~/Thesis_H_titia_seasonal_polyphenism_evolution")
phen<-read.csv("data/MASTER_titia_phenotypes_v2.0_processed.csv")
peak.trough.season <- read.csv("data/peak_trough_season_dates.csv")



lt<-LETTERS[1:8]
tips.a.to.h<-c("7050000010-NB0103",  "7050002710-ZANAa05","7050054670-BALB05", "7050849600-TXRSa04",  "7050045480-ACSFc04","7050041410-RCJRa01","7050822250-CYa04" , "7050048620-CVa06" ) 


for(i in 1:length(lt)){
	lt.sub<-subset(phen,CLUSTER==lt[i])
	lt.peak.sub<-subset(lt.sub,Julian.date>=peak.trough.season[i,2],Julian.date<=peak.trough.season[i,3])
	hw.trans.peak = log( (lt.peak.sub$proportion.pigmented-0.01) / (1-((lt.peak.sub$proportion.pigmented-0.01))))
	
	if(i ==1){
	
		out.mat = cbind(lt[i],hw.trans.peak,"peak")
	
	} else {
	
		out.mat = rbind(out.mat, cbind(lt[i],hw.trans.peak,"peak"))
	}
	
	
	lt.trough.sub<-subset(lt.sub,Julian.date>=287 | Julian.date<=91)
	hw.trans.trough = log( (lt.trough.sub$proportion.pigmented-0.01) / (1-((lt.trough.sub$proportion.pigmented-0.01))))
		out.mat = rbind(out.mat, cbind(lt[i],hw.trans.trough,"trough"))
}




write.csv(out.mat,file="rphylopars_data.csv")