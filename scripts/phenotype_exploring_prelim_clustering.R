#This script plots phenotypic data in different ways as a way of exploring different approaches to clustering

require(data.table)
require(ggplot2)

phen<-read.csv("~/Desktop/MASTER_titia_phenotypes.csv")

dt<-data.table(phen)
newd<-dt[,list(count=sum(!is.na(draft.estimate.prop.pigment)),sites=paste(unique(Site_code_revised),collapse="_")),by="HYBAS_ID_5,ddRAD_hydro_ID_5"]

dfn<-as.data.frame(newd)

dft<-subset(dfn,ddRAD_hydro_ID_5)

for(i in 1:dim(dft)[1]){
	p<-ggplot(subset(phen,HYBAS_ID_5==dft$HYBAS_ID_5[i]),aes(x=Julian.date,y=draft.estimate.prop.pigment,col=best.estimate.source))+geom_point()+ggtitle(dft$sites[i])+xlim(c(0,365))+ylim(c(0,1))
	plot(p)
	readline(prompt = "Pause. Press <Enter> to continue...")
}


#> table(phen$best.estimate.source)
#    iNat.AI iNat.Manual      Std.AI  Std.Manual 
#       1552         458         465        3628 



all.dist = rdist.earth(phen[,c(10,11)],miles=FALSE,R=6371) #change R according to #from: https://www.r-bloggers.com/2010/11/great-circle-distance-calculations-in-r/

furthest.within<-c() 
closest.between<-c()


for(i in 1:dim(dft)[1]){

	basin = subset(phen,HYBAS_ID_5==dft$HYBAS_ID_5[i])
	basin.IDs=which(phen$HYBAS_ID_5==dft$HYBAS_ID_5[i])
	other.IDs=which(phen$HYBAS_ID_5!=dft$HYBAS_ID_5[i])
	furthest.within = c(furthest.within,max(all.dist[basin.IDs,basin.IDs]))
	closest.between = c(closest.between,min(all.dist[other.IDs,basin.IDs]))

}



##looking at hydrobasin 3


require(data.table)
require(ggplot2)

phen<-read.csv("~/Desktop/MASTER_titia_phenotypes.csv")

dt<-data.table(phen)
newd<-dt[,list(count=sum(!is.na(draft.estimate.prop.pigment)),sites=paste(unique(Site_code_revised),collapse="_")),by="HYBAS_ID_3,ddRAD_hydro_ID_3"]

dfn<-as.data.frame(newd)

dft<-subset(dfn,ddRAD_hydro_ID_3)

for(i in 1:dim(dft)[1]){
	#p<-ggplot(subset(phen,HYBAS_ID_3==dft$HYBAS_ID_3[i]),aes(x=Julian.date,y=draft.estimate.prop.pigment,col=best.estimate.source))+geom_point()+ggtitle(dft$sites[i])+xlim(c(0,365))+ylim(c(0,1))
	
	#p<-ggplot(subset(phen,HYBAS_ID_3==dft$HYBAS_ID_3[i]),aes(x=Julian.date,y=asin(sqrt(draft.estimate.prop.pigment)),col=best.estimate.source))+geom_point()+ggtitle(dft$sites[i])+xlim(c(0,365))
	#p<-ggplot(subset(phen,HYBAS_ID_3==dft$HYBAS_ID_3[i]),aes(x=Julian.date,y=log(draft.estimate.prop.pigment/(1-draft.estimate.prop.pigment)),col=best.estimate.source))+geom_point()+ggtitle(dft$sites[i])+xlim(c(0,365))
	p<-ggplot(subset(phen,HYBAS_ID_3==dft$HYBAS_ID_3[i]),aes(x=Julian.date,y=asin(sqrt(draft.estimate.prop.pigment))))+geom_point()+geom_smooth(method="loess")+ggtitle(dft$sites[i])+xlim(c(0,365))

	plot(p)
	readline(prompt = "Pause. Press <Enter> to continue...")
}

#playing around with clustering by lon/lat
library(NbClust)
clust.out<-NbClust(phen[,10:11],method="complete",index="all")
save(clust.out,file="clust.out.RData")



dt<-data.table(phen)
newd<-dt[,list(count=sum(!is.na(draft.estimate.prop.pigment)),sites=paste(unique(Site_code_revised),collapse="_")),by="HYBAS_ID_4,ddRAD_hydro_ID_4"]

dfn<-as.data.frame(newd)

dft<-subset(dfn,ddRAD_hydro_ID_4)

for(i in 1:dim(dft)[1]){
	#p<-ggplot(subset(phen,HYBAS_ID_3==dft$HYBAS_ID_3[i]),aes(x=Julian.date,y=draft.estimate.prop.pigment,col=best.estimate.source))+geom_point()+ggtitle(dft$sites[i])+xlim(c(0,365))+ylim(c(0,1))
	
	#p<-ggplot(subset(phen,HYBAS_ID_3==dft$HYBAS_ID_3[i]),aes(x=Julian.date,y=asin(sqrt(draft.estimate.prop.pigment)),col=best.estimate.source))+geom_point()+ggtitle(dft$sites[i])+xlim(c(0,365))
	#p<-ggplot(subset(phen,HYBAS_ID_3==dft$HYBAS_ID_3[i]),aes(x=Julian.date,y=log(draft.estimate.prop.pigment/(1-draft.estimate.prop.pigment)),col=best.estimate.source))+geom_point()+ggtitle(dft$sites[i])+xlim(c(0,365))
	p<-ggplot(subset(phen,HYBAS_ID_4==dft$HYBAS_ID_4[i]),aes(x=Julian.date,y=asin(sqrt(draft.estimate.prop.pigment))))+geom_point()+geom_smooth(method="loess")+ggtitle(dft$sites[i])+xlim(c(0,365))

	plot(p)
	readline(prompt = "Pause. Press <Enter> to continue...")
}

