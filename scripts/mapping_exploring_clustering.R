#SCRIPT used to plot hydrobasins at different levels while exploring different clustering approaches

library(sf)
library(ggplot2)
library(ggmap)
library(elevatr)
#install.packages("elevatr")
library(terra)

phen<-read.csv("data/MASTER_titia_phenotypes_v1.0.csv")

worldmap <- read_sf("data/world_data/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
e <- c(xmin = -130, xmax = -60, ymin = 0, ymax = 45)
worldmap <- st_crop(worldmap, e)

hydrobasins03 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev03_v1c.shp")
hydrobasins03.sub <- hydrobasins03[hydrobasins03$ENDO==0,]
hydrobasins03_geo <- st_geometry(hydrobasins03.sub)
pdf(file="~/Downloads/hydrobasins03.pdf")
plot(hydrobasins03_geo)
dev.off()

hydrobasins04 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev04_v1c.shp")
hydrobasins04.sub <- hydrobasins04[hydrobasins04$ENDO==0,]
hydrobasins04_geo <- st_geometry(hydrobasins04.sub)
pdf(file="~/Downloads/hydrobasins04.pdf")
plot(hydrobasins04_geo)
dev.off()


hydrobasins05 <- sf::st_read("data/hydrosheds/hybas_na_lev01-06_v1c/hybas_na_lev05_v1c.shp")
hydrobasins05.sub <- hydrobasins05[hydrobasins05$ENDO==0,]
hydrobasins05_geo <- st_geometry(hydrobasins05.sub)
pdf(file="~/Downloads/hydrobasins05.pdf")
plot(hydrobasins05_geo)
dev.off()

ggplot() + geom_sf(data = hydrobasins04_geo) #this version fills the countries

require(maps)

NAm_map<-map_data("world",region=c("Canada","USA","Mexico","Belize","Guatemala","Honduras","El Salvador","Nicaragua","Costa Rica","Panama"))

ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins04.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_4,HYBAS_ID%in%phen$HYBAS_ID_4[which(phen$ddRAD_hydro_ID_4)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

pdf("plots/cluster_exploration/hydrobasins04_sample_overlap.pdf")
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins04.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_4,HYBAS_ID%in%phen$HYBAS_ID_4[which(phen$ddRAD_hydro_ID_4)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+geom_sf_text(data = hydrobasins04.sub, aes(label = HYBAS_ID),size=1)+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")
dev.off()

ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins05.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_5,HYBAS_ID%in%phen$HYBAS_ID_5[which(phen$ddRAD_hydro_ID_5)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins05.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_5,HYBAS_ID%in%phen$HYBAS_ID_5[which(phen$ddRAD_hydro_ID_5)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

pdf("plots/cluster_exploration/hydrobasins05_sample_overlap.pdf")
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins05.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_5,HYBAS_ID%in%phen$HYBAS_ID_5[which(phen$ddRAD_hydro_ID_5)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+geom_sf_text(data = hydrobasins05.sub, aes(label = HYBAS_ID),size=1)+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")
dev.off()

ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins03.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_3,HYBAS_ID%in%phen$HYBAS_ID_3[which(phen$ddRAD_hydro_ID_3)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

pdf("plots/cluster_exploration/hydrobasins03_sample_overlap.pdf")
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins03.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_3,HYBAS_ID%in%phen$HYBAS_ID_3[which(phen$ddRAD_hydro_ID_3)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+geom_sf_text(data = hydrobasins03.sub, aes(label = HYBAS_ID),size=1)+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")
dev.off()

# Function written by Jonathan not sure what its for
plot_basin<-function(level=5,basin=7050000010){

if(level==5){
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins05.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_5,HYBAS_ID==basin)))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

}

if(level==4){
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins04.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_4,HYBAS_ID==basin)))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

}


}


p<-ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

col_pal <- rainbow(13)

p+geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID==7050000010),aes(fill=col_pal[1])) +
	geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050004220, 7050005420)),aes(fill=col_pal[2]))+
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050047850, 7050736840, 7050047980, 7050048620)),aes(fill=col_pal[3])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050668250, 7050741770, 7050682360)),aes(fill=col_pal[4])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050822250, 7050823030)),aes(fill=col_pal[5])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050050000, 7050050930, 7050849600)),aes(fill=col_pal[6])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050050640)),aes(fill=col_pal[7])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050054660, 7050054670)),aes(fill=col_pal[8])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050043860, 7050043870, 7050045480)),aes(fill=col_pal[9])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050045900)),aes(fill=col_pal[10])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050041410, 7050453940, 7050494360)),aes(fill=col_pal[11])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050609780)),aes(fill=col_pal[12])) +
  geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050511950, 7050531210)),aes(fill=col_pal[13]))