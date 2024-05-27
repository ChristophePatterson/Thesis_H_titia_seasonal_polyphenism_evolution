#SCRIPT used to plot hydrobasins at different levels while exploring different clustering approaches

library(sf)
library(ggplot2)
library(ggmap)
library(elevatr)
#install.packages("elevatr")
library(terra)
worldmap <- read_sf("~/Downloads/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
e <- c(xmin = -130, xmax = -60, ymin = 0, ymax = 45)
worldmap <- st_crop(worldmap, e)



hydrobasins03 <- sf::st_read("~/Downloads/hybas_na_lev01-12_v1c/hybas_na_lev03_v1c.shp")
hydrobasins03.sub <- hydrobasins03[hydrobasins03$ENDO==0,]
hydrobasins03_geo <- st_geometry(hydrobasins03.sub)
pdf(file="~/Downloads/hydrobasins03.pdf")
plot(hydrobasins03_geo)
dev.off()

hydrobasins04 <- sf::st_read("~/Downloads/hybas_na_lev01-12_v1c/hybas_na_lev04_v1c.shp")
hydrobasins04.sub <- hydrobasins04[hydrobasins04$ENDO==0,]
hydrobasins04_geo <- st_geometry(hydrobasins04.sub)
pdf(file="~/Downloads/hydrobasins04.pdf")
plot(hydrobasins04_geo)
dev.off()


hydrobasins05 <- sf::st_read("~/Downloads/hybas_na_lev01-12_v1c/hybas_na_lev05_v1c.shp")
hydrobasins05.sub <- hydrobasins05[hydrobasins05$ENDO==0,]
hydrobasins05_geo <- st_geometry(hydrobasins05.sub)
pdf(file="~/Downloads/hydrobasins05.pdf")
plot(hydrobasins05_geo)
dev.off()

ggplot() + geom_polygon(data = hydrobasins04_geo, aes(x=long, y = lat)) #this version fills the countries

require(maps)

NAm_map<-map_data("world",region=c("Canada","USA","Mexico","Belize","Guatemala","Honduras","El Salvador","Nicaragua","Costa Rica","Panama"))

ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins04.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_4,HYBAS_ID%in%phen$HYBAS_ID_4[which(phen$ddRAD_hydro_ID_4)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

pdf()
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins04.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_4,HYBAS_ID%in%phen$HYBAS_ID_4[which(phen$ddRAD_hydro_ID_4)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+geom_sf_text(data = hydrobasins04.sub, aes(label = HYBAS_ID),size=1)+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")
dev.off()


ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins05.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_5,HYBAS_ID%in%phen$HYBAS_ID_5[which(phen$ddRAD_hydro_ID_5)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins05.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_5,HYBAS_ID%in%phen$HYBAS_ID_5[which(phen$ddRAD_hydro_ID_5)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

pdf()
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins05.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_5,HYBAS_ID%in%phen$HYBAS_ID_5[which(phen$ddRAD_hydro_ID_5)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+geom_sf_text(data = hydrobasins05.sub, aes(label = HYBAS_ID),size=1)+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")
dev.off()

ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins03.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_3,HYBAS_ID%in%phen$HYBAS_ID_3[which(phen$ddRAD_hydro_ID_3)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

pdf()
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins03.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_3,HYBAS_ID%in%phen$HYBAS_ID_3[which(phen$ddRAD_hydro_ID_3)])))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+geom_sf_text(data = hydrobasins03.sub, aes(label = HYBAS_ID),size=1)+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")
dev.off()

plot_basin<-function(level=5,basin=7050000010){

if(level==5){
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins05.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_5,HYBAS_ID==basin)))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

}

if(level==4){
ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+geom_sf(data=hydrobasins04.sub,aes(fill=interaction(HYBAS_ID%in%phen$HYBAS_ID_4,HYBAS_ID==basin)))+scale_fill_manual(values=c("NA","grey90","blue"))+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

}


}


p<-ggplot() + geom_polygon(data = NAm_map, aes(x=long, y = lat,group=group),fill=NA,col="black",lwd=0.25)+xlim(c(-130,-60))+ylim(c(7,45))+theme_bw()+ theme(legend.position = "none")

p+geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID==7050000010),aes(fill=rainbow(8)[2])) +
	geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050002700, 7050002710)),aes(fill=rainbow(8)[1]))+
	geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050007380, 7050007000, 7050004670, 7050004210, 7050006430, 7050006420, 7050005420, 7050004220)),aes(fill=rainbow(8)[1]))+
	geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050054670, 7050054660, 7050054310)),aes(fill=rainbow(8)[3]))+
	geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050050930, 7050050650, 7050050640, 7050050000)),aes(fill=rainbow(8)[4]))+
	geom_sf(data=subset(hydrobasins04.sub,HYBAS_ID%in%c(7040050920)),aes(fill=rainbow(8)[4]))+
	geom_sf(data=subset(hydrobasins04.sub,HYBAS_ID%in%c(7040049990, 7040049280, 7040049270)),aes(fill=rainbow(8)[5]))+
	geom_sf(data=subset(hydrobasins04.sub,HYBAS_ID%in%c(7040048410, 7040048420, 7040048300, 7040047850, 7040047970,7040047980, 7040047070, 7040046750, 7040046600, 7040741760, 7040686450, 7040686540, 7040048310,7040741770, 7040047060,7040047860)),aes(fill=rainbow(8)[6]))+
	geom_sf(data=subset(hydrobasins05.sub,HYBAS_ID%in%c(7050045480, 7050043860, 7050043870)),aes(fill=rainbow(8)[7]))+
	geom_sf(data=subset(hydrobasins04.sub,HYBAS_ID%in%c(7040612640,7040569650,7040392650,7040402430,7040041410, 7040041400)),aes(fill=rainbow(8)[8]))

