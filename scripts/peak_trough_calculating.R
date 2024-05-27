## This script explores different ways of modelling the seasonal polyphenism (lines 3-120)

## It also plots the phenotypic data in peak and trough as calculated using double logistic curve fitting approach (lines 121 and beyond)

## For further details see Planning phylogenetic analyses of titia polyphenism.docx


require(ggplot2)
library(data.table)
require(gridExtra)
phen<-read.csv("data/MASTER_titia_phenotypes_v1.0.csv")


ggplot(data=subset(phen,CLUSTER!="NA"),aes(y= asin(sqrt(draft.estimate.prop.pigment)),x=Julian.date,group=CLUSTER,col=CLUSTER))+geom_point()+ stat_smooth(method = "lm", formula = y ~ x + I(x^2))+facet_wrap(~CLUSTER)


phen$Julian.date2<-phen$Julian.date^2
m.l.global<-lm(asin(sqrt(draft.estimate.prop.pigment))~Julian.date+Julian.date2,data=subset(phen, draft.estimate.prop.pigment<=1))

#peak = -b/2a in ax2 + bx + c, so peak = 195.35 in this global model


m.l.A<-lm(asin(sqrt(draft.estimate.prop.pigment))~Julian.date+Julian.date2,data=subset(phen, draft.estimate.prop.pigment<=1 & CLUSTER=="A"))


ggplot(data=subset(phen,CLUSTER!="NA" & draft.estimate.prop.pigment<=1),aes(y=draft.estimate.prop.pigment,x=Julian.date,group=CLUSTER,col=CLUSTER))+geom_point()+ stat_smooth(method = "glm", formula = y ~ x + I(x^2), method.args = list(family='binomial'))+facet_wrap(~CLUSTER)


dt<-data.table(phen)
out<-dt[,list(day.mean=mean(draft.estimate.prop.pigment)),by="Julian.date,CLUSTER"]

out<-out[order(Julian.date),]
out<-out[order(CLUSTER),]
frollmean(out[,day.mean],7,na.rm=TRUE,hasNA=TRUE)->rolling_mean
out<-as.data.frame(out)
out$rolling_mean<-rolling_mean
ggplot(data=subset(out,CLUSTER!="<NA>"),aes(y=rolling_mean,x=Julian.date,group=CLUSTER,col=CLUSTER))+geom_point()+ geom_smooth()+facet_wrap(~CLUSTER)





ACI.curve <- function(P,data){ # inputs are parameters (P) and dataset
  
  # First, scale parameter values, this makes it easier for optim, as you feed it
  # parameters on the same scale
  (d.star<-P[1]/10)
  (sigma<-P[2]/20)
  (alpha1<-P[3]/1200)
  (alpha2<-P[4]/200000)
  (phi<-P[5]/10000)
  (z<-P[6]/250)
  
  if(alpha2>=alpha1){alpha2=alpha1-0.0001}
  
  lambda<-alpha2+(alpha1-alpha2)*exp(-(abs((c(data$date)-d.star)/sigma)^z)) # This is equation for function
  
  if (length(which(lambda==1))>0 | length(which(lambda==0))>0) LL.ind<-rep(-1000000,length(data$date)) else {
    
    a<-lambda/phi
    b<-(1-lambda)/phi
    
    LL.ind<-lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(data$scaledACI+0.0001)+(b-1)*log(1-data$scaledACI+0.0001) # This calculates likelihood
  }
  
  return(sum(LL.ind)) # Sums likelihood across each value
}




####dt<-data.table(inat)
#data<-dt[,list(mean.hw=mean(perched.hw2,na.rm=T)),by="julian.date"]
data<-phen[,c(10,40)]
data<-data.frame(data)
    # Remove any NAs
data<-data[which(is.na(data[,2])==F),]
data<-data[which(is.na(data[,1])==F),]
colnames(data)[1]<-"date"
data$scaledACI<-(data[,2]-min(data[,2]))/(max(data[,2])-min(data[,2]))
    # Sequence of dates
date.seq<-seq(min(data[,1]),max(data[,1]),1)
#data<-data[order(data[,1]),]
    
    ## Fit curve ##
    initial<-c(1300,1000,1000,1000,1000,1000) # Initial parameter values, all on roughly same scale
    low<-c(1,200,100,1,1,250) # Min allowed parameter values
    high<-c(2000,2000,1200,10000,10000,5000) # Max allowed parameter values
    result<-optim(par=initial,fn=ACI.curve,data=data,method="L-BFGS-B",lower=low,upper=high,control=list(trace=1,fnscale=-1,maxit=1000)) # Run optim
    P.new<-c(result$par[1]/10,result$par[2]/20,result$par[3]/1200,result$par[4]/200000,result$par[5]/10000,result$par[6]/250) # Rescale parameters
    lambda<-P.new[4]+(P.new[3]-P.new[4])*exp(-(abs((date.seq-P.new[1])/(P.new[2]))^(P.new[6]))) # Fit function with fitted parameter values
    
    ### Differentiating ###
    first.deriv<-(P.new[6]*exp(-((P.new[1]-date.seq)^2/P.new[2]^2)^(P.new[6]/2)) * ((P.new[1]-date.seq)^2/P.new[2]^2)^(P.new[6]/2)) / (P.new[1]-date.seq)
    second.deriv<-( P.new[6]*exp(-((P.new[1]-date.seq)^2/P.new[2]^2)^(P.new[6]/2)) * ((P.new[1]-date.seq)^2/P.new[2]^2)^((P.new[6]-2)/2) * (P.new[6]*(((P.new[1]-date.seq)^2/P.new[2]^2)^(P.new[6]/2) -1) +1) ) / (P.new[1]^2)
    deriv.data<-as.data.frame(cbind(first=first.deriv,second=second.deriv,date=date.seq))
    #plot(deriv.data$first~deriv.data$date,type="l") # Plot first derivative if you're interested
    #plot(deriv.data$second~deriv.data$date,type="l") # Plot second derivative if you're interested
    # Date of first max 2nd deritave
    max.second1<-which(deriv.data[deriv.data$date<=P.new[1],'second'] %in% max(deriv.data[deriv.data$date<=P.new[1],'second']))
    (date.max.second1<-deriv.data[deriv.data$date<=P.new[1],'date'][max.second1])
    # Date of second max 2nd deritave
    max.second2<-which(deriv.data[deriv.data$date>=P.new[1],'second'] %in% max(deriv.data[deriv.data$date>=P.new[1],'second']))
    (date.max.second2<-deriv.data[deriv.data$date>=P.new[1],'date'][max.second2])
    # Date of first max 2nd deritave
    max.first<-max(deriv.data[deriv.data$date<=P.new[1],'first'])
    max.first1<-which(deriv.data[deriv.data$date<=P.new[1],'first'] %in% max.first)
    (date.max.first<-deriv.data[deriv.data$date<=P.new[1],'date'][max.first1])
    (lambda.max.first<-P.new[4]+(P.new[3]-P.new[4])*exp(-(abs((date.max.first-P.new[1])/(P.new[2]))^(P.new[6]))))
    # Date of second max 2nd deritave
    min.first<-min(deriv.data[deriv.data$date>=P.new[1],'first'])
    min.first2<-which(deriv.data[deriv.data$date>=P.new[1],'first'] %in% min.first)
    (date.min.first<-deriv.data[deriv.data$date>=P.new[1],'date'][min.first2])
    (lambda.min.first<-P.new[4]+(P.new[3]-P.new[4])*exp(-(abs((date.min.first-P.new[1])/(P.new[2]))^(P.new[6]))))
    
    ### Plotting ###
    with(data,plot(scaledACI~date,type="n",ylab="perched hindwing",xlab="date",xlim=c(0,365),main=names(data)[2])) # Plot
    ## Filling area under curve
    lambda.new<-P.new[4]+(P.new[3]-P.new[4])*exp(-(abs((seq(date.max.second1,date.max.second2,0.1)-P.new[1])/(P.new[2]))^(P.new[6])))
    lambda.one<-P.new[4]+(P.new[3]-P.new[4])*exp(-(abs((date.max.second1-P.new[1])/(P.new[2]))^(P.new[6])))
    lambda.growing<-P.new[4]+(P.new[3]-P.new[4])*exp(-(abs((seq(round(date.max.second1,0),round(date.max.second2,0),1)-P.new[1])/(P.new[2]))^(P.new[6])))
    x<-c(seq(date.max.second1,date.max.second2,0.1),rev(seq(date.max.second1,date.max.second2,0.1)))
    y<-c(lambda.new,rep(0,length(lambda.new)))
    pol <- cbind(x, y)
    polygon(x,y,col="light grey",border=NA)
    lines(lambda~date.seq,col=1,lwd=2)
    ## Add points
    with(data,points(scaledACI~date,col="dark grey"))
    ## Annotating
    # duration
    segments(date.max.second1,0,date.max.second1,.9,lty=2,col=2); text(date.max.second1,.93,"start",cex=0.9)
    segments(date.max.second2,0,date.max.second2,.9,lty=2,col=2); text(date.max.second2,.93,"end",cex=0.9)
    text(P.new[1],.85,"duration",cex=0.9)
    arrows(P.new[1]-15,.85,date.max.second1,.85,length=.1,col=2)
    arrows(P.new[1]+15,.85,date.max.second2,.85,length=.1,col=2)
    # sum
    text(P.new[1],.4,"sum",cex=0.9)
    # peak
    arrows(35,max(lambda),P.new[1],max(lambda),length=.1,col=2)
    text(25,max(lambda),"peak",cex=0.9)
    # rate
    arrows(35,lambda.max.first,date.max.first,lambda.max.first,length=.1,col=2)
    text(22,lambda.max.first,"max rate",cex=.9)
    arrows(225,lambda.min.first,date.min.first,lambda.min.first,length=.1,col=2)
    text(238,lambda.max.first,"min rate",cex=.9)




###PLOTTING 
phen_mean<-read.csv("data/peak_trough_v1.0.csv")

clA<-subset(phen,CLUSTER=="A")
clA.peak.start = 132
clA.peak.end = 160
clA$col<-ifelse(clA$Julian.date >= 287 | clA$Julian.date <= 91,"trough",ifelse(clA$Julian.date>=clA.peak.start & clA$Julian.date<=clA.peak.end,"peak","neutral"))

clA_plot_A<-ggplot(data=clA,aes(x=Julian.date,y= asin(sqrt(draft.estimate.prop.pigment)),col=col))+geom_point()+scale_colour_manual(values=c("grey90","firebrick4","steelblue1"))+xlim(c(0,365))+ylim(c(0,1.6))+theme_bw()+ theme(legend.position = "none")+geom_segment(aes(x=0, xend=91,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=287, xend=365,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=clA.peak.start, xend=clA.peak.end,y=1.6,yend=1.6),col="firebrick4",lwd=1.5)

clA_plot_B<-ggplot(data=subset(phen_mean,CLUSTER=="A"),aes(x=time,y=mean,col=time))+geom_point(cex=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0,lwd=2)+scale_colour_manual(values=c("firebrick4","steelblue1"))+theme_bw()+ylim(c(0,1.6))+theme(legend.position="none")+ylab("")

clB<-subset(phen,CLUSTER=="B")
clB.peak.start = 186
clB.peak.end = 214
clB$col<-ifelse(clB$Julian.date >= 287 | clB$Julian.date <= 91,"trough",ifelse(clB$Julian.date>=clB.peak.start & clB$Julian.date<=clB.peak.end,"peak","neutral"))

clB_plot_A<-ggplot(data=clB,aes(x=Julian.date,y= asin(sqrt(draft.estimate.prop.pigment)),col=col))+geom_point()+scale_colour_manual(values=c("grey90","firebrick4","steelblue1"))+xlim(c(0,365))+ylim(c(0,1.6))+theme_bw()+ theme(legend.position = "none")+geom_segment(aes(x=0, xend=91,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=287, xend=365,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=clB.peak.start, xend=clB.peak.end,y=1.6,yend=1.6),col="firebrick4",lwd=1.5)
clB_plot_B<-ggplot(data=subset(phen_mean,CLUSTER=="B"),aes(x=time,y=mean,col=time))+geom_point(cex=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0,lwd=2)+scale_colour_manual(values=c("firebrick4","steelblue1"))+theme_bw()+ylim(c(0,1.6))+theme(legend.position="none")+ylab("")

clC<-subset(phen,CLUSTER=="C")
clC.peak.start = 171
clC.peak.end = 199
clC$col<-ifelse(clC$Julian.date >= 287 | clC$Julian.date <= 91,"trough",ifelse(clC$Julian.date>=clC.peak.start & clC$Julian.date<=clC.peak.end,"peak","neutral"))

clC_plot_A<-ggplot(data=clC,aes(x=Julian.date,y= asin(sqrt(draft.estimate.prop.pigment)),col=col))+geom_point()+scale_colour_manual(values=c("grey90","firebrick4","steelblue1"))+xlim(c(0,365))+ylim(c(0,1.6))+theme_bw()+ theme(legend.position = "none")+geom_segment(aes(x=0, xend=91,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=287, xend=365,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=clC.peak.start, xend=clC.peak.end,y=1.6,yend=1.6),col="firebrick4",lwd=1.5)
clC_plot_B<-ggplot(data=subset(phen_mean,CLUSTER=="C"),aes(x=time,y=mean,col=time))+geom_point(cex=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0,lwd=2)+scale_colour_manual(values=c("firebrick4","steelblue1"))+theme_bw()+ylim(c(0,1.6))+theme(legend.position="none")+ylab("")

clD<-subset(phen,CLUSTER=="D")
clD.peak.start = 186
clD.peak.end = 214
clD$col<-ifelse(clD$Julian.date >= 287 | clD$Julian.date <= 91,"trough",ifelse(clD$Julian.date>=clD.peak.start & clD$Julian.date<=clD.peak.end,"peak","neutral"))

clD_plot_A<-ggplot(data=clD,aes(x=Julian.date,y= asin(sqrt(draft.estimate.prop.pigment)),col=col))+geom_point()+scale_colour_manual(values=c("grey90","firebrick4","steelblue1"))+xlim(c(0,365))+ylim(c(0,1.6))+theme_bw()+ theme(legend.position = "none")+geom_segment(aes(x=0, xend=91,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=287, xend=365,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=clD.peak.start, xend=clD.peak.end,y=1.6,yend=1.6),col="firebrick4",lwd=1.5)
clD_plot_B<-ggplot(data=subset(phen_mean,CLUSTER=="D"),aes(x=time,y=mean,col=time))+geom_point(cex=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0,lwd=2)+scale_colour_manual(values=c("firebrick4","steelblue1"))+theme_bw()+ylim(c(0,1.6))+theme(legend.position="none")+ylab("")

clE<-subset(phen,CLUSTER=="E")
clE.peak.start = 174
clE.peak.end = 202
clE$col<-ifelse(clE$Julian.date >= 287 | clE$Julian.date <= 91,"trough",ifelse(clE$Julian.date>=clE.peak.start & clE$Julian.date<=clE.peak.end,"peak","neutral"))

clE_plot_A<-ggplot(data=clE,aes(x=Julian.date,y= asin(sqrt(draft.estimate.prop.pigment)),col=col))+geom_point()+scale_colour_manual(values=c("grey90","firebrick4","steelblue1"))+xlim(c(0,365))+ylim(c(0,1.6))+theme_bw()+ theme(legend.position = "none")+geom_segment(aes(x=0, xend=91,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=287, xend=365,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=clE.peak.start, xend=clE.peak.end,y=1.6,yend=1.6),col="firebrick4",lwd=1.5)
clE_plot_B<-ggplot(data=subset(phen_mean,CLUSTER=="E"),aes(x=time,y=mean,col=time))+geom_point(cex=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0,lwd=2)+scale_colour_manual(values=c("firebrick4","steelblue1"))+theme_bw()+ylim(c(0,1.6))+theme(legend.position="none")+ylab("")

clF<-subset(phen,CLUSTER=="F")
clF.peak.start = 170
clF.peak.end = 198
clF$col<-ifelse(clF$Julian.date >= 287 | clF$Julian.date <= 91,"trough",ifelse(clF$Julian.date>=clF.peak.start & clF$Julian.date<=clF.peak.end,"peak","neutral"))

clF_plot_A<-ggplot(data=clF,aes(x=Julian.date,y= asin(sqrt(draft.estimate.prop.pigment)),col=col))+geom_point()+scale_colour_manual(values=c("grey90","firebrick4","steelblue1"))+xlim(c(0,365))+ylim(c(0,1.6))+theme_bw()+ theme(legend.position = "none")+geom_segment(aes(x=0, xend=91,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=287, xend=365,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=clF.peak.start, xend=clF.peak.end,y=1.6,yend=1.6),col="firebrick4",lwd=1.5)
clF_plot_B<-ggplot(data=subset(phen_mean,CLUSTER=="F"),aes(x=time,y=mean,col=time))+geom_point(cex=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0,lwd=2)+scale_colour_manual(values=c("firebrick4","steelblue1"))+theme_bw()+ylim(c(0,1.6))+theme(legend.position="none")+ylab("")

clG<-subset(phen,CLUSTER=="G")
clG.peak.start = 186
clG.peak.end = 214
clG$col<-ifelse(clG$Julian.date >= 287 | clG$Julian.date <= 91,"trough",ifelse(clG$Julian.date>=clG.peak.start & clG$Julian.date<=clG.peak.end,"peak","neutral"))

clG_plot_A<-ggplot(data=clG,aes(x=Julian.date,y= asin(sqrt(draft.estimate.prop.pigment)),col=col))+geom_point()+scale_colour_manual(values=c("grey90","firebrick4","steelblue1"))+xlim(c(0,365))+ylim(c(0,1.6))+theme_bw()+ theme(legend.position = "none")+geom_segment(aes(x=0, xend=91,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=287, xend=365,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=clG.peak.start, xend=clG.peak.end,y=1.6,yend=1.6),col="firebrick4",lwd=1.5)
clG_plot_B<-ggplot(data=subset(phen_mean,CLUSTER=="G"),aes(x=time,y=mean,col=time))+geom_point(cex=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0,lwd=2)+scale_colour_manual(values=c("firebrick4","steelblue1"))+theme_bw()+ylim(c(0,1.6))+theme(legend.position="none")+ylab("")

clH<-subset(phen,CLUSTER=="H")
clH.peak.start = 186
clH.peak.end = 214
clH$col<-ifelse(clH$Julian.date >= 287 | clH$Julian.date <= 91,"trough",ifelse(clH$Julian.date>=clH.peak.start & clH$Julian.date<=clH.peak.end,"peak","neutral"))

clH_plot_A<-ggplot(data=clH,aes(x=Julian.date,y= asin(sqrt(draft.estimate.prop.pigment)),col=col))+geom_point()+scale_colour_manual(values=c("grey90","firebrick4","steelblue1"))+xlim(c(0,365))+ylim(c(0,1.6))+theme_bw()+ theme(legend.position = "none")+geom_segment(aes(x=0, xend=91,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=287, xend=365,y=0,yend=0),col="steelblue1",lwd=1.5)+geom_segment(aes(x=clH.peak.start, xend=clH.peak.end,y=1.6,yend=1.6),col="firebrick4",lwd=1.5)
clH_plot_B<-ggplot(data=subset(phen_mean,CLUSTER=="H"),aes(x=time,y=mean,col=time))+geom_point(cex=2)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0,lwd=2)+scale_colour_manual(values=c("firebrick4","steelblue1"))+theme_bw()+ylim(c(0,1.6))+theme(legend.position="none")+ylab("")

grid.arrange(clA_plot_A,clA_plot_B, clB_plot_A,clB_plot_B,clC_plot_A,clC_plot_B, clD_plot_A,clD_plot_B,ncol=2,widths=c(2.75,1))
grid.arrange(clE_plot_A,clE_plot_B, clF_plot_A,clF_plot_B,clG_plot_A,clG_plot_B, clH_plot_A,clH_plot_B,ncol=2,widths=c(2.75,1))


peak.trough.season <- rbind(c("A", clA.peak.start, clA.peak.end),
      c("B", clB.peak.start, clB.peak.end),
      c("C", clC.peak.start, clC.peak.end),
      c("D", clD.peak.start, clD.peak.end),
      c("E", clE.peak.start, clE.peak.end),
      c("F", clF.peak.start, clF.peak.end),
      c("G", clG.peak.start, clG.peak.end),
      c("H", clH.peak.start, clH.peak.end))

peak.trough.season <- data.frame(peak.trough.season)
colnames(peak.trough.season) <- c("cluster", "peak.start", "peak.end")
write.table(peak.trough.season, file = "data/peak_trough_season_dates.csv", sep = ",")

