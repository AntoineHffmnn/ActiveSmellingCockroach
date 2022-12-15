###########################################################################################################################
#                               Active smelling in the American cockroach: Fig. 3 & Fig. S1B                              #
###########################################################################################################################

# Author: Antoine Hoffmann
# R version 4.2.0



#--------------------------------------------------------------------------------------------------------------------------
# Settings
#--------------------------------------------------------------------------------------------------------------------------


# Folder structure: -Analysis folder
#                     -Data (store data files)
#                     -Main (store .R script files)
#                       -Plots (save plots)





### Packages and functions ###
# run this before any analysis to make sure all necessary packages and custom functions are loaded


require("dplyr")
require("grDevices")
require("viridis")
require("brms")
require("cmdstanr")
require("bayesplot")
require("ggpubr")
require("reshape2")



#Colors
color_scheme_set("purple")

dodgerblue3<-as.vector(col2rgb(col =c("dodgerblue3"),alpha=F ))
firebrick3<-as.vector(col2rgb(col =c("firebrick3"),alpha=F ))


# Function for the opposite of "%in%".
# Returns TRUE if the tested element IS NOT present in the comparison vector
"%!in%" <- function(x,y){
  !("%in%"(x,y))
}






#--------------------------------------------------------------------------------------------------------------------------
# Center odor: Responsive/Non-responsive antennae metrics
#--------------------------------------------------------------------------------------------------------------------------


#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)





########################################
### Metrics for pooled center stream ###
########################################

### data
data<-AllData

distr.metrics.list<-list()
for(i in 1:3){
  
  j<-c("Colony","Butanol","Linalool")[i]
  
  dat<-dplyr::filter(data,
                     Dataset=="3D",
                     Protocol %in% c("ProtocolE","ProtocolB","ProtocolGb"),
                     Odor==j,
                     Count>=-4 & Count<=2)
  dat<-dat[,c("Test","Count","TestPhase","Lant.X","Lant.Y","Lant.Z","Rant.X","Rant.Y","Rant.Z","head.X","head.Y","head.Z")]
  colnames(dat)[4:9]<-c("Ant.X","Ant.Y","Ant.Z","Ant.X","Ant.Y","Ant.Z")
  #restructuring to have an antenna variable
  dat<-rbind(data.frame(dat[,c(1:6,10:12)],"Ant"=rep("L",nrow(dat))),
             data.frame(dat[,c(1:3,7:12)],"Ant"=rep("R",nrow(dat)))) 
  
  dat<-dplyr::filter(dat,Test %!in% c("T618","T662","T786","T787","T944","T945",
                                      "T843","T946"))  # these barely move at all and heavily influences the density plots
  
  dat$Ant.Z[dat$Ant.Z<0]<-0 #change everything below zero to zero, because 0 is supposed to be the ground
  
  
  # Normalization of coordinates for each trial & distribution metrics
  datlist<-split(dat,f=dat$Test)  #separate the df into a list by trial id
  datlist<-lapply(datlist,
                  FUN=function(df){
                    
                    #centering coordinates around the head
                    df$Ant.X<-df$Ant.X-df$head.X
                    df$Ant.Y<-df$Ant.Y-df$head.Y
                    
                    #max.X before norm
                    df$max.X<-max(abs(df$Ant.X),na.rm=T)
                    
                    # normalizing by maximum coordinates within trial
                    df$Ant.X<-df$Ant.X/max(abs(df$Ant.X),na.rm=T)
                    df$Ant.Y<-df$Ant.Y/max(abs(df$Ant.Y),na.rm=T)
                    df$Ant.Z<-df$Ant.Z/max(abs(df$Ant.Z),na.rm=T)
                    
                    
                    return(df)
                    
                  })
  dat<-dplyr::bind_rows(datlist) #put df back together
  dat<-dat[,-c(7:9)] #remove head coordinates
  
  
  
  
  # For each trial get distribution metrics for pooled L&R 
  datlist<-split(dat,f=dat$Test)  
  datlist<-lapply(datlist,
                  FUN=function(df){
                    
                    #left
                    
                    # change during PS as control comparison
                    d<-dplyr::filter(df,TestPhase=="PS",
                                     Count<(-2), #ps1
                                     Ant=="L")
                    lps1width<-max(d$Ant.X,na.rm=T)-min(d$Ant.X,na.rm=T)
                    lps1dist<-mean(abs(d$Ant.X),na.rm=T)
                    
                    d<-dplyr::filter(df,TestPhase=="PS",
                                     Count>=-2, #ps2
                                     Ant=="L")
                    lps2width<-max(d$Ant.X,na.rm=T)-min(d$Ant.X,na.rm=T)
                    lps2dist<-mean(abs(d$Ant.X),na.rm=T)
                    
                    lctrlwidthdiff<-lps2width-lps1width
                    lctrldistdiff<-lps2dist-lps1dist
                    
                    
                    # S change
                    ldps<-dplyr::filter(df,TestPhase=="PS",
                                        Count>=-2,
                                        Ant=="L")
                    lpswidth<-max(ldps$Ant.X,na.rm=T)-min(ldps$Ant.X,na.rm=T)
                    lpsdist<-mean(abs(ldps$Ant.X),na.rm=T)
                    
                    
                    lds<-dplyr::filter(df,TestPhase=="S",Ant=="L")
                    lswidth<-max(lds$Ant.X,na.rm=T)-min(lds$Ant.X,na.rm=T)
                    lsdist<-mean(abs(lds$Ant.X),na.rm=T)

                    lwidthdiff.full<-lswidth-lpswidth
                    ldistdiff<-lsdist-lpsdist
                    
                    
                    
                    
                    #right
                    
                    # change during PS as control comparison
                    d<-dplyr::filter(df,TestPhase=="PS",
                                     Count<(-2), #ps1
                                     Ant=="R")
                    rps1width<-max(d$Ant.X,na.rm=T)-min(d$Ant.X,na.rm=T)
                    rps1dist<-mean(abs(d$Ant.X),na.rm=T)
                    
                    d<-dplyr::filter(df,TestPhase=="PS",
                                     Count>=-2, #ps2
                                     Ant=="R")
                    rps2width<-max(d$Ant.X,na.rm=T)-min(d$Ant.X,na.rm=T)
                    rps2dist<-mean(abs(d$Ant.X),na.rm=T)
                    
                    rctrlwidthdiff<-rps2width-rps1width
                    rctrldistdiff<-rps2dist-rps1dist
                    
                    
                    # S change
                    rdps<-dplyr::filter(df,TestPhase=="PS",
                                        Count>=-2,
                                        Ant=="R")
                    rpswidth<-max(rdps$Ant.X,na.rm=T)-min(rdps$Ant.X,na.rm=T)
                    rpsdist<-mean(abs(rdps$Ant.X),na.rm=T)
                    
                    
                    rds<-dplyr::filter(df,TestPhase=="S",Ant=="R")
                    rswidth<-max(rds$Ant.X,na.rm=T)-min(rds$Ant.X,na.rm=T)
                    rsdist<-mean(abs(rds$Ant.X),na.rm=T)
                    
                    
                    
                    rwidthdiff.full<-rswidth-rpswidth
                    rdistdiff<-rsdist-rpsdist
                    
                    
                    
                    # full range
                    
                    #ps ctrl
                    d<-dplyr::filter(df,TestPhase=="PS",Count<(-2))
                    ps1fullwidth<-abs(max(d$Ant.X)-min(d$Ant.X))
                    
                    d<-dplyr::filter(df,TestPhase=="PS",Count>=-2)
                    ps2fullwidth<-abs(max(d$Ant.X)-min(d$Ant.X))
                    
                    ctrlfullwidthdiff<-ps2fullwidth-ps1fullwidth
                    
                    
                    
                    # S change
                    dps<-dplyr::filter(df,TestPhase=="PS",Count>=-2)
                    psfullwidth<-abs(max(dps$Ant.X)-min(dps$Ant.X))
                    
                    ds<-dplyr::filter(df,TestPhase=="S")
                    sfullwidth<-abs(max(ds$Ant.X)-min(ds$Ant.X))
                    
                    
                    fullwidthdiff<-sfullwidth-psfullwidth
                    
                    
                    
                    return(data.frame("Odor"=rep(j,4),
                                      "Test"=rep(unique(df$Test),4),
                                      "TestPhase"=c(rep("PS",2),rep("S",2)),
                                      "Ant"=rep(c("L","R"),2),
                                      "width"=c(lpswidth,rpswidth,lswidth,rswidth),
                                      "width.change"=c(lctrlwidthdiff,rctrlwidthdiff,lwidthdiff.full,rwidthdiff.full),
                                      "fullwidth"=c(psfullwidth,psfullwidth,sfullwidth,sfullwidth),
                                      "fullwidth.change"=c(ctrlfullwidthdiff,ctrlfullwidthdiff,fullwidthdiff,fullwidthdiff),
                                      "meancenterdist"=c(lpsdist,rpsdist,lsdist,rsdist),
                                      "meancenterdist.change"=c(lctrldistdiff,rctrldistdiff,ldistdiff,rdistdiff)
                                      )
                           )
                    
                  })
  
  dat.distr.metrics<-dplyr::bind_rows(datlist)
  
  
  # retrieve animal ID
  animal<-data[data$Test %in% dat.distr.metrics$Test,c("Animal","Test")] #select correct animals based on test ids in dat
  animal<-animal[!duplicated(animal),] #remove duplcated rows to have 1 iteration of Test-Animal pairs
  animal<-rep(animal$Animal,each=4) #4 iterations per trial
  
  dat.distr.metrics$Animal<-animal #add animal ID to our model data
  dat.distr.metrics<-dat.distr.metrics[,c(1,11,2:10)]#rearranging
  
  
  distr.metrics.list[[i]]<-dat.distr.metrics
  names(distr.metrics.list)[i]<-j
  
}

data.centerpool<-dplyr::bind_rows(distr.metrics.list)










############################
### "Responsive" antenna ###
############################

### Select the data for the antenna with the MIN dist change:

data.min.dist<-aggregate(meancenterdist.change~Animal+Test+Odor+TestPhase,data=data.centerpool,min)
data.min.dist<-dplyr::filter(data.min.dist,TestPhase=="S") 

data.min<-data.centerpool
data.min<-split(data.min,f=data.min$Test)
data.min<-lapply(data.min,
                 FUN=function(df){
                   
                   #choose the right trial
                   d<-data.min.dist[data.min.dist$Test==unique(df$Test),] 
                   #choose the antenna with the biggest decrease in distance in S
                   ant<-df$Ant[df$meancenterdist.change %in% d$meancenterdist.change] 
                   #select the data for the correct antenna
                   d<-dplyr::filter(df,Ant==ant)
                   
                   return(d)
                   
                 })

data.resp<-dplyr::bind_rows(data.min)





### save Test+Antenna ID for later in order to create the responsive antenna density plot
responsiveAntenna.dist<-data.resp[,1:5]
responsiveAntenna.dist.col<-dplyr::filter(responsiveAntenna.dist,
                                          Odor=="Colony")
responsiveAntenna.dist.col<-aggregate(Ant ~ Test, data=responsiveAntenna.dist.col, unique)





### getting PSctrl (non-antenna-matched for dist)

data.min.ctrl.dist<-aggregate(meancenterdist.change~Animal+Test+Odor+TestPhase,data=data.centerpool,min)
data.min.ctrl.dist<-dplyr::filter(data.min.ctrl.dist,TestPhase=="PS")

# data.resp has the antenna matching control for dist and range. But we want a non-matching control for dist.
# find indices to match data.resp with data.min.ctrl.dist (non-antenna matching control)
data.resp$id<-paste(data.resp$Test,data.resp$TestPhase,data.resp$Odor,sep="_")
data.min.ctrl.dist$id<-paste(data.min.ctrl.dist$Test,data.min.ctrl.dist$TestPhase,data.min.ctrl.dist$Odor,sep="_")
ind<-match(data.min.ctrl.dist$id,data.resp$id) 
# replace into data.resp, the new, non-matching control
data.resp[ind,"meancenterdist.change"]<-data.min.ctrl.dist$meancenterdist.change

# get overall range
overallrange.dat<-data.centerpool[,c("Animal","Test","Odor","TestPhase","fullwidth.change")] #add full range
overallrange.dat<-overallrange.dat[!duplicated(overallrange.dat),]#remove duplicates
#rearrange odor order to combine
data.resp<-arrange(data.resp,Odor)
overallrange.dat<-arrange(overallrange.dat,Odor)
data.resp<-cbind(data.resp[,c("Animal","Test","Ant","Odor","TestPhase","meancenterdist.change","width.change")],
                 overallrange.dat[,"fullwidth.change"]) #combine the data
colnames(data.resp)[8]<-"fullwidth.change"



### final min data
data.model.min<-data.resp 



### save responsive antenna ID for other analyses
RespAntDist<-data.model.min[,1:3]
RespAntDist<-RespAntDist[!duplicated(RespAntDist),] #remove duplicates
# Saving the data
write.table(data,paste("../Data/RespAntDist.csv",sep=""),row.names=F,sep=";",dec=".")








################################
### "Non-Responsive" antenna ###
################################

### Select the data for the antenna with the MAX dist change (= opposite of responsive):

data.max.dist<-aggregate(meancenterdist.change~Animal+Test+Odor+TestPhase,data=data.centerpool,max)
data.max.dist<-dplyr::filter(data.max.dist,TestPhase=="S") 

data.max<-data.centerpool
data.max<-split(data.max,f=data.max$Test)
data.max<-lapply(data.max,
                 FUN=function(df){
                   
                   #choose the right trial
                   d<-data.max.dist[data.max.dist$Test==unique(df$Test),] 
                   #choose the antenna with the biggest decrease in distance in S
                   ant<-df$Ant[df$meancenterdist.change %in% d$meancenterdist.change] 
                   #select the data for the correct antenna
                   d<-dplyr::filter(df,Ant==ant)
                   
                   return(d)
                   
                 })

data.nonresp<-dplyr::bind_rows(data.max)





### getting PSctrl (non-antenna-matched for dist)

data.max.ctrl.dist<-aggregate(meancenterdist.change~Animal+Test+Odor+TestPhase,data=data.centerpool,max)
data.max.ctrl.dist<-dplyr::filter(data.max.ctrl.dist,TestPhase=="PS")

# data.nonresp has the antenna matching control for dist and range. But we want a non-matching control for dist.
# find indices to match data.nonresp with data.max.ctrl.dist (non-antenna matching control)
data.nonresp$id<-paste(data.nonresp$Test,data.nonresp$TestPhase,data.nonresp$Odor,sep="_")
data.max.ctrl.dist$id<-paste(data.max.ctrl.dist$Test,data.max.ctrl.dist$TestPhase,data.max.ctrl.dist$Odor,sep="_")
ind<-match(data.max.ctrl.dist$id,data.nonresp$id) 
# replace into data.nonresp, the new, non-matching control
data.nonresp[ind,"meancenterdist.change"]<-data.max.ctrl.dist$meancenterdist.change

# get overall range (not necessary here but avoids having to re-run the resp code to check it)
overallrange.dat<-data.centerpool[,c("Animal","Test","Odor","TestPhase","fullwidth.change")] #add full range
overallrange.dat<-overallrange.dat[!duplicated(overallrange.dat),]#remove duplicates
#rearrange odor order to combine
data.nonresp<-arrange(data.nonresp,Odor)
overallrange.dat<-arrange(overallrange.dat,Odor)
data.nonresp<-cbind(data.nonresp[,c("Animal","Test","Ant","Odor","TestPhase","meancenterdist.change","width.change")],
                 overallrange.dat[,"fullwidth.change"]) #combine the data
colnames(data.nonresp)[8]<-"fullwidth.change"



### final max data
data.model.max<-data.nonresp 











#########################
### Modelling metrics ###
#########################
# run this for each metric and for each group (resp/non-resp) by changing the parameters


# dat.model<-dplyr::filter(data.model.min,
#                        TestPhase=="S")
dat.model<-dplyr::filter(data.model.max,
                       TestPhase=="S")

#overall distribution
hist(dat.model$meancenterdist.change)
hist(dat.model$width.change)
hist(dat.model$fullwidth.change)




### Fit
set.seed(110321)
brmsmod<-brm(meancenterdist.change ~ Odor + (1|Animal),
             # width.change ~ Odor + (1|Animal),
             # fullwidth.change ~ Odor + (1|Animal),
             data=dat.model,
             family=student(),
             iter=10000, chains=4, cores=4,
             backend="cmdstanr", control=list(adapt_delta=0.99,max_treedepth=20))



### Checks

# pp check
pp<-pp_check(brmsmod,ndraws=200) + xlim(-1,1)
pp

pp_check(brmsmod,type="stat",stat="mean",ndraws=200)
pp_check(brmsmod,type="stat",stat="sd",ndraws=200)


#convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce 



#estimate vs real mean: checking for good correlation and stable variance with increasing mean
estim<-ce$Odor #estimates

predict_vs_real<-merge(estim[,c("Odor","estimate__")],dat.model)
plot(meancenterdist.change~estimate__,data=predict_vs_real)

mean_reald<-aggregate(meancenterdist.change~Odor,data=dat.model,FUN="mean") #mean per group
predict_vs_real_avg<-merge(estim[,c("Odor","estimate__")],mean_reald) #avg estimate vs real
points(meancenterdist.change~estimate__,data=predict_vs_real_avg,pch=20,col="forestgreen",cex=2)
abline(a=0,b=1)




### Probabilities

# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))

#ctrls (resp)
col.ctrl.dist<-mean(data.model.min$meancenterdist.change[data.model.min$Odor=="Colony" & data.model.min$TestPhase=="PS"])
but.ctrl.dist<-mean(data.model.min$meancenterdist.change[data.model.min$Odor=="Butanol" & data.model.min$TestPhase=="PS"])
lin.ctrl.dist<-mean(data.model.min$meancenterdist.change[data.model.min$Odor=="Linalool" & data.model.min$TestPhase=="PS"])

col.ctrl.range<-mean(data.model.min$width.change[data.model.min$Odor=="Colony" & data.model.min$TestPhase=="PS"])
but.ctrl.range<-mean(data.model.min$width.change[data.model.min$Odor=="Butanol" & data.model.min$TestPhase=="PS"])
lin.ctrl.range<-mean(data.model.min$width.change[data.model.min$Odor=="Linalool" & data.model.min$TestPhase=="PS"])

#ctrls (nonresp)
col.ctrl.dist<-mean(data.model.max$meancenterdist.change[data.model.max$Odor=="Colony" & data.model.max$TestPhase=="PS"])
but.ctrl.dist<-mean(data.model.max$meancenterdist.change[data.model.max$Odor=="Butanol" & data.model.max$TestPhase=="PS"])
lin.ctrl.dist<-mean(data.model.max$meancenterdist.change[data.model.max$Odor=="Linalool" & data.model.max$TestPhase=="PS"])

col.ctrl.range<-mean(data.model.max$width.change[data.model.max$Odor=="Colony" & data.model.max$TestPhase=="PS"])
but.ctrl.range<-mean(data.model.max$width.change[data.model.max$Odor=="Butanol" & data.model.max$TestPhase=="PS"])
lin.ctrl.range<-mean(data.model.max$width.change[data.model.max$Odor=="Linalool" & data.model.max$TestPhase=="PS"])



# Colony
ColS.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColS.quants
print(paste("P(Col<0) = ",round(mean((simval[,1]+simval[,2])<0),2),sep=""))
print(paste("P(Col<ctrl) = ",round(mean((simval[,1]+simval[,2])<col.ctrl.dist),2),sep="")) #change the control depending on metric
# Butanol
ButS.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ButS.quants
print(paste("P(But<0) = ",round(mean((simval[,1])<0),2),sep=""))
print(paste("P(But<ctrl) = ",round(mean((simval[,1])<but.ctrl.dist),2),sep="")) #change the control depending on metric
# Linalool
LinS.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
LinS.quants
print(paste("P(Lin<0) = ",round(mean((simval[,1]+simval[,3])<0),2),sep=""))
print(paste("P(Lin<ctrl) = ",round(mean((simval[,1]+simval[,3])<lin.ctrl.dist),2),sep="")) #change the control depending on metric









### Effect plot
# values have been specified by hand from the various model outputs

pdf(file=paste("./Plots/Max&Minchanges_DistRangeOverallrange_PooledS2s_EffectBarplot.pdf",sep=""),
    width=9,height=7,bg="white",useDingbats=F)


par(mfrow=c(1,3),lwd=5)


# dist to stream
m.fit<-matrix(ncol=2,nrow=3)
colnames(m.fit)<-c("max","min")
rownames(m.fit)<-c("Col","But","Lin")
m.fit[,1]<-c(0.02,0.05,0.1) #max
m.fit[,2]<-c(-0.17,-0.2,-0.17) #min

m.lwr<-matrix(ncol=2,nrow=3)
colnames(m.lwr)<-c("max","min")
rownames(m.lwr)<-c("Col","But","Lin")
m.lwr[,1]<-c(-0.02,0.02,0.06) #max
m.lwr[,2]<-c(-0.22,-0.24,-0.21) #min

m.upr<-matrix(ncol=2,nrow=3)
colnames(m.upr)<-c("max","min")
rownames(m.upr)<-c("Col","But","Lin")
m.upr[,1]<-c(0.06,0.08,0.14) #max
m.upr[,2]<-c(-0.13,-0.16,-0.14) #min

ylim=c(-0.3,0.25)

barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,las=1,
                    cex.axis=2,ylab="Delta",main="max & min change 
in distance to stream")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr)
axis(1,at=barcenters,labels=rep(c("Col","But","Lin"),2),lty=0,cex.axis=1)

points(barcenters,c(0.07,0.12,0.07,-0.12,-0.09,-0.1),pch=20,col="forestgreen",cex=0.5) #control




# antennal range
m.fit<-matrix(ncol=2,nrow=3)
colnames(m.fit)<-c("max","min")
rownames(m.fit)<-c("Col","But","Lin")
m.fit[,1]<-c(0.07,0.08,0.05) #max
m.fit[,2]<-c(0.13,0.27,0.22) #min

m.lwr<-matrix(ncol=2,nrow=3)
colnames(m.lwr)<-c("max","min")
rownames(m.lwr)<-c("Col","But","Lin")
m.lwr[,1]<-c(0,0.02,-0.02) #max
m.lwr[,2]<-c(0.04,0.21,0.15) #min

m.upr<-matrix(ncol=2,nrow=3)
colnames(m.upr)<-c("max","min")
rownames(m.upr)<-c("Col","But","Lin")
m.upr[,1]<-c(0.13,0.14,0.12) #max
m.upr[,2]<-c(0.21,0.34,0.28) #min

ylim=c(-0.15,0.4)

barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,las=1,
                    cex.axis=2,ylab="Delta",main="max & min change 
in antennal range")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr)
axis(1,at=barcenters,labels=rep(c("Col","But","Lin"),2),lty=0,cex.axis=1)

points(barcenters,c(0.01,0.02,0.04,0.03,-0.06,-0.01),pch=20,col="forestgreen",cex=0.5) #control





# overall range
m.fit<-matrix(c(-0.01,0.01,0))
rownames(m.fit)<-c("Col","But","Lin")
m.lwr<-matrix(c(-0.07,-0.03,-0.04))
rownames(m.lwr)<-c("Col","But","Lin")
m.upr<-matrix(c(0.03,0.06,0.05))
rownames(m.upr)<-c("Col","But","Lin")

ylim=c(-0.2,0.2)


par(lwd=5)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    col=rgb(0,0,0,100,maxColorValue=255),
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="Delta",main="overall range change")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,)

axis(1,at=barcenters,labels=c("Col","But","Lin"),lty=0,cex.axis=2)



dev.off()

















#--------------------------------------------------------------------------------------------------------------------------
# Side odor: Left/Right antennae metrics
#--------------------------------------------------------------------------------------------------------------------------


#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)





###############################
### Metrics for side stream ###
###############################

### data
data<-AllData

dat<-dplyr::filter(data,
                   Protocol=="ProtocolGa",
                   Count>=-4 & Count<=4)
dat<-dat[,c("Test","Protocol","Odor","Count","TestPhase",
            "Lant.X","Lant.Y","Rant.X","Rant.Y","head.X","head.Y")]

#stream coordinates
psstream<-rep(30,400)
staticstream<-rep(30,400)

#window size
winsize=4

distr.metrics.list<-list()
for(i in 1:3){
  
  j<-c("Colony","Butanol","Linalool")[i]
  
  d<-dplyr::filter(dat,Odor==j)
  
  
  # Normalization of coordinates for each trial & distribution metrics
  datlist<-split(d,f=d$Test)  #separate the df into a list by trial id
  datlist<-lapply(datlist,
                  FUN=function(df){
                    
                    #centering coordinates around the head
                    df$Lant.X<-df$Lant.X-df$head.X
                    df$Rant.X<-df$Rant.X-df$head.X
                    
                    # normalizing by maximum coordinates within trial
                    df$Lant.X<-df$Lant.X/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                    df$Rant.X<-df$Rant.X/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                    
                    #normalizing stream coordinates the same way
                    psstream<-psstream/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                    staticstream<-staticstream/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                    
                    
                    return(df)
                    
                  })
  d<-dplyr::bind_rows(datlist) #put df back together
  
  
  
  # Getting avg distances to the odor
  datlist<-split(d,f=d$Test) 
  datlist<-lapply(datlist,
                  FUN=function(df){
                    
                    
                    #Right dist & range
                    
                    #S
                    rps<-mean(abs(psstream[1:(winsize*100)]-df$Rant.X[df$Count<=0 & df$Count>(-winsize)]),na.rm=T)
                    rstatic<-mean(abs(staticstream[1:(winsize*100)]-df$Rant.X[df$Count>0 & df$Count<=winsize]),na.rm=T)
                    rstaticdiff<-rstatic-rps
                    
                    rpswidth<-max(df$Rant.X[df$Count<=0 & df$Count>(-winsize)],na.rm=T)-min(df$Rant.X[df$Count<=0 & df$Count>(-winsize)],na.rm=T)
                    rswidth<-max(df$Rant.X[df$Count>0 & df$Count<=winsize],na.rm=T)-min(df$Rant.X[df$Count>0 & df$Count<=winsize],na.rm=T)
                    rswidthdiff.full<-rswidth-rpswidth
                    
                    
                    
                    #Left dist & range
                    
                    # S
                    lps<-mean(abs(psstream[1:(winsize*100)]-df$Lant.X[df$Count<=0 & df$Count>(-winsize)]),na.rm=T)
                    lstatic<-mean(abs(staticstream[1:(winsize*100)]-df$Lant.X[df$Count>0 & df$Count<=winsize]),na.rm=T)
                    lstaticdiff<-lstatic-lps
                    
                    lpswidth<-max(df$Lant.X[df$Count<=0 & df$Count>(-winsize)],na.rm=T)-min(df$Lant.X[df$Count<=0 & df$Count>(-winsize)],na.rm=T)
                    lswidth<-max(df$Lant.X[df$Count>0 & df$Count<=winsize],na.rm=T)-min(df$Lant.X[df$Count>0 & df$Count<=winsize],na.rm=T)
                    lswidthdiff.full<-lswidth-lpswidth
                    
                    
                    
                    #overall range 
                    
                    # S
                    dps<-dplyr::filter(df,TestPhase=="PS", Count>(-winsize))
                    psfullwidth<-abs(max(c(dps$Lant.X,dps$Rant.X))-min(c(dps$Lant.X,dps$Rant.X)))
                    
                    ds<-dplyr::filter(df,TestPhase=="S", Count<=winsize)
                    sfullwidth<-max(c(ds$Lant.X,ds$Rant.X))-min(c(ds$Lant.X,ds$Rant.X))
                    
                    fullwidthdiff<-sfullwidth-psfullwidth
                    
                    
                    
                    
                    
                    dist.df<-data.frame("Odor"=rep(unique(df$Od),2),
                                        "Test"=rep(unique(df$Test),2),
                                        "Protocol"=rep(unique(df$Protocol),2),
                                        "TestPhase"=c("PS","Static"),
                                        "Rdist"=c(rps,rstatic),
                                        "Rdist.change"=c(NA,rstaticdiff),
                                        "Ldist"=c(lps,lstatic),
                                        "Ldist.change"=c(NA,lstaticdiff),
                                        "Rwidth.change"=c(NA,rswidthdiff.full),
                                        "Lwidth.change"=c(NA,lswidthdiff.full),
                                        "fullwidth.change"=c(NA,fullwidthdiff)
                    )
                    
                    return(dist.df)
                    
                  })
  dat.distr.metrics<-dplyr::bind_rows(datlist)
  
  
  # retrieve animal ID
  animal<-data[data$Test %in% dat.distr.metrics$Test,c("Animal","Test")] #select correct animals based on test ids in dat
  animal<-animal[!duplicated(animal),] #remove duplcated rows to have 1 iteration of Test-Animal pairs
  animal<-rep(animal$Animal,each=2) #2 iterations per trial
  
  dat.distr.metrics$Animal<-animal #add animal ID to our model data
  dat.distr.metrics<-dat.distr.metrics[,c(1,12,2:11)]#rearranging
  
  
  distr.metrics.list[[i]]<-dat.distr.metrics
  names(distr.metrics.list)[i]<-j
  
}

data.side<-dplyr::bind_rows(distr.metrics.list)














#########################
### Modelling metrics ###
#########################


### Distance model data
dat.dist<-data.side[,c("Animal","Odor","Test","TestPhase","Rdist.change","Ldist.change")]
# reshaping to have a response variable (delta) and an antenna variable
dat.dist<-reshape2::melt(id.vars=1:4,data=dat.dist,variable.name="Ant",value.name="delta")
dat.dist$Ant<-as.character(dat.dist$Ant)
# keep odor delta
dat.model.dist<-dplyr::filter(dat.dist,TestPhase=="Static")


### Range model data
dat.range<-data.side[,c("Animal","Odor","Test","TestPhase","Rwidth.change","Lwidth.change")]
# reshaping to have a response variable (delta) and an antenna variable
dat.range<-reshape2::melt(id.vars=1:4,data=dat.range,variable.name="Ant",value.name="delta")
dat.range$Ant<-as.character(dat.range$Ant)
# keep odor delta
dat.model.range<-dplyr::filter(dat.range,TestPhase=="Static")


### Overall range model data
dat.model.overallrange<-data.side[,c("Animal","Odor","Test","TestPhase","fullwidth.change")]
# Get odor delta
dat.model.overallrange<-dplyr::filter(dat.model.overallrange,TestPhase=="Static")
colnames(dat.model.overallrange)[5]<-"delta"







### Model

# choose data
dat.model<-dat.model.dist
# dat.model<-dat.model.range
# dat.model<-dat.model.overallrange

# look at distribution
hist(dat.model$delta)


# Fit
set.seed(110321)
brmsmod<-brm(delta ~ Odor + Ant + Odor:Ant + (1|Animal),
             # delta ~ Odor + (1|Animal), #for overall range
             data=dat.model, 
             family=student(),
             iter=10000, chains=4, cores=4,
             backend="cmdstanr", control=list(adapt_delta=0.99,max_treedepth=20))







### Checks

# pp check
pp<-pp_check(brmsmod,ndraws=200) + xlim(-1,1)
pp

pp_check(brmsmod,type="stat",stat="mean",ndraws=200)
pp_check(brmsmod,type="stat",stat="sd",ndraws=200)


#convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce 



#estimate vs real mean: checking for good correlation and stable variance with increasing mean
#change the variable to check each 
estim<-ce$Odor #estimates

predict_vs_real<-merge(estim[,c("Odor","estimate__")],dat.model)
plot(delta~estimate__,data=predict_vs_real)

mean_reald<-aggregate(delta~Odor,data=dat.model,FUN="mean") #mean per group
predict_vs_real_avg<-merge(estim[,c("Odor","estimate__")],mean_reald) #avg estimate vs real
points(delta~estimate__,data=predict_vs_real_avg,pch=20,col="forestgreen",cex=2)
abline(a=0,b=1)




### Probabilities

# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))


# Colony Left
ColSL.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColSL.quants
print(paste("P(Col L < 0) = ",round(mean((simval[,1]+simval[,2])<0),2),sep=""))
# Colony Right
ColSR.quants<-round(quantile((simval[,1]+simval[,2]+simval[,4]+simval[,5]),probs=c(0.5,0.025,0.975)),2)
ColSR.quants
print(paste("P(Col R < 0) = ",round(mean((simval[,1]+simval[,2]+simval[,4]+simval[,5])<0),2),sep=""))

# Butanol Left
ButSL.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ButSL.quants
print(paste("P(But L < 0) = ",round(mean((simval[,1])<0),2),sep=""))
# Butanol Right
ButSR.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
ButSR.quants
print(paste("P(But R < 0) = ",round(mean((simval[,1]+simval[,4])<0),2),sep=""))

# Linalool Left
LinSL.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
LinSL.quants
print(paste("P(Lin L < 0) = ",round(mean((simval[,1]+simval[,3])<0),2),sep=""))
# Linalool Right
LinSR.quants<-round(quantile((simval[,1]+simval[,3]+simval[,4]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
LinSR.quants
print(paste("P(Lin R < 0) = ",round(mean((simval[,1]+simval[,3]+simval[,4]+simval[,6])<0),2),sep=""))








### Effect plot
# values have been specified by hand from the various model outputs

#colors
lcol<-rgb(firebrick3[1],firebrick3[2],firebrick3[3],100,maxColorValue=255)
rcol<-rgb(dodgerblue3[1],dodgerblue3[2],dodgerblue3[3],100,maxColorValue=255)



pdf(file=paste("./Plots/Range&Dist&Overallrange_SideStream_EffectBarplot.pdf",sep=""),
    width=9,height=7,bg="white",useDingbats=F)


par(mfrow=c(1,3),lwd=5)


# distance %
m.fit<-matrix(ncol=2,nrow=3)
colnames(m.fit)<-c("L","R")
rownames(m.fit)<-c("Col","But","Lin")
m.fit[,1]<-c(-0.09,-0.05,-0.18) #L
m.fit[,2]<-c(-0.03,-0.04,-0.05) #R

m.lwr<-matrix(ncol=2,nrow=3)
colnames(m.lwr)<-c("L","R")
rownames(m.lwr)<-c("Col","But","Lin")
m.lwr[,1]<-c(-0.18,-0.14,-0.26) #L
m.lwr[,2]<-c(-0.1,-0.13,-0.13) #R

m.upr<-matrix(ncol=2,nrow=3)
colnames(m.upr)<-c("L","R")
rownames(m.upr)<-c("Col","But","Lin")
m.upr[,1]<-c(0,0.03,-0.09) #L
m.upr[,2]<-c(0.05,0.04,0.03) #R

ylim=c(-0.3,0.2)


barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    col=rep(c(lcol,rcol),each=3),border=rep(c("firebrick3","dodgerblue3"),each=3),
                    xaxt="n",ylim=ylim,las=1,
                    cex.axis=2,ylab="Delta",main="distance ")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col=rep(c("firebrick3","dodgerblue3"),each=3))
axis(1,at=barcenters,labels=rep(c("Col","But","Lin"),2),lty=0,cex.axis=1)




# range %
m.fit<-matrix(ncol=2,nrow=3)
colnames(m.fit)<-c("L","R")
rownames(m.fit)<-c("Col","But","Lin")
m.fit[,1]<-c(0.1,0.22,0.15) #L
m.fit[,2]<-c(-0.02,-0.03,-0.07) #R

m.lwr<-matrix(ncol=2,nrow=3)
colnames(m.lwr)<-c("L","R")
rownames(m.lwr)<-c("Col","But","Lin")
m.lwr[,1]<-c(-0.01,0.1,0.04) #L
m.lwr[,2]<-c(-0.13,-0.14,-0.17) #R

m.upr<-matrix(ncol=2,nrow=3)
colnames(m.upr)<-c("L","R")
rownames(m.upr)<-c("Col","But","Lin")
m.upr[,1]<-c(0.22,0.33,0.26) #L
m.upr[,2]<-c(0.1,0.08,0.04) #R


ylim=c(-0.25,0.4)


barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    col=rep(c(lcol,rcol),each=3),border=rep(c("firebrick3","dodgerblue3"),each=3),
                    xaxt="n",ylim=ylim,las=1,
                    cex.axis=2,ylab="Delta",main="antennal range ")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col=rep(c("firebrick3","dodgerblue3"),each=3))
axis(1,at=barcenters,labels=rep(c("Col","But","Lin"),2),lty=0,cex.axis=1)




# overall range %
m.fit<-matrix(c(0,0.02,-0.05))
rownames(m.fit)<-c("Col","But","Lin")
m.lwr<-matrix(c(-0.07,-0.04,-0.12))
rownames(m.lwr)<-c("Col","But","Lin")
m.upr<-matrix(c(0.07,0.08,0.01))
rownames(m.upr)<-c("Col","But","Lin")

ylim=c(-0.2,0.4)


par(lwd=5)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    col=rgb(0,0,0,100,maxColorValue=255),
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="Delta",main="overall range ")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,)
axis(1,at=barcenters,labels=c("Col","But","Lin"),lty=0,cex.axis=2)




dev.off()
















#--------------------------------------------------------------------------------------------------------------------------
# Antenna-specific presence density maps
#--------------------------------------------------------------------------------------------------------------------------


#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)



# Parameters
#/////////////
Od="Colony" # "Colony", "Butanol" or "Linalool"
Prot=c("ProtocolB","ProtocolGb","ProtocolE") # For centre stream, Prot = c("ProtocolB","ProtocolGb","ProtocolE"). For side stream, Prot = c("ProtocolGa")
nbins=24   #number of bins in x and y 
#/////////////

# raw data

dat<-dplyr::filter(AllData,
                   Dataset=="3D",
                   Protocol %in% Prot,
                   Odor==Od,
                   Count>=-2 & Count<=2)

dat.sub<-dat[,c("Test","Count","TestPhase","Lant.X","Lant.Y","Lant.Z",
                "Rant.X","Rant.Y","Rant.Z","head.X","head.Y","head.Z")]

#restructuring to have an antenna variable
dat<-dat.sub
colnames(dat)[4:9]<-c("Ant.X","Ant.Y","Ant.Z","Ant.X","Ant.Y","Ant.Z")
dat<-rbind(data.frame(dat[,c(1:6,10:12)],"Ant"=rep("L",nrow(dat))),
           data.frame(dat[,c(1:3,7:12)],"Ant"=rep("R",nrow(dat)))) 

# these barely move at all and heavily influences the density plots
dat<-dplyr::filter(dat,Test %!in% c("T618","T662","T786","T787","T944","T945",
                                    "T843","T946"))  


# Normalization of coordinates for each trial & distribution metrics
datlist<-split(dat,f=dat$Test)  #separate the df into a list by trial id
datlist<-lapply(datlist,
                FUN=function(df){
                  
                  #centering coordinates around the head
                  df$Ant.X<-df$Ant.X-df$head.X
                  df$Ant.Y<-df$Ant.Y-df$head.Y
                  
                  # normalizing by maximum coordinates of all trials
                  df$Ant.X<-df$Ant.X/max(c(max(dat$Ant.X-dat$head.X,na.rm=T),abs(min(dat$Ant.X-dat$head.X,na.rm=T))))
                  df$Ant.Y<-df$Ant.Y/max(c(max(dat$Ant.Y-dat$head.Y,na.rm=T),abs(min(dat$Ant.Y-dat$head.Y,na.rm=T))))
                  
                  return(df)
                  
                })
dat<-dplyr::bind_rows(datlist) #put df back together
dat<-dat[,-c(7:9)] #remove head coordinates




# LEFT & RIGHT
dat.L<-dplyr::filter(dat,Ant=="L")
dat.R<-dplyr::filter(dat,Ant=="R")





# Presence density matrices 

#LEFT
Presence.list<-list()
Hist.list<-list()
index=1
# for(i in 1:3){
for(i in 1:2){
  
  phase<-c("PS","S","PstS")[i]
  subdat<-dplyr::filter(dat.L,
                        TestPhase==phase)
  
  Hist.list[[index]]<-subdat
  names(Hist.list)[index]<-phase
  
  
  x.bin<-seq(-1,1,length=nbins)  # bin borders between min and max coordinates (normalized), equal in size
  y.bin<-seq(-1,1,length=nbins)
  
  
  presence.df<-as.data.frame(table(findInterval(subdat$Ant.X,x.bin),findInterval(subdat$Ant.Y,y.bin)))

  
  presence.df[,1]<-as.numeric(presence.df[,1])
  presence.df[,2]<-as.numeric(presence.df[,2])
  # Somehow, findInterval labels the intervals backwards? 
  # Need the biggest bin to always be max(nbins), because of how I normalise the coordinates. 
  # For the left antenna, only the Y-bins need to be fixed, because the X-bins are already all the way left
  presence.df[,2]<-presence.df[,2]+(nbins-max(presence.df[,2])) 
  
  spacemat<-matrix(data=0,nrow=nbins,ncol=nbins)
  spacemat[cbind(presence.df[,1],presence.df[,2])]<-presence.df[,3]
  
  spacemat<-spacemat/sum(spacemat)
  
  Presence.list[[index]]<-spacemat
  names(Presence.list)[index]<-phase
  index<-index+1
  
}

Presence.list.L<-Presence.list





#RIGHT
Presence.list<-list()
Hist.list<-list()
index=1

for(i in 1:2){
  
  phase<-c("PS","S")[i]
  subdat<-dplyr::filter(dat.R,
                        TestPhase==phase)
  
  Hist.list[[index]]<-subdat
  names(Hist.list)[index]<-phase
  
  
  x.bin<-seq(-1,1,length=nbins)  # bin borders between min and max coordinates (normalized), equal in size
  y.bin<-seq(-1,1,length=nbins)
  
  

  presence.df<-as.data.frame(table(findInterval(subdat$Ant.X,x.bin),findInterval(subdat$Ant.Y,y.bin)))

  
  presence.df[,1]<-as.numeric(presence.df[,1])
  presence.df[,2]<-as.numeric(presence.df[,2])
  # Somehow, findInterval labels the intervals backwards? 
  # Need the biggest bin to always be max(nbins), because of how I normalise the coordinates. 
  # For the right antenna, both axes need to be fixed
  presence.df[,1]<-presence.df[,1]+(nbins-max(presence.df[,1])) 
  presence.df[,2]<-presence.df[,2]+(nbins-max(presence.df[,2])) 
  
  spacemat<-matrix(data=0,nrow=nbins,ncol=nbins)
  spacemat[cbind(presence.df[,1],presence.df[,2])]<-presence.df[,3]
  
  spacemat<-spacemat/sum(spacemat,na.rm=T)
  
  Presence.list[[index]]<-spacemat
  names(Presence.list)[index]<-phase
  index<-index+1
  
  
}

Presence.list.R<-Presence.list





# Plotting maps 

# color scales
alpha=0.9

darkred<-col2rgb("darkred")
red3<-col2rgb("red3")
red2<-col2rgb("red2")
darkdarkred.hex<-rgb(100,darkred[2],darkred[3],alpha*255,maxColorValue=255)
red3.hex<-rgb(red3[1],red3[2],red3[3],alpha*255,maxColorValue=255)
red2light.hex<-rgb(red2[1],red2[2],red2[3],0.05*255,maxColorValue=255)

darkblue<-col2rgb("darkblue")
blue3<-col2rgb("blue3")
blue2<-col2rgb("blue2")
darkdarkblue.hex<-rgb(darkblue[1],darkblue[2],100,alpha*255,maxColorValue=255)
blue3.hex<-rgb(blue3[1],blue3[2],blue3[3],alpha*255,maxColorValue=255)
blue2light.hex<-rgb(blue2[1],blue2[2],blue2[3],0.05*255,maxColorValue=255)


cols.L<-c(red2light.hex,red3.hex,darkdarkred.hex)
colfunc.L<-colorRampPalette(cols.L,alpha=T)

cols.R<-c(blue2light.hex,blue3.hex,darkdarkblue.hex)
colfunc.R<-colorRampPalette(cols.R,alpha=T)





# pdf(file=paste("./Plots/PresenceDensity_Pooled2s_REDBLUE_topview_",Od,Prot,".pdf",sep=""),
#     width=7,height=4,bg="white",useDingbats=F)
png(filename=paste("./Plots/PresenceDensity_Pooled2s_REDBLUE_topview_",Od,Prot,".png",sep=""),
    width=7,height=4,units="in",bg="white",res=300)


par(mfrow=c(1,2))


x=x.bin
y=y.bin
ylab="y"


for(i in 1:2){
  
  
  ### LEFT
  d<-Presence.list.L[[i]]
  d[d==0]<-NA
  
  image(x,y,d,col=colfunc.L(100),
        breaks=(seq(0,max(d,na.rm=T),length.out=101)),
        xaxs="i",yaxs="i",xlab="",
        axes=F, #for a version without abox around
        sub=names(Presence.list.L)[i],cex.sub=1.5,font.sub=2)

  
  ### RIGHT
  d<-Presence.list.R[[i]]
  d[d==0]<-NA
  
  image(x,y,d,col=colfunc.R(100),
        breaks=(seq(0,max(d,na.rm=T),length.out=101)),
        add=T)

}


dev.off()







# COLOR BAR
pdf(file=paste("./Plots/PresenceDensity_Pooled2s_REDBLUE_topview_",Od,Prot,"_colorbar.pdf",sep=""),
    width=3,height=8,bg="white",useDingbats=F)


par(mfrow=c(2,2))


d<-Presence.list.L$PS
zmat<-matrix(ncol=100,nrow=1,data=seq(min(d),max(d),length.out=100))
image(x=1,y=seq(min(d),max(d),length.out=100),zmat,col=colfunc.L(100),
      ylab="",xlab="",xaxt="n",main="PS")

d<-Presence.list.R$PS
zmat<-matrix(ncol=100,nrow=1,data=seq(min(d),max(d),length.out=100))
image(x=1,y=seq(min(d),max(d),length.out=100),zmat,col=colfunc.R(100),
      ylab="",xlab="",xaxt="n",main="PS")


d<-Presence.list.L$S
zmat<-matrix(ncol=100,nrow=1,data=seq(min(d),max(d),length.out=100))
image(x=1,y=seq(min(d),max(d),length.out=100),zmat,col=colfunc.L(100),
      ylab="",xlab="",xaxt="n",main="S")

d<-Presence.list.R$S
zmat<-matrix(ncol=100,nrow=1,data=seq(min(d),max(d),length.out=100))
image(x=1,y=seq(min(d),max(d),length.out=100),zmat,col=colfunc.R(100),
      ylab="",xlab="",xaxt="n",main="S")



dev.off()







# Xcoord density curves

pdf(file=paste("./Plots/XcoordDensity_PGa_REDBLUE_",Od,Prot,".pdf",sep=""),
    width=5,height=5,bg="white",useDingbats=F)


dps.l<-dplyr::filter(dat.L,TestPhase=="PS")
dps.r<-dplyr::filter(dat.R,TestPhase=="PS")
ds.l<-dplyr::filter(dat.L,TestPhase=="S")
ds.r<-dplyr::filter(dat.R,TestPhase=="S")

dens<-density(dps.l$Ant.X,adjust=2)
plot(dens$x,dens$y,type="l",lwd=2,lty=2,
     col="firebrick3",
     ylab="Density",xlab="X",yaxs="i",xaxs="i",
     xlim=c(-1.2,1.2),
     ylim=c(0,2.5))
dens<-density(dps.r$Ant.X,adjust=2)
lines(dens$x,dens$y,lwd=2,lty=2,
      col="dodgerblue3")

dens<-density(ds.l$Ant.X,adjust=2)
lines(dens$x,dens$y,lwd=5,lty=1,
      col="firebrick3")
dens<-density(ds.r$Ant.X,adjust=2)
lines(dens$x,dens$y,lwd=5,lty=1,
      col="dodgerblue3")


dev.off()














### Density curves for resp/nonresp antennae
# Run the code in the section getting the responsive antennae IDs FIRST to have the correct data (starting line 60)
# Then run the code to get the antenna-specific densities

pdf(file=paste("./Plots/XcoordDensity_PooledS2s_Resp",Od,Prot,".pdf",sep=""),
    width=5,height=5,bg="white",useDingbats=F)


dps<-Hist.list$PS
dps$TestAnt<-paste(dps$Test,dps$Ant,sep="_")
responsiveAntenna.dist.col$TestAnt<-paste(responsiveAntenna.dist.col$Test,responsiveAntenna.dist.col$Ant,sep="_")
dps.resp<-dps[dps$TestAnt %in% responsiveAntenna.dist.col$TestAnt,] #select the antenna that is in the responsive list
dps.nonresp<-dps[dps$TestAnt %!in% responsiveAntenna.dist.col$TestAnt,] #select the antenna that is NOT in the responsive

ds<-Hist.list$S
ds$TestAnt<-paste(ds$Test,ds$Ant,sep="_")
ds.resp<-ds[ds$TestAnt %in% dps.resp$TestAnt,] #transfer the identity to S
ds.nonresp<-ds[ds$TestAnt %in% dps.nonresp$TestAnt,] #transfer the identity to S


adj=2


dps.resp$Ant.X[dps.resp$Ant=="L"]<-(-dps.resp$Ant.X[dps.resp$Ant=="L"])
dps.nonresp$Ant.X[dps.nonresp$Ant=="L"]<-(-dps.nonresp$Ant.X[dps.nonresp$Ant=="L"])
ds.resp$Ant.X[ds.resp$Ant=="L"]<-(-ds.resp$Ant.X[ds.resp$Ant=="L"])
ds.nonresp$Ant.X[ds.nonresp$Ant=="L"]<-(-ds.nonresp$Ant.X[ds.nonresp$Ant=="L"])


dens<-density(dps.resp$Ant.X,adjust=adj)
plot(dens$x,dens$y,type="l",lwd=5,
     col="dodgerblue3",
     ylab="Density",xlab="Distance to the center",yaxs="i",xaxs="i",xlim=c(-0.5,1.2),
     ylim=c(0,2))
#PS non-resp
dens<-density(dps.nonresp$Ant.X,adjust=adj)
lines(dens$x,dens$y,lwd=5,col="dodgerblue3",lty=2)
#S resp
dens<-density(ds.resp$Ant.X,adjust=adj)#
lines(dens$x,dens$y,lwd=5,col="firebrick3")
#S non-resp
dens<-density(ds.nonresp$Ant.X,adjust=adj)
lines(dens$x,dens$y,lwd=5,col="firebrick3",lty=2)



dev.off()








