###########################################################################################################################
#                               Active smelling in the American cockroach: Fig. 5 & Fig. S3                               #
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

require("grDevices")
require("brms")
require("cmdstanr")
require("bayesplot")
require("ggpubr")
require("dplyr")
require("reshape2")
require("ghibli")
require("WaveletComp")
require("abind")


# colors
color_scheme_set("purple")











#--------------------------------------------------------------------------------------------------------------------------
# Generate wavelet data
#--------------------------------------------------------------------------------------------------------------------------
# Run only once to generate the wavelet data, then use the created table for subsequent analyses.
# Or just import the wavelet data to skip this, it's long.


#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)




####################
# Wavelet analysis #
####################

# Note: 
# The power calculation is on an absolute scale, so they are comparable across time, frequency and trials.
# HOWEVER, for the plots, the colors are set by the quantiles! So red in one trial is NOT equal to red in another.
#//////////
data=dplyr::filter(AllData, Protocol=="ProtocolB")
keep=TRUE  # if FALSE, doesn't save the spectrograms (in case you just want to plot)
Plot=FALSE  # Long! Creates the raw spectrograms for each trial
#//////////

# raw data and plots for each trial ----
Powerdatlist<-split(data,f=data$Test)
Powerdatlist<-lapply(Powerdatlist,
                     FUN=function(dat){
                       
                       # Note to analyze.wavelet():
                       # dt= number of seconds per observation (here: 100 fr/sec, so 0.01)
                       # dj= period axis resoltution
                       
                       #left horizontal
                       dat2<-dplyr::select(dat,"AngleLeft.H")
                       wavelet<-analyze.wavelet(dat2,"AngleLeft.H",loess.span=0,dt=1/100,dj=1/100,lowerPeriod=0.06,upperPeriod=64,make.pva=F,verbose=F)
                       Powerdat<-as.data.frame(wavelet$Power) # columns=frames, lines=periods
                       Powerdat<-cbind("Freq"=(1/wavelet$Period),Powerdat) #adding frequency values
                       PowerdatLH<-Powerdat[order(nrow(Powerdat):1),] #invert row order
                       #right horizontal
                       dat2<-dplyr::select(dat,"AngleRight.H")
                       wavelet<-analyze.wavelet(dat2,"AngleRight.H",loess.span=0,dt=1/100,dj=1/100,lowerPeriod=0.06,upperPeriod=64,make.pva=F,verbose=F)
                       Powerdat<-as.data.frame(wavelet$Power) 
                       Powerdat<-cbind("Freq"=(1/wavelet$Period),Powerdat) 
                       PowerdatRH<-Powerdat[order(nrow(Powerdat):1),] 
                       
                       
                       # This dataset contains early trials where 3D tracking wasn't available yet, so no vertical angles
                       
                       if(isTRUE(unique(dat$Dataset)=="3D")){ # If 3D data is available
                         
                         #left vertical
                         dat2<-dplyr::select(dat,"AngleLeft.V")
                         wavelet<-analyze.wavelet(dat2,"AngleLeft.V",loess.span=0,dt=1/100,dj=1/100,lowerPeriod=0.06,upperPeriod=64,make.pva=F,verbose=F)
                         Powerdat<-as.data.frame(wavelet$Power) 
                         Powerdat<-cbind("Freq"=(1/wavelet$Period),Powerdat) 
                         PowerdatLV<-Powerdat[order(nrow(Powerdat):1),] 
                         #right vertical
                         dat2<-dplyr::select(dat,"AngleRight.V")
                         wavelet<-analyze.wavelet(dat2,"AngleRight.V",loess.span=0,dt=1/100,dj=1/100,lowerPeriod=0.06,upperPeriod=64,make.pva=F,verbose=F)
                         Powerdat<-as.data.frame(wavelet$Power) 
                         Powerdat<-cbind("Freq"=(1/wavelet$Period),Powerdat) 
                         PowerdatRV<-Powerdat[order(nrow(Powerdat):1),] 
                         
                         #df with H & V
                         Antenna<-rep(c(rep("L",nrow(Powerdat)),rep("R",nrow(Powerdat))),2)
                         Animal<-rep(unique(dat$Animal),4*nrow(Powerdat))
                         Test<-rep(unique(dat$Test),4*nrow(Powerdat))
                         Odor<-rep(unique(dat$Odor),4*nrow(Powerdat))
                         Protocol<-rep(unique(dat$Protocol),4*nrow(Powerdat))
                         Plane<-rep(c(rep("H",nrow(Powerdat)*2),rep("V",nrow(Powerdat)*2)))
                         Set<-rep("3D",4*nrow(Powerdat))

                         
                         Pdat<-data.frame(Set,
                                          Animal,
                                          Test,
                                          Odor,
                                          Protocol,
                                          Antenna,
                                          Plane,
                                          rbind(PowerdatLH,PowerdatRH,PowerdatLV,PowerdatRV)) #columns=frames
                         # THE ORDER OF RBINDING MATTERS AND HAS TO MATCH THE ANTENNA AND PLANE VECTORS
                         
                         
                       } else { # if only horizontal data is available
                         
                         #df with H only
                         Antenna<-c(rep("L",nrow(Powerdat)),rep("R",nrow(Powerdat)))
                         Animal<-rep(unique(dat$Animal),2*nrow(Powerdat))
                         Test<-rep(unique(dat$Test),2*nrow(Powerdat))
                         Odor<-rep(unique(dat$Odor),2*nrow(Powerdat))
                         Protocol<-rep(unique(dat$Protocol),2*nrow(Powerdat))
                         Plane<-rep("H",nrow(Powerdat)*2)
                         Set<-rep("2D",2*nrow(Powerdat))
                         
                         
                         Pdat<-data.frame(Set,
                                          Animal,
                                          Test,
                                          Odor,
                                          Protocol,
                                          Antenna,
                                          Plane,
                                          rbind(PowerdatLH,PowerdatRH)) #columns=frames
                         # THE ORDER OF RBINDING MATTERS AND HAS TO MATCH THE ANTENNA AND PLANE VECTORS
                         
                       }
                       
                       

                       # Plots 
                       
                       if(Plot==T){
                         
                         print("plotting...")
                         
                         timescale<-c(-7,11,18)
                         stimoff<-6
                         test<-unique(Pdat$Test)
                         
                         png(filename = paste("./Plots/RawWavelets_",test,"_Colony.png",sep=""),width=20,height=20,units="in",bg="white",res=300) 
                         
                         
                         par(mfrow=c(2,2),mar=c(4,4,2,2))
                         
                         #horizontal
                         image(x=dat$Count,y=c(PowerdatLH[,1]),z=t(PowerdatLH[,-1]),col=rev(rainbow(100,start=0,end=0.7)),
                               breaks=quantile(as.matrix(PowerdatLH[,-1]),probs=seq(0,1,length.out=100+1)),
                               log="y",xaxp=timescale,ylim=c(0.05,15),xlab="Times (s)",ylab="Frequency (Hz)",
                               main=paste(test,unique(Odor),unique(Protocol),"- Power spectrum - Left horizontal"))
                         abline(h=c(0.25,0.5,1,2,3,4,5,10),col="white")
                         abline(v=0,lty=2,lwd=2,col="black")
                         abline(v=stimoff,lty=2,lwd=2,col="black")
                         
                         image(x=dat$Count,y=c(PowerdatRH[,1]),z=t(PowerdatRH[,-1]),col=rev(rainbow(100,start=0,end=0.7)),
                               breaks=quantile(as.matrix(PowerdatRH[,-1]),probs=seq(0,1,length.out=100+1)),
                               log="y",xaxp=timescale,ylim=c(0.05,15),xlab="Times (s)",ylab="Frequency (Hz)",
                               main="Right horizontal")
                         abline(h=c(0.25,0.5,1,2,3,4,5,10),col="white")
                         abline(v=0,lty=2,lwd=2,col="black")
                         abline(v=stimoff,lty=2,lwd=2,col="black")
                         
                         
                         if(isTRUE(unique(dat$Dataset)=="3D")){
                           
                           #vertical
                           image(x=dat$Count,y=c(PowerdatLV[,1]),z=t(PowerdatLV[,-1]),col=rev(rainbow(100,start=0,end=0.7)),
                                 breaks=quantile(as.matrix(PowerdatLV[,-1]),probs=seq(0,1,length.out=100+1)),
                                 log="y",xaxp=timescale,ylim=c(0.05,15),xlab="Time (s)",ylab="Frequency (Hz)",
                                 main="vertical")
                           abline(h=c(0.25,0.5,1,2,3,4,5,10),col="white")
                           abline(v=0,lty=2,lwd=2,col="black")
                           abline(v=stimoff,lty=2,lwd=2,col="black")
                           
                           image(x=dat$Count,y=c(PowerdatRV[,1]),z=t(PowerdatRV[,-1]),col=rev(rainbow(100,start=0,end=0.7)),
                                 breaks=quantile(as.matrix(PowerdatRV[,-1]),probs=seq(0,1,length.out=100+1)),
                                 log="y",xaxp=timescale,ylim=c(0.05,15),xlab="Time (s)",ylab="Frequency (Hz)",
                                 main="vertical")
                           abline(h=c(0.25,0.5,1,2,3,4,5,10),col="white")
                           abline(v=0,lty=2,lwd=2,col="black")
                           abline(v=stimoff,lty=2,lwd=2,col="black")
                           
                         }
                         
                         
                         dev.off()
                         
                       }

                       
                       
                       
                       
                       # Binned wavelets (time & frequency bins) & removing edges
                       meta=8 #number of metadata columns
                       
                       for(p in unique(Pdat$Plane)){ #for each available plane
                         
                         # removing edge-effect of the wavelets
                         tstart<-100+meta+1 #removing the first second (keeping metadata)
                         tend<-ncol(Pdat)-100 #removing the last second
                         
                         Pdat_avg<-Pdat[Pdat$Plane==p,c(1:meta,tstart:tend)] #removing edges but keeping metadata
                         Pdat_avg<-Pdat_avg[Pdat_avg$Freq>0.15,] # removing frequencies <= 0.15
                         Pdat_avg<-Pdat_avg[Pdat_avg$Freq<10,] # removing frequencies >= 10
                         
                         #timebins
                         binnedPdat_avg<-Pdat_avg[,1:meta] #initializing the new df with the metadata to then cbind the means in bins 
                         midframe<-seq(-6.5,((ncol(Pdat_avg)-meta)/100)-7.5,1) #middle of the bins always go from -6.5 in increments of 1, with 0= stim onset 
                         k=meta+1
                         j=1
                         while((k<=ncol(Pdat_avg)) & j<=length(midframe)){ #number of seconds of the recording without borders
                           
                           bin<-Pdat_avg[,k:(k+100-1)] #bins of 1s
                           bin<-apply(bin,1,mean,na.rm=T) #mean
                           binnedPdat_avg<-cbind(binnedPdat_avg,bin)
                           
                           k=k+100 #increments of 1s
                           j=j+1
                         }
                         colnames(binnedPdat_avg)[-(1:meta)]<-midframe
                         Pdat_avg<-binnedPdat_avg #replacing the old data
                         
                         #frequency bins
                         bins<-c(0.15,0.5,seq(1,9,1),10) # intervals for frequency bins: 0.15 0.50  1.00  2.00  3.00  4.00  5.00  6.00  7.00  8.00  9.00  10.00
                         Pdat_avg$freqband<-cut(Pdat_avg$Freq,breaks=bins) #label rows with defined freqbands (extreme breaks need to be outside of the min and max of the data)
                         Pdat_avg<-Pdat_avg[,c(1:meta,ncol(Pdat_avg),(meta+1):(ncol(Pdat_avg)-1))] #rearranging columns
                         frames<-colnames(Pdat_avg)[-(1:(meta+1))]
                         
                         # Averaging in freqbands 
                         datlist<-list()
                         for(k in (meta+2):ncol(Pdat_avg)){ # start after the metadata 
                           df<-aggregate(Pdat_avg[,k]~Set+Animal+Test+Antenna+Plane+Protocol+Odor+freqband,data=Pdat_avg,FUN=mean,na.rm=T) 
                           #average in freqbands for every frame and keeping the "metadata"
                           colnames(df)[ncol(df)]<-colnames(Pdat_avg)[k]
                           
                           if(k==(meta+2)){
                             datlist[[k-meta]]<-df #for the first frame, keep the metadata columns
                           } else {
                             datlist[[k-meta]]<-as.data.frame(df[,ncol(df)]) #for the next, no need for the metadata
                           }
                         }
                         Pdat_avg<-dplyr::bind_cols(datlist) #cbind the list
                         colnames(Pdat_avg)[(meta+1):ncol(Pdat_avg)]<-frames
                         
                         #saving 
                         if(p=="H"){
                           Pdat_H<-Pdat_avg
                         } else if(p=="V"){
                           Pdat_V<-Pdat_avg
                         }
                         
                       }
                       
                       
                       #combine H and V back together and replace the old, unbinned data
                       if(isTRUE(unique(dat$Dataset)=="3D")){
                         Pdat<-rbind(Pdat_H,Pdat_V) 
                       } else{
                         Pdat<-rbind(Pdat_H) 
                       }
                       
                       
                       #saving in the list
                       if(keep==T){
                         return(Pdat) 
                       }
                       
                     })


# Creating data.frame of wavelet power
Powerdat_all_bins<-dplyr::bind_rows(Powerdatlist)




#saving table for later
write.table(Powerdat_all_bins,"../Data/Powerdat_all_bins_activesmelling.csv",row.names=F,sep=",",dec=".") #with a comma here









#--------------------------------------------------------------------------------------------------------------------------
# Average Changes in Power
#--------------------------------------------------------------------------------------------------------------------------


#load data set
Powerdat_all_bins<-read.csv("../Data/Powerdat_all_bins_activesmelling.csv",sep=",",dec=".",header=T,check.names=F) #with a comma here






#####################
# Average Power map #
#####################
# scaled change in power in frequency bands

#//////////
data=Powerdat_all_bins
# prebinned=TRUE # is the data already binned in time and frequency?
meta=8 #number of meta data columns
#//////////

# Data
Powerdat.norm<-data
# Scaling changes in power
# Compare power change in each frequency band
tempdf<-Powerdat.norm[,(meta+1):ncol(Powerdat.norm)] #using only the power data
tempdf<-tempdf[,-c(1,2)] #removing first 2s of PS

zerocol=5 # last column before stimulus
#applying the scaling to each row
tempdf<-as.data.frame(t(apply(tempdf,MARGIN=1, #each row = each unique combination of Animal+Test+Antenna+Plane+Protocol+Odor+freqband
                              FUN=function(x){
                                x<-(x-mean(x[1:zerocol],na.rm=T))/mean(x[1:zerocol],na.rm=T) #normalize by mean(PS) = similar to Delta F/F
                              })))
#reassembling with metadata
Powerdat.norm<-cbind(Powerdat.norm[,1:meta],tempdf)



### averaging power within conditions ###

#creating unique ids for the combinations of plane+odor
Powerdat.norm$id<-paste(Powerdat.norm$Plane,Powerdat.norm$Odor,sep=".") 
Powerdat.norm<-Powerdat.norm[,c(ncol(Powerdat.norm),1:(ncol(Powerdat.norm)-1))] #rearranging

# REMOVING FREQ<0.5
Powerdat.norm<-dplyr::filter(Powerdat.norm,freqband!="(0.15,0.5]")

#getting frequency band ids
freqbands<-unique(Powerdat.norm$freqband)

index=1 #initializing the index to fill the final list
AvgWavelet.list<-list()
ids<-unique(Powerdat.norm$id)

for(i in ids){
  
  print(i)
  dat<-subset(Powerdat.norm,id==i)
  
  dat$ATid<-paste(dat$Test,dat$Antenna,sep="") #unique test+antenna ID
  ATids<-unique(dat$ATid)
  dat<-dat[,-(1:9)] #only ATid and power data
  dat<-dat[,c(ncol(dat),1:(ncol(dat)-1))] #rearrange 
  
  #creating a list with individual traces 
  datlist<-list()
  for(j in 1:length(ATids)){ #looping through unique ids of test+antenna
    dat2<-subset(dat,ATid==ATids[j],select=-ATid) #only keeping the power values
    datlist[[j]]<-dat2
    names(datlist)[j]<-ATids[j] # keep Ant+Test ids
  }
  
  
  # 3D array with all the traces to then average each cell over the third dimension
  array<-abind(datlist,along=3)
  
  #average regardless of antenna identity
  dat3<-as.data.frame(apply(array,1:2,mean,na.rm=T)) #averaging the cells in the third dim (=trials)
  dat3$freqband<-freqbands
  dat3<-dat3[,c(ncol(dat3),1:(ncol(dat3)-1))] #rearranging
  
  #filling the final list
  AvgWavelet.list[[index]]<-dat3 
  names(AvgWavelet.list)[index]<-i
  index=index+1 #incrementing the index for the list
  
}









### Scaled power maps ###

for(i in 1:length(AvgWavelet.list)){

  
png(filename=paste("./Plots/AvgWavelet_Scaled_timebins_",names(AvgWavelet.list)[i],".png",sep=""),width=10,height=5,units="in",bg="white",res=300)


  layout(matrix(c(1,2,3),1,3,byrow=T),widths=c(10,1,0.2),heights=c(1,1),respect=F)


  dat<-AvgWavelet.list[[i]]

  dat1<-as.matrix(dat[,-1])
  dat1<-log(dat1+1)   # Transforming the same as for the modelling later on 
  dat1[dat1>=log(1000+1)]<-log(1000+1) #clip everything above 1000 for better color scaling
  top<-max(dat1)
  bottom<-min(dat1)

  dat1<-as.data.frame(dat1)
  dat<-cbind(dat[,1],dat1)
  dat[,1]<-as.character(dat[,1])
  colnames(dat)[1]<-"freqband"

  timescale<-seq(-4.5,9.5,1)
  ticks<-c(-5,10,length(timescale))
  stimoff<-6


  cols<-hcl.colors(100, "YlOrRd", rev = TRUE)
  cols[1]<-"aliceblue"


  par(mai=c(1.02,0.92,0.82,0))
  image(x=timescale,y=0:10,z=t(dat[,-1]),useRaster=T,
        col=cols,
        breaks=c(bottom,seq(log(1),top,length.out=100)),
        yaxt="n",xaxp=ticks,xlab="Time [s]",ylab="Frequency bands [Hz]",
        main=paste(names(AvgWavelet.list)[i],"- mean change in wavelet power"),
        cex.main=2,cex.lab=2,cex.axis=2)

  axis(side=2,at=(0:9)+0.5,labels=dat$freqband,las=2,cex.lab=3)
  box()
  abline(v=0,lty=2,lwd=3,col="black")
  abline(v=stimoff,lty=2,lwd=3,col="black")


  #legend
  par(mai=c(1.02,0.2,0.82,0.3))
  image(x=1,
        y=c(bottom,seq(log(1),top,length.out=100)),z=t(c(bottom,seq(log(1),top,length.out=100))),
        col=cols,
        breaks=c(bottom,seq(log(1),top,length.out=100)),
        axes=F,xlab="",ylab="")
  axis(side=4,
       at=log(c(0,1,10,100,1000)+1),labels=c(0,1,10,100,1000),
       tick=T,las=2,cex.lab=3,cex.axis=1,tcl=-0.5)



  
  dev.off()
  
  
}












#################
### Modelling ###
#################

#//////////
Od="Colony" #choose odor
plane="H" #choose plane
#/////////


# First: restructuring the data created above (Powerdat.norm)
Powerdat.model<-dplyr::filter(Powerdat.norm,
                              Odor==Od,
                              Plane==plane)

# Time as a variable
Powerdat.model<-reshape2::melt(id.vars=colnames(Powerdat.model)[1:9],data=Powerdat.model,variable.name="timepoint",value.name="pow") #meta=5 cols
Powerdat.model$timepoint<-as.numeric(as.character(Powerdat.model$timepoint))
Powerdat.model$TestPhase[Powerdat.model$timepoint<0]<-"PS"
Powerdat.model$TestPhase[Powerdat.model$timepoint>0 & Powerdat.model$timepoint<6]<-"S" 
Powerdat.model$TestPhase[Powerdat.model$timepoint>6]<-"PstS"


# Averaging to have one value per testphase+freqband+test+antenna
#Bigger frequency bands
Powerdat.model$Newfreqband<-as.character(Powerdat.model$freqband)
Powerdat.model<-dplyr::filter(Powerdat.model,freqband!="(0.15,0.5]") #removing freqs<0.5
Powerdat.model$Newfreqband[Powerdat.model$Newfreqband=="(1,2]" |
                             Powerdat.model$Newfreqband=="(2,3]"] <-"(1,3]"
Powerdat.model$Newfreqband[Powerdat.model$Newfreqband=="(3,4]" |
                                  Powerdat.model$Newfreqband=="(4,5]"] <-"(3,5]"
Powerdat.model$Newfreqband[Powerdat.model$Newfreqband=="(5,6]" |
                                  Powerdat.model$Newfreqband=="(6,7]"|
                                  Powerdat.model$Newfreqband=="(7,8]"|
                                  Powerdat.model$Newfreqband=="(8,9]"|
                                  Powerdat.model$Newfreqband=="(9,10]"] <-"(5,10]"
# average per new freqband
Powerdat.model<-aggregate(pow~Test+Antenna+TestPhase+Newfreqband,FUN=mean,data=Powerdat.model)

#removing PS, because it's 0 anyway
Powerdat.model<-Powerdat.model[Powerdat.model$TestPhase!="PS",]
# making groups for the model
Powerdat.model$group<-paste(Powerdat.model$TestPhase,Powerdat.model$Newfreqband)





# data transformation
trans<-function(x){ log(x+1) }
backtrans<-function(x){ exp(x)-1 }


hist(trans(Powerdat.model$pow))



# Fit
set.seed(210621)
mod<-brm(trans(pow) ~ group + (1|Test/Antenna),
         data=Powerdat.model,
         family=student(),
         iter=10000, chains=4, cores=4, backend="cmdstanr", 
         control=list(adapt_delta=0.99,max_treedepth=20))




# pp check
pp<-pp_check(mod,ndraws=200) + xlim(-20,20)
pp

pp_check(mod,type="stat",stat="mean",ndraws=200)
pp_check(mod,type="stat",stat="min",ndraws=200)
pp_check(mod,type="stat",stat="max",ndraws=200)


# convergence check
plot(mod)



#quick effect plot
ce<-conditional_effects(mod)
ce


#avg estimate vs real means 
estim<-ce$group #estimates
predict_vs_real<-merge(estim[,c("group","estimate__")],Powerdat.model) #full data
plot(trans(pow)~estimate__,data=predict_vs_real)

mean_reald<-aggregate(trans(pow)~group,data=Powerdat.model,FUN="mean") #mean per group
predict_vs_real_avg<-merge(estim[,c("group","estimate__")],mean_reald)
points(`trans(pow)`~estimate__,data=predict_vs_real_avg,pch=20,col="forestgreen",cex=2)
abline(a=0,b=1)



#quick back-transformed effect plot
ce<-conditional_effects(mod)
ce$group$estimate__<-backtrans(ce$group$estimate__)
ce$group$lower__<-backtrans(ce$group$lower__)
ce$group$upper__<-backtrans(ce$group$upper__)
ce






# MCMC samples
simval<-as.data.frame(as_draws_df(mod))


# Estimates, 95% CrI, and probabilities
quants.psts051<-round(quantile(backtrans(simval[,1]),probs=c(0.5,0.025,0.975)),2)
quants.psts051
print(paste("P(PstS0.51 > 0) = ",round(mean(backtrans(simval[,1])>0),2),sep=""))

quants.psts13<-round(quantile(backtrans(simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
quants.psts13
print(paste("P(PstS13 > 0) = ",round(mean(backtrans(simval[,1]+simval[,2])>0),2),sep=""))

quants.psts35<-round(quantile(backtrans(simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2) 
quants.psts35
print(paste("P(PstS35 > 0) = ",round(mean(backtrans(simval[,1]+simval[,3])>0),2),sep=""))

quants.psts510<-round(quantile(backtrans(simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
quants.psts510
print(paste("P(PstS510 > 0) = ",round(mean(backtrans(simval[,1]+simval[,4])>0),2),sep=""))


quants.s051<-round(quantile(backtrans(simval[,1]+simval[,5]),probs=c(0.5,0.025,0.975)),2) 
quants.s051
print(paste("P(S0.51 > 0) = ",round(mean(backtrans(simval[,1]+simval[,5])>0),2),sep=""))

quants.s13<-round(quantile(backtrans(simval[,1]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
quants.s13
print(paste("P(S13 > 0) = ",round(mean(backtrans(simval[,1]+simval[,6])>0),2),sep=""))

quants.s35<-round(quantile(backtrans(simval[,1]+simval[,7]),probs=c(0.5,0.025,0.975)),2) 
quants.s35
print(paste("P(S35 > 0) = ",round(mean(backtrans(simval[,1]+simval[,7])>0),2),sep=""))

quants.s510<-round(quantile(backtrans(simval[,1]+simval[,8]),probs=c(0.5,0.025,0.975)),2) 
quants.s510
print(paste("P(S510 > 0) = ",round(mean(backtrans(simval[,1]+simval[,8])>0),2),sep=""))




# Effect plots

#specifying colors with transparency
cols<-col2rgb(ghibli_palettes$LaputaMedium[2:5],alpha=T)
cols[4,]<-rep(100,4) #alpha
cols<-rgb(cols[1,],cols[2,],cols[3,],cols[4,],maxColorValue=255)



pdf(file=paste("./Plots/WaveletPower_Effectplot_",Od,plane,".pdf",sep=""),
  width=7,height=7,bg="white",useDingbats=F)


par(lwd=3)


m.fit<-matrix(ncol=2,nrow=4)
colnames(m.fit)<-c("S","PstS")
rownames(m.fit)<-c("0.5-1","1-3","3-5","5-10")
m.fit[,1]<-c(quants.s051[1],quants.s13[1],quants.s35[1],quants.s510[1]) 
m.fit[,2]<-c(quants.psts051[1],quants.psts13[1],quants.psts35[1],quants.psts510[1]) 

m.lwr<-matrix(ncol=2,nrow=4)
colnames(m.lwr)<-c("S","PstS")
rownames(m.lwr)<-c("0.5-1","1-3","3-5","5-10")
m.lwr[,1]<-c(quants.s051[2],quants.s13[2],quants.s35[2],quants.s510[2]) 
m.lwr[,2]<-c(quants.psts051[2],quants.psts13[2],quants.psts35[2],quants.psts510[2]) 

m.upr<-matrix(ncol=2,nrow=4)
colnames(m.upr)<-c("S","PstS")
rownames(m.upr)<-c("0.5-1","1-3","3-5","5-10")
m.upr[,1]<-c(quants.s051[3],quants.s13[3],quants.s35[3],quants.s510[3]) 
m.upr[,2]<-c(quants.psts051[3],quants.psts13[3],quants.psts35[3],quants.psts510[3]) 


ylim=c(min(m.lwr),max(m.upr))
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,las=1,
                    col=cols, border=ghibli_palettes$LaputaMedium[2:5],
                    cex.axis=1,ylab="Delta P / P",
                    main=paste(Od,plane))
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col=ghibli_palettes$LaputaMedium[2:5])
axis(1,at=barcenters,labels=rep(c("0.5-1","1-3","3-5","5-10"),2),lty=0,cex.axis=1)




dev.off()














