###########################################################################################################################
#                                     Active smelling in the American cockroach: Figure 1                                 #
###########################################################################################################################

# Author: Antoine Hoffmann
# R version 4.2.0


#--------------------------------------------------------------------------------------------------------------------------
# Settings
#--------------------------------------------------------------------------------------------------------------------------


# Folder structure: -Analysis folder
#                       -Data (store data files)
#                       -Main (store .R script files)
#                           -Plots (save plots)






### Packages and functions ###
# run this before any analysis to make sure all necessary packages and custom functions are loaded


require("dplyr")
require("seewave")
require("grDevices")



# colors
gray40<-as.vector(col2rgb(col =c("gray40"),alpha=F ))
FränziPurple<-"#542788FF"
FränziYellow<-"#F1A340FF"
FränziPurplergb<-as.vector(col2rgb(col =c(FränziPurple),alpha=F ))
FränziYellowrgb<-as.vector(col2rgb(col =c(FränziYellow),alpha=F ))




# Function for the opposite of "%in%".
# Returns TRUE if the tested element IS NOT present in the comparison vector
"%!in%" <- function(x,y){
  !("%in%"(x,y))
}


# Speed of antennae (can be tips or angles)
ant.speed<-function(x){
  
  x<-as.numeric(x)
  speed<-abs(diff(x)/diff(1:length(x))) # absolute first derivative = mm/0.01s
  speed<-speed*100 # mm/s
  speed<-c(NA,speed) #make it the same length as the original vector
  
  return(speed)
  
}


# Function to replace NAs with last value
# Fills up NAs in the beginning with the first appearing value
repeat.before<-function(x) {   
  
  ind = which(!is.na(x))  # get positions of NON-MISSING values
  
  y<-rep(x[ind], times = diff( c(ind, length(x) + 1) ))  # repeat the values at these indices
  # diffing the indices + length yields how often they need to be repeated
  
  
  if(length(y)<length(x)){ #if the first indices of x are NAs, y will be missing these indices
    
    y<-rev(y) # reverse y to put the NAs at the end
    
    ind = which(!is.na(y)) # do the same as before on the reversed vector
    y<-rep(y[ind], times = diff( c(ind, length(x) + 1) ))
    
    y<-rev(y) #then put it back in normal order 
    
  }
  
  return(y)
}


# Function to compute a correlation/covariance matrix
# mat2 = matrix or data.frame where columns are the variales on which to compute the pairwise correlation
# if mat2 exists, compute the correlation matrix across. DIMENSIONS oF MAT1 AND MAT2 NEED TO BE EXACTLY THE SAME!
# methods are "pearson", "kendall" or "spearman". See ?cor
# returns a matrix with rows = columns = variables of the input data
cor.mat<-function(mat1,mat2=NULL,meth,covar=FALSE){
  
  cormat<-matrix(NA,ncol=ncol(mat1),nrow=ncol(mat1))
  
  if(is.null(mat2)==TRUE){ #for a single matrix
    
    if(covar==FALSE){
      
      for(i in 1:ncol(mat1)){
        for(j in 1:ncol(mat1)){
          cormat[i,j]<-cor(mat1[,i],mat1[,j],method=meth)
          cormat[j,i]<-cor(mat1[,i],mat1[,j],method=meth)
        }
      }
      
    } else {
      
      for(i in 1:ncol(mat1)){
        for(j in 1:ncol(mat1)){
          cormat[i,j]<-cov(mat1[,i],mat1[,j],method=meth)
          cormat[j,i]<-cov(mat1[,i],mat1[,j],method=meth)
        }
      }
      
    }
    
    
  } else { #to compare corresponding variables of 2 matrices
    
    if(covar==FALSE){
      
      for(i in 1:ncol(mat1)){
        for(j in 1:ncol(mat1)){
          cormat[i,j]<-cor(mat1[,i],mat2[,j],method=meth)
          cormat[j,i]<-cor(mat1[,i],mat2[,j],method=meth)
        }
      }
      
    } else {
      
      for(i in 1:ncol(mat1)){
        for(j in 1:ncol(mat1)){
          cormat[i,j]<-cov(mat1[,i],mat2[,j],method=meth)
          cormat[j,i]<-cov(mat1[,i],mat2[,j],method=meth)
        }
      }
      
    }
    
  }
  
  return(cormat)
  
}











#--------------------------------------------------------------------------------------------------------------------------
# Odour stream width (PID)
#--------------------------------------------------------------------------------------------------------------------------


# Load data
StreamWidthPID<-read.csv("../Data/StreamWidthPID.csv",sep=";",dec=".",header=T,check.names=F)


# mean and sd
StreamWidthPID.mean<-aggregate(DV~CenterDist,data=StreamWidthPID,"mean")
StreamWidthPID.sd<-aggregate(DV~CenterDist,data=StreamWidthPID,"sd")

colnames(StreamWidthPID.mean)[2]<-"mean"
StreamWidthPID.mean$sd<-StreamWidthPID.sd$DV





# Plot

# make vectors 
x=c(rev(-StreamWidthPID.mean$CenterDist),StreamWidthPID.mean$CenterDist) 
y=c(rev(StreamWidthPID.mean$mean),StreamWidthPID.mean$mean)

ysd0=c(rev(StreamWidthPID.mean$mean),StreamWidthPID.mean$mean)-c(rev(StreamWidthPID.mean$sd),StreamWidthPID.mean$sd)
ysd1=c(rev(StreamWidthPID.mean$mean),StreamWidthPID.mean$mean)+c(rev(StreamWidthPID.mean$sd),StreamWidthPID.mean$sd)


pdf(file=paste("./Plots/StreamWidthPID_points.pdf",sep=""),width=5,height=5,bg="white",useDingbats=F)

plot(x,y,
     type="l",lwd=2,
     ylim=c(0,0.35),xlim=c(-10,10),
     xaxs="i",yaxs="i",
     ylab="PID voltage",xlab="Distance from the center position [mm]",main="Odour stream width (n=5, mean +- sd)")
points(x,y,pch=20,cex=1.5)

segments(x0=x,y0=ysd0,x1=x,y1=ysd1,lwd=2)

dev.off()







#--------------------------------------------------------------------------------------------------------------------------------
# Odour pulse profile (PID)
#--------------------------------------------------------------------------------------------------------------------------------


# load raw PID trace (1kHz sample rate)
PulseShape<-read.csv("../Data/PIDpulses.csv",sep=";",dec=".",header=T,check.names=F)


# visualizing raw trace
plot(PID~Time,data=PulseShape,type="l")
lines(PulseShape$Time,(PulseShape$Trigger*0.05),col="red")


# trigger as binary
PulseShape[PulseShape$Trigger>1,"Trigger"]<-1
PulseShape[PulseShape$Trigger<1,"Trigger"]<-0

# segmenting label
#label trigger starts
PulseShape$triggerstarts<-c(diff(PulseShape$Trigger)/diff(1:length(PulseShape$Trigger)),0) #deriv
PulseShape[PulseShape$triggerstarts<=0,"triggerstarts"]<-0
PulseShape[PulseShape$triggerstarts!=0,"triggerstarts"]<-1
#lines for the starts
starts<-which(PulseShape$triggerstarts==1)
#trim df again
PulseShape<-PulseShape[starts[1]:(starts[6]-1),]
starts<-which(PulseShape$triggerstarts==1)
#label trials with the "starts" lines
PulseShape$trial<-NA
PulseShape[1:(starts[2]-1),"trial"]<-1
PulseShape[1:(starts[2]-1),"Time"]<-seq(0,5.999,0.001)
PulseShape[starts[2]:(starts[3]-1),"trial"]<-2
PulseShape[starts[2]:(starts[3]-1),"Time"]<-seq(0,5.999,0.001)
PulseShape[starts[3]:(starts[4]-1),"trial"]<-3
PulseShape[starts[3]:(starts[4]-1),"Time"]<-seq(0,5.999,0.001)
PulseShape[starts[4]:(starts[5]-1),"trial"]<-4
PulseShape[starts[4]:(starts[5]-1),"Time"]<-seq(0,5.999,0.001)
PulseShape[starts[5]:nrow(PulseShape),"trial"]<-5
PulseShape[starts[5]:nrow(PulseShape),"Time"]<-seq(0,5.999,0.001)


# mean voltage per time stamp
PulseShape.mean<-aggregate(PID~Time,data=PulseShape,"mean")
PulseShape.sd<-aggregate(PID~Time,data=PulseShape,"sd")

#Prolong the start by 1 second for the aligned baseline
m<-mean(PulseShape.mean[1:100,"PID"])
s<-sd(PulseShape.mean[1:100,"PID"])
meanappend<-rnorm(n=1000,mean=m,sd=s)
meanappend<-data.frame("Time"=seq(-1,-0.001,0.001),
                       "PID"=meanappend)
PulseShape.mean<-rbind(meanappend,PulseShape.mean)

m<-mean(PulseShape.sd[1:100,"PID"])
s<-sd(PulseShape.sd[1:100,"PID"])
sappend<-rnorm(n=1000,mean=m,sd=s)
sappend<-data.frame("Time"=seq(-1,-0.001,0.001),
                    "PID"=sappend)
PulseShape.sd<-rbind(sappend,PulseShape.sd)

PulseShape.mean$sd<-PulseShape.sd$PID

# getting the original Trigger
a<-aggregate(Trigger~Time,data=PulseShape,"mean")
PulseShape.mean$Trigger<-c(rep(0,1000),a$Trigger)



# Plot
pdf(file=paste("./Plots/AvgOdourPulse_2s.pdf",sep=""),width=6,height=5,bg="white",useDingbats=F)

plot(PID~Time,data=PulseShape.mean,
     type="l",lwd=2,xaxs="i",yaxs="i",ylim=c(0,0.35),
     xlab="Time [s]",ylab="PID Voltage [V]",main="Average 2s odour pulse (n=5, mean+-sd)")
lines(PulseShape.mean$Time,(PulseShape.mean$Trigger*10)-1,col="red")
#sd
polygon(x=c(rev(PulseShape.mean$Time),PulseShape.mean$Time),
        y=c(rev(PulseShape.mean$PID+PulseShape.mean$sd),PulseShape.mean$PID-PulseShape.mean$sd),
        col=rgb(gray40[1],gray40[2],gray40[3],alpha=70,names=NULL,maxColorValue=255),border=NA)

dev.off()








#--------------------------------------------------------------------------------------------------------------------------
# Antennal speed
#--------------------------------------------------------------------------------------------------------------------------


#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)



# settings
Prot=c("ProtocolB","ProtocolGb") #PB and PGb = center odor
anglethresh=10  #any angle below 10° will not be counted
speedthresh=1 #threshold for angular speed values in °/s

# data
data<-AllData
data<-dplyr::filter(data,
                    Protocol %in% Prot,
                    Count>=-4 & Count<=4) #Count = frame count (100fps)

data<-data[,c("Test","Protocol","Count","TestPhase","AngleLeft.V","AngleRight.V","AngleLeft.H","AngleRight.H")]

# excluding trials with no movement 
data<-dplyr::filter(data,Test %!in% c("T618","T662","T786","T787","T944","T945","T843","T946"))  

# getting angular speed per trial
datlist<-split(data,f=data$Test)
datlist<-lapply(datlist,function(df){
  
  # print(unique(df$Test))
  
  df$lspeedH<-ant.speed(df$AngleLeft.H)
  df$rspeedH<-ant.speed(df$AngleRight.H)
  df$lspeedV<-ant.speed(df$AngleLeft.V)
  df$rspeedV<-ant.speed(df$AngleRight.V)
  
  #ignoring low speed values, only interested in actual sweeps
  df$lspeedH[df$lspeedH<speedthresh]<-NA
  df$rspeedH[df$rspeedH<speedthresh]<-NA
  
  df$lspeedV[df$AngleLeft.V<anglethresh]<-NA # only interested in strokes
  df$rspeedV[df$AngleRight.V<anglethresh]<-NA
  df$lspeedV[df$lspeedV<speedthresh]<-NA
  df$rspeedV[df$rspeedV<speedthresh]<-NA
  
  return(df)
  
})

data<-dplyr::bind_rows(datlist)




# distributions of V vs H speeds

pdf(file="./Plots/HvsVspeedDensity_center_allodors.pdf",
    width=5,height=6,bg="white")

adj=2
bw=12

dens.H<-density(na.omit(c(data$lspeedH,data$rspeedH)),adjust=adj,bw=bw)
plot(dens.H$x,dens.H$y,type="l",lwd=5,
     col=FränziYellow,
     ylab="Density",xlab="Speed [°/s]",main="H vs V angular speed",
     yaxs="i",xaxs="i",xlim=c(min(dens.H$x),1000),ylim=c(0,0.007))

dens.V<-density(na.omit(c(data$lspeedV,data$rspeedV)),adjust=adj,bw=bw)
lines(dens.V$x,dens.V$y,lwd=5,col=FränziPurple)

dev.off()







#--------------------------------------------------------------------------------------------------------------------------
# Fast Fourier Transform
#--------------------------------------------------------------------------------------------------------------------------


#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)


# settings
Prot=c("ProtocolB","ProtocolGb") # center odor

# data 

data<-AllData
data<-dplyr::filter(data,
                    Protocol %in% Prot,
                    Count>=-4 & Count<=4)

data<-data[,c("Test","Protocol","Count","TestPhase","AngleLeft.V","AngleRight.V","AngleLeft.H","AngleRight.H")]

# excluding trials with no movement 
data<-dplyr::filter(data,Test %!in% c("T618","T662","T786","T787","T944","T945",
                                      "T843","T946"))  

datlist<-split(data,f=data$Test)
datlist<-lapply(datlist,function(df){
  
  # print(unique(df$Test))
  
  H<-as.data.frame(seewave::spec(df$AngleLeft.H,f=100,flim=c(0,0.005),at=NULL,norm=T,scaled=F,alab="Power",plot=F))
  V<-as.data.frame(seewave::spec(df$AngleLeft.V,f=100,flim=c(0,0.005),at=NULL,norm=T,scaled=F,alab="Power",plot=F))
  
  L<-data.frame("Test"=rep(unique(df$Test),nrow(H)),"Ant"=rep("L",nrow(H)),
                "freq"=H$x*1000,"H.power"=H$y,"V.power"=V$y)
  
  
  H<-as.data.frame(seewave::spec(df$AngleRight.H,f=100,flim=c(0,0.005),at=NULL,norm=T,scaled=F,alab="Power",plot=F))
  V<-as.data.frame(seewave::spec(df$AngleRight.V,f=100,flim=c(0,0.005),at=NULL,norm=T,scaled=F,alab="Power",plot=F))
  
  R<-data.frame("Test"=rep(unique(df$Test),nrow(H)),"Ant"=rep("R",nrow(H)),
                "freq"=H$x*1000,"H.power"=H$y,"V.power"=V$y)
  
  spectrograms<-rbind(L,R)
  
  return(spectrograms)
  
})

data<-dplyr::bind_rows(datlist)





# Plot

n<-length(unique(data$Test))

# H mean and SEM vectors
H<-aggregate(H.power~freq,data=data,mean)
Hsd<-aggregate(H.power~freq,data=data,sd)
Hsem<-(Hsd$H.power/sqrt(n))
H.polygon.x<-c(H$freq,rev(H$freq))
H.polygon.y1<-c(H$H.power+Hsem,rev(H$H.power))
H.polygon.y2<-c(H$H.power-Hsem,rev(H$H.power))
H.polygon.y2[H.polygon.y2<=0]<-1e-5

# V mean and SEM vectors
V<-aggregate(V.power~freq,data=data,mean)
Vsd<-aggregate(V.power~freq,data=data,sd)
Vsem<-(Vsd$V.power/sqrt(n))
V.polygon.x<-c(V$freq,rev(V$freq))
V.polygon.y1<-c(V$V.power+Vsem,rev(V$V.power))
V.polygon.y2<-c(V$V.power-Vsem,rev(V$V.power))
V.polygon.y2[V.polygon.y2<=0]<-1e-5




pdf(file=paste("./Plots/avgFFT_HvsV_PBPGb_allodors_SEM_loglog.pdf",sep=""),
    width=5,height=5,bg="white",useDingbats=F)



plot(H.power~freq,data=H,type="l",col=FränziYellow,
     xaxs="i",
     yaxs="i",
     ylim=c(1e-5,1),
     log="xy",
     ylab="Norm. power")

par(lwd=2)

polygon(x=H.polygon.x,
        y=H.polygon.y1,
        col=rgb(FränziYellowrgb[1],FränziYellowrgb[2],FränziYellowrgb[3],alpha=100,names=NULL,maxColorValue=255),
        border=NA)
polygon(x=H.polygon.x,
        y=H.polygon.y2,
        col=rgb(FränziYellowrgb[1],FränziYellowrgb[2],FränziYellowrgb[3],alpha=100,names=NULL,maxColorValue=255),
        border=NA)
lines(H.power~freq,data=H,type="l",
      col=FränziYellow)


polygon(x=V.polygon.x,
        y=V.polygon.y1,
        col=rgb(FränziPurplergb[1],FränziPurplergb[2],FränziPurplergb[3],alpha=100,names=NULL,maxColorValue=255),
        border=NA)
polygon(x=V.polygon.x,
        y=V.polygon.y2,
        col=rgb(FränziPurplergb[1],FränziPurplergb[2],FränziPurplergb[3],alpha=100,names=NULL,maxColorValue=255),
        border=NA)
lines(V.power~freq,data=V,type="l",
      col=FränziPurple)




dev.off()












#--------------------------------------------------------------------------------------------------------------------------
# Pairwise leg correlation (tripod gait)
#--------------------------------------------------------------------------------------------------------------------------


data<-AllData  
data<-data[,c("Lfrontleg.Y","Lmidleg.Y","Lhindleg.Y","Rfrontleg.Y","Rmidleg.Y","Rhindleg.Y",
              "Lfrontleg.X","Lmidleg.X","Lhindleg.X","Rfrontleg.X","Rmidleg.X","Rhindleg.X")]


#fill up NAs
data<-apply(data,MARGIN=2,FUN=repeat.before)
#1D speed
data<-apply(data,MARGIN=2,FUN=function(x){ diff(x)/diff(1:length(x))})


#find rows with no 0
r<-apply(data,1,function(row){ all(row!=0)}) 
data<-data[r,] #keep only those


#only Y speed (forward motion)
data<-data[,1:6] 
cmat<-cor.mat(mat1=data,meth="pearson") #correlation matrix
xl<-"correlation matrix legs Y speed"


# color scale
colfunc<-colorRampPalette(c("darkblue","darkblue","white","firebrick4","firebrick4"))

# plot
image(0:ncol(cmat),0:ncol(cmat),cmat,xaxs="i",yaxs="i",
      col=colfunc(50),breaks=seq(-1,1,length.out=51),
      main="",axes=F,ylab="",xlab=xl,cex.lab=1.5)
axis(side=3,labels=colnames(data),at=seq(0.5,5.5,1))
axis(side=2,labels=colnames(data),at=seq(0.5,5.5,1))
box(which="plot",lty="solid")


# color scale
zmat<-matrix(ncol=100,nrow=1,data=seq(-1,1,length.out=100))
image(x=1,y=seq(-1,1,length.out=100),zmat,col=colfunc(100),
      ylab="Pearson Coefficient",xlab="",xaxt="n")








