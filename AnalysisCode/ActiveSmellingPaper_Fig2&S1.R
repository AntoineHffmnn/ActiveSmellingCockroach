###########################################################################################################################
#                                Active smelling in the American cockroach: Fig. 2 & Fig. S1                              #
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
require("reshape2")
require("brms")
require("cmdstanr")
require("bayesplot")
require("ggpubr")
require("grDevices")
require("viridis")


color_scheme_set("purple")



FränziOrange<-"#E66101FF"
FränziPurple<-"#542788FF"
FränziYellow<-"#F1A340FF"
dodgerblue1<-as.vector(col2rgb(col=c("dodgerblue1"),alpha=F))
dodgerblue2<-as.vector(col2rgb(col=c("dodgerblue2"),alpha=F))
dodgerblue3<-as.vector(col2rgb(col=c("dodgerblue3"),alpha=F))
dodgerblue4<-as.vector(col2rgb(col=c("dodgerblue4"),alpha=F))
firebrick4<-as.vector(col2rgb(col=c("firebrick4"),alpha=F))
firebrick3<-as.vector(col2rgb(col=c("firebrick3"),alpha=F))
firebrick2<-as.vector(col2rgb(col=c("firebrick2"),alpha=F))
firebrick1<-as.vector(col2rgb(col=c("firebrick1"),alpha=F))





# Euclidean distance between two points in space
EuclidDist<-function(x1,y1,z1,x2,y2,z2){
  
  dist<-((((x1-x2)^2)+((y1-y2)^2)+((z1-z2)^2))^(1/2))
  return(dist)
  
}


# Function to replace NAs with last value - from StackOverflow
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


# Moving mean, median, sum, min, max or sd
# returns a df with center index of the bins & value
# if new.length != NULL, directly create a vector that can be added to a df instead of a separate output df
MovingFUN<-function(x,winsize,shift,method,FUN,new.length=NULL){
  
  if(exists("repeat.before")==F){ #if repeat.before() is not in the environment
    
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
    
    
    message("Function 'repeat.before' had to be compiled.")
    
  }
  
  
  if(anyNA(x)==TRUE){
    x<-repeat.before(x)
  }
  
  
  fun.out<-c()
  index<-c()
  
  
  if(method=="centered"){
    i=1
    j=1
    while(i<(length(x)-(winsize/2))){
      
      if(FUN=="mean"){
        fun.out[j]<-mean(x[i:(i+(winsize-1))])  
      } else if(FUN=="sum"){
        fun.out[j]<-sum(x[i:(i+(winsize-1))])  
      } else if(FUN=="min"){
        fun.out[j]<-min(x[i:(i+(winsize-1))])
      } else if(FUN=="max"){
        fun.out[j]<-max(x[i:(i+(winsize-1))])
      } else if(FUN=="median"){
        fun.out[j]<-median(x[i:(i+(winsize-1))])
      } else if(FUN=="sd"){
        fun.out[j]<-sd(x[i:(i+(winsize-1))])
      }
      
      if(is.nan(fun.out[j])==T){ fun.out[j]<-NA }
      
      index[j]<-(i+(winsize/2))  
      i=i+shift
      j=j+1
      
    }
    
  } else if(method=="forward"){
    
    i=1
    j=1
    while(i<(length(x)-(winsize))){
      
      if(FUN=="mean"){
        fun.out[j]<-mean(x[i:(i+(winsize-1))])  
      } else if(FUN=="sum"){
        fun.out[j]<-sum(x[i:(i+(winsize-1))])  
      } else if(FUN=="min"){
        fun.out[j]<-min(x[i:(i+(winsize-1))])
      } else if(FUN=="max"){
        fun.out[j]<-max(x[i:(i+(winsize-1))])
      } else if(FUN=="median"){
        fun.out[j]<-median(x[i:(i+(winsize-1))])
      } else if(FUN=="sd"){
        fun.out[j]<-sd(x[i:(i+(winsize-1))])
      }
      
      if(is.nan(fun.out[j])==T){ fun.out[j]<-NA }
      
      index[j]<-(i+winsize)  
      i=i+shift
      j=j+1
      
    }
    
  } else if(method=="back"){
    
    i=winsize
    j=1
    while(i<=length(x)){
      
      if(FUN=="mean"){
        fun.out[j]<-mean(x[(i-winsize):i])  
      } else if(FUN=="sum"){
        fun.out[j]<-sum(x[(i-winsize):i])  
      } else if(FUN=="min"){
        fun.out[j]<-min(x[(i-winsize):i])
      } else if(FUN=="max"){
        fun.out[j]<-max(x[(i-winsize):i])
      } else if(FUN=="median"){
        fun.out[j]<-median(x[(i-winsize):i])
      } else if(FUN=="sd"){
        fun.out[j]<-sd(x[(i-winsize):i])
      }
      
      if(is.nan(fun.out[j])==T){ fun.out[j]<-NA }
      
      index[j]<-(i+winsize/2)  
      i=i+shift
      j=j+1
      
    }
    
  }
  
  
  
  df<-data.frame(index,fun.out) # make df
  
  
  if(is.null(new.length)==T){ #if we want a separate df
    
    return(df)
    
  } else { #if we want to make an appropriate vector
    
    ind<-df$index[is.na(df$fun.out)==F]  #indexes of the non-missing values
    vect<-rep(NA,new.length)  #create empty vector of the length new.length 
    vect[ind]<-as.numeric(na.omit(df$fun.out)) #put the non-missing values at the right index in the new vector
    # vect<-repeat.before(vect) 
    
    return(vect)
    
  }
  
  
}


# Mean with na.rm=TRUE
Mean<-function(x){
  base::mean(x,na.rm=TRUE,trim=0)
}


# Function that gives back the locations of the peaks in a time series 
# gives a vector with sequence of 0, +1  and -1 (+1 = positive peak, -1 = negative peak)
# time = time vector, var = vector on which to detect peaks
# halfwin = the number of frames in each direction of the peak of interest to calculate the mean and define a positive or negative peak.
# thresh = the minimum difference between the peak of interest and the surrounding mean to count it as a real peak
peaks<-function(time,var,halfwin,thresh){
  
  peaklocations<-rep(0,length(time))
  deriv<-diff(var)/diff(time)
  
  for (i in 1:(length(time)-1)){
    
    a<-deriv[i]
    b<-deriv[i+1]
    ratio<-a/b
    
    if (isTRUE(ratio<0)==TRUE){
      peaklocations[i+1]<-1
    }
    
  }
  
  
  # threshold to avoid detecting noise and classify the peaks as positive or negative
  d<-which(peaklocations==1)
  
  for(i in min(which(d>=halfwin)):length(d)){ #starting from the first peak where peaklocation-halfwin exists (>=0)
    
    m<-mean(var[(d[i]-halfwin):(d[i]+halfwin)]) # surrounding average (+- halfwin)
    
    if(isTRUE(m<var[d[i]])==TRUE #if the surrounding average is below the point of interest
       & isTRUE(abs(max(c(m,var[d[i]]))-min(c(m,var[d[i]])))>=thresh)==TRUE){  #AND the difference is above threshold
      
      peaklocations[d[i]]<-1 #positive peak
      
      
    } else if(isTRUE(m>var[d[i]])==TRUE #if the surrounding average is above the point of interest
              & isTRUE(abs(max(c(m,var[d[i]]))-min(c(m,var[d[i]])))>=thresh)==TRUE){ #AND the difference is above threshold
      
      peaklocations[d[i]]<--1 #negative peak
      
      
    } else {
      
      peaklocations[d[i]]<-0 #no peak
      
    }
    
  }
  
  
  return (as.numeric(peaklocations))
  
}


# Function to calculate the sweep frequency instananeously, in bins or sliding windows 
# gives back a vector with frequency values for each frame
# for instantaneous frequency, the first and last frame are counted as events.
# Arguments: peakvect (vector with 0s and 1s for "peak" or "no peak), method ("sliding window", "bins" or "instantaneous"), winsize
rate.calc<-function(peakvect,method,winsize){
  
  if (method=="sliding window"){
    
    timeunit<-winsize*0.01
    fsweeps<-c(rep(NA,length(peakvect)))
    start<-(winsize/2)+1
    end<-length(peakvect)-(winsize/2)
    
    for (i in start:end){
      
      win<-peakvect[(i-(winsize/2)):(i+((winsize/2)-1))]
      nsweeps<-(sum(win==1)-1)
      fsweeps[i]<-nsweeps/timeunit
      
    }
    
    fsweeps[which(fsweeps<0)]<-0
    return(fsweeps)
    
    
    
    
    
  } else if (method=="bins"){
    
    timeunit<-winsize*0.01
    fsweeps<-c(rep(NA,length(peakvect)))
    
    i<-1
    while (i<=length(peakvect)){
      
      win<-peakvect[i:(i+winsize)]
      nsweeps<-(sum(win==1)-1)
      fsweeps[i:(i+winsize)]<-nsweeps/timeunit
      i<-i+winsize+1
      
    }
    
    fsweeps[which(fsweeps<0)]<-0
    return(fsweeps)
    
    
    
    
    
    
  } else if (method=="instantaneous"){
    
    b<-which(peakvect==1)
    b<-c(1,b,length(peakvect)) #adding the first frame and the last frame
    c<-c(rep(NA,length(peakvect)))
    
    i<-1
    j<-b[1]
    while ((i<length(b))&(j<length(c[0:max(b)]))){
      
      d<-b[i+1]-b[i]
      c[(j+1):(j+d)]<-d
      i<-i+1
      j<-j+d
      
    }
    
    fsweeps<-(100/c)
    
    if(b[1]>1){
      fsweeps[1:b[1]]<-0
    }
    
    if(b[length(b)]<length(fsweeps)){
      fsweeps[b[length(b)]:length(fsweeps)]<-0
    }
    
    
    return(fsweeps)
    
  }
}


# Function to measure the amplitude of every sweep or in a sliding window or in bins
# Returns a vector of the length of the angle vector
# Arguments: peakvect, anglevect, method ("sweeps", "sliding", "bins"), winsize
sweep.amplitude<-function(method,peakvect,anglevect,winsize){
  
  if (method=="sweeps"){
    
    a<-which(peakvect==1)
    amp<-rep(NA,length(peakvect))
    
    i<-1
    j<-a[1]
    while ((i<=length(a))&(j<length(amp[0:max(a)]))){
      
      c<-abs((max(anglevect[a[i]:a[i+1]]))-(min(anglevect[a[i]:a[i+1]])))
      d<-(a[i+1]-a[i])
      amp[(j):(j+d-1)]<-c
      
      i<-i+1
      j<-j+d
      
    }
    
    return(amp)
    
  } else if (method=="sliding"){
    
    amp<-c(rep(NA,length(anglevect)))
    start<-(winsize/2)+1
    end<-length(anglevect)-(winsize/2)
    
    for (i in start:end){
      
      win<-anglevect[(i-(winsize/2)):(i+((winsize/2)-1))]
      amp[i]<-(max(win)-min(win))
      
    }
    
    return(amp)
    
  } else if (method=="bins"){
    
    amp<-c(rep(NA,length(anglevect)))
    j<-1
    k<-1
    
    while ((j<=(length(amp)-winsize)) & (k<=(length(anglevect)-winsize))){
      
      amp[j:(j+winsize)]<-(abs(max(anglevect[k:(k+winsize)]))-abs(min(anglevect[k:(k+winsize)])))
      k<-k+winsize
      j<-j+winsize
      
    }
    
    return(amp)
    
  }
  
}


# Function for the opposite of "%in%".
# Returns TRUE if the tested element IS NOT present in the comparison vector
"%!in%" <- function(x,y){
  !("%in%"(x,y))
}


# Function to correct unrealistic tracking jumps (both directions) and apply a low-pass filter if needed.
# Needs the package "signal" to be installed
# vector= coordinate vector (needs to be free of NAs)
# jumps= apply the jump correction?
# thresh= critical second derivative value to detect a jump (10-20 seems like a good start). We are detecting when the derivative changes too drastically.
# win= number of frames that will be replaced with the last value before the jump (doesn't take into account the actual length of the jump)
# butter.low= apply the Butterworth low-pass filter?
# cutoff= cut-off frequency for the filter as a fraction [0:1] of the Nyquist Frequency (= 0.5*samprate). For samprate=100Hz, NF=50HZ.
# Can give very different results with different thresholds. See case-by-case with AbsDeriv2().
trackingjumps<-function(vector,jumps,thresh,win,butter.low,cutoff){ 
  
  vector<-as.numeric(vector)
  
  na=FALSE #default
  if(anyNA(vector)==TRUE){  na=TRUE  }
  
  # jump correction
  if(jumps==TRUE){
    
    message('Applying jump correction.')
    if(na==TRUE){  message('WARNING: NAs detected. Might cause errors during filtering if butter.low = TRUE.')  }
    
    deriv<-diff(vector)/diff(seq(1,length(vector),1)) # first derivative 
    deriv<-abs(diff(deriv)/diff(seq(1,length(deriv),1)))  # absolute second derivative (considers both upwards and downwards jumps)
    
    jumpindexes<-which(deriv>=thresh) 
    #indexes where the derivative is above thresh (includes only the starts of the jumps)
    
    if(length(jumpindexes)!=0){
      
      for(i in 1:length(jumpindexes)){ #looping through the selected indexes
        ind<-jumpindexes[i]
        
        if(i==1){
          
          vector[ind:(ind+win)]<-vector[ind]
          #replacing the values from the start of the jump until start+win with the start value
          
        } else if((ind %in% jumpindexes[i-1]:(jumpindexes[i-1]+win))==FALSE){ 
          #if the current index is not within the last replacement window (and if it's not the first index of the list)
          #this avoids counting a jump twice and replacing too much or having to apply the correction twice.
          
          vector[ind:(ind+win)]<-vector[ind]
          
        } 
        
      }
      
    } else {  message('WARNING: No jumps detected with current threshold.')  }
    
  } else {  message('No jump correction.')  }
  
  # low-pass filtering
  if(butter.low==TRUE){
    
    message('Applying low-pass filter.')
    if(na==TRUE){  stop("NAs detected. Can't apply Butterworth filter.")  }
    message('   Note°1 : Filtering induces a significant edge effect! Take this into account in the data analysis!')
    message('   Note°2 : If the filter is applied to only a SUBSET of the original data, the filtered result will not match the original data.')
    message('   Note°3 : The filter is applied from both sides (i.e. twice in a row) to avoid a shift in the data. This will reduce the amplitude quite a bit. Take this into account when choosing a cut-off.')
    
    require("signal")
    
    lpfilter<-butter(1,W=cutoff,type="low") # W is a fraction of the Nyquist Frequency. 
    # The NF is the highest frequency that can be represented at a given sampling rate. NF= 0.5*samprate =50Hz  in our case.
    vector<-as.numeric(signal::filter(lpfilter,vector)) 
    vector<-as.numeric(signal::filter(lpfilter,rev(vector))) # applying filter from the other side too by flipping the vector
    vector<-rev(vector)# flipping it back into normal time order.
    
  } else {  message('No low-pass filter.')  }
  
  return(vector)
  
}

# absolute second order derivative to choose a threshold
AbsDeriv2<-function(vector){
  
  vector<-as.numeric(vector)
  
  deriv<-diff(vector)/diff(seq(1,length(vector),1)) # first derivative 
  deriv<-abs(diff(deriv)/diff(seq(1,length(deriv),1)))  # absolute second derivative (considers both upwards and downwards jumps)
  
  return(deriv)
}


# Function to find local maxima in any vector
# creates a vector with 0s and 1s = local maximum
localmax<-function(vect){
  
  index<-seq(1,length(vect),1)
  maxlocations<-c(rep(0,length(index)))
  deriv<-diff(vect)/diff(index)
  for (i in 1:(length(index)-1)){
    a<-deriv[i]
    b<-deriv[i+1]
    ratio<-a/b
    if ((isTRUE(ratio<=0)==TRUE)&(a>=0)){
      maxlocations[i+1]<-1
    }
  }
  
  # adding a threshold to avoid counting noise as peaks
  # d<-which(maxlocations==1)
  # for(i in 2:length(d)){
  #   if((isTRUE(abs(abs(vect[d[i]])-abs(vect[d[i-1]]))>=1)==TRUE)){
  #     maxlocations[d[i]]<-1
  #   } else {
  #     maxlocations[d[i]]<-0
  #   }
  # }
  
  return (as.vector(maxlocations))
}


# Function for the angle between the antenna and the midline
# parameters: antenna end points XYZ, reference point XYZ for the midline 
# reference can be a single point or a vector of same length as the antenna coordinates
# returns a vector of the same length as antenna coordinates vector with the angles to the reference midline
# Two versions: one for when Y increases towards the front, one for when it goes the other way (depends on the tracking system)
ant.angle.up<-function(plane,antX,antY,antZ,refX,refY,refZ){
  
  if(plane=="horizontal"){
    
    if((length(refX)==length(refY)) & 
       (((length(refX)==1) & (length(refY)==1)) | ((length(refX)==length(antX)) & (length(refY)==length(antX))))){
      #if refX and Y have the same length AND if they are either of length 1 or the length of the antenna coordinates
      
      if((length(refX)==1) | (length(refY)==1)){
        refX<-rep(refX,length(antX)) # making vectors with the point coordinates for the reference
        refY<-rep(refY,length(antY))
      }
      
      dev<-((((refX-antX)^2)+((antY-antY)^2))^(1/2)) 
      #deviation = euclidean distance between the antenna end (antX & antY) and its projection on the midline from the reference (refX & refY)
      ant.length<-((((refX-antX)^2)+((refY-antY)^2))^(1/2)) 
      #antenna length = euclidean distance between the antenna end (antX & antY) and the reference point (refX & refY)
      
      angle<-asin(dev/ant.length) #absolute angle (radians)
      angle<-((angle*180)/pi) #change to degrees
      
      for(i in 1:length(angle)){
        if(isTRUE(refX[i]<antX[i])==TRUE){
          angle[i]<-angle[i] #positive angle if the antenna is to the right of the reference midline
          if(isTRUE(antY[i]<refY[i])==TRUE){
            angle[i]<-(90+abs(90-angle[i])) #if the antenna tip is behind the head, the angle is more than 90°
          }
        } else if(isTRUE(refX[i]>antX[i])==TRUE){
          angle[i]<-(-angle[i]) #negative angle if it's to the left 
          if(isTRUE(antY[i]<refY[i])==TRUE){
            angle[i]<-((-90)-abs(90-abs(angle[i]))) #if the antenna tip is behind the head, the angle is more than 90°in the negative direction
          }
        }
      }
      
      return(angle)
      
    } else { 
      stop("Reference vectors have to be either of length 1 or as long as the antenna coordinates vector")
    }
    
  } 
  else if(plane=="vertical"){
    
    if((length(refY)==length(refZ)) & (length(refZ)==length(refX)) & (((length(refY)==1) & (length(refZ)==1)) | 
                                                                      ((length(refY)==length(antY)) & (length(refZ)==length(antZ)) & (length(refX)==length(antX))))){
      
      if((length(refY)==1) | (length(refZ)==1) | (length(refX)==1)){
        refX<-rep(refX,length(antX))
        refY<-rep(refY,length(antY)) 
        refZ<-rep(refZ,length(antZ))
      }
      
      dev<-(((antX-antX)^2)+(((antY-antY)^2)+((refZ-antZ)^2))^(1/2)) #deviation from the horizontal plane at the level of the head 
      ant.length<-(((refX-antX)^2)+(((refY-antY)^2)+((refZ-antZ)^2))^(1/2)) # need length in 3D
      
      angle<-asin(dev/ant.length)
      angle<-((angle*180)/pi)  #degrees
      
      for(i in 1:length(angle)){
        if(isTRUE(antZ[i]<refZ[i])==TRUE){ # and if the tip is lower than the head
          angle[i]<-(-angle[i])  #negative angle
        }
      }
      
      return(angle)
      
    } else { 
      stop("Reference vectors have to be either of length 1 or as long as the antenna coordinates vector")
    }
    
  }
  
}











#--------------------------------------------------------------------------------------------------------------------------
# Run first: Getting raw stepping rate and heading direction 
#--------------------------------------------------------------------------------------------------------------------------



#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)





##########################################################
#   Getting raw head direction and stepping rate first   #
##########################################################

# run this first to generate the data necessary for the subsequent analyses

data<-AllData
datlist.3D<-split(data,f=data$Test)


#calculating step frequency for the front and middle legs
#scaling Freq for X and Y within each trial and leg
#calculating heading angle 

# X (sideways)
datlist.3D<-lapply(datlist.3D,
                   FUN=function(df){
                     
                     print(unique(df$Test))
                     
                     
                     #Lfrontleg
                     print("Lfrontleg X")
                     
                     if(sum(is.na(df$Lfrontleg.X))==nrow(df)){ #if the column only contains NAs
                       
                       df$LFrontSteps.X<-rep(0,nrow(df))
                       df$LFrontFreq.X<-rep(0,nrow(df))
                       
                     } else {
                       
                       df$Lfrontleg.X<-repeat.before(df$Lfrontleg.X)
                       a<-as.numeric(trackingjumps(df$Lfrontleg.X,jumps=T,thresh=5,win=20,butter.low=T,cutoff=0.1))
                       
                       if(length(a)>=nrow(df)){
                         
                         df$LFrontSteps.X<-a[1:nrow(df)] #in some cases, the filtering adds frames. I don't understand why, but it's minor
                       
                      } else {
                         
                         df$LFrontSteps.X<-c(a,rep(0,(nrow(df)-length(a))))
                         
                       }
                       
                       df$LFrontSteps.X<-peaks(time=df$Count,var=df$LFrontSteps.X,thresh=0.25,halfwin=15) 
                       df$LFrontSteps.X[df$LFrontSteps.X!=1]<-0 #only keep positive peaks
                       
                       if(sum(df$LFrontSteps.X)>=2){ #if there are at least two peaks (can't calculate anything with just 0 or 1)
                         
                         df$LFrontFreq.X<-rate.calc(df$LFrontSteps.X,method="instantaneous") #freq per step
                         df$LFrontFreq.X[df$LFrontFreq.X>10]<-NA #remove frequency values above 10, because that is too fast for legs and is due to tracking errors
                         df$LFrontFreq.X<-repeat.before(df$LFrontFreq.X) #replace with previous frequency value
                         df$LFrontFreq.X<-df$LFrontFreq.X/max(df$LFrontFreq.X,na.rm=T) #scaling
                         
                       } else {
                         
                         df$LFrontFreq.X<-rep(0,nrow(df)) 
                         
                       }
                       
                     }
                     
                     
                     
                     #Rfrontleg
                     print("Rfrontleg X")
                     if(sum(is.na(df$Rfrontleg.X))==nrow(df)){ #if the column only contains NAs
                       
                       df$RFrontSteps.X<-rep(0,nrow(df))
                       df$RFrontFreq.X<-rep(0,nrow(df))
                       
                     } else {
                       
                       df$Rfrontleg.X<-repeat.before(df$Rfrontleg.X)
                       a<-as.numeric(trackingjumps(df$Rfrontleg.X,jumps=T,thresh=5,win=20,butter.low=T,cutoff=0.1))
                       
                       if(length(a)>=nrow(df)){
                         
                         df$RFrontSteps.X<-a[1:nrow(df)] 
                         
                       } else {
                         
                         df$RFrontSteps.X<-c(a,rep(0,(nrow(df)-length(a))))
                         
                       }
                       
                       df$RFrontSteps.X<-peaks(time=df$Count,var=df$RFrontSteps.X,thresh=0.25,halfwin=15) 
                       df$RFrontSteps.X[df$RFrontSteps.X!=1]<-0 
                       
                       if(sum(df$RFrontSteps.X)>=2){
                         
                         df$RFrontFreq.X<-rate.calc(df$RFrontSteps.X,method="instantaneous") 
                         df$RFrontFreq.X[df$RFrontFreq.X>10]<-NA
                         df$RFrontFreq.X<-repeat.before(df$RFrontFreq.X)
                         df$RFrontFreq.X<-df$RFrontFreq.X/max(df$RFrontFreq.X,na.rm=T)
                         
                       } else {
                         
                         df$RFrontFreq.X<-rep(0,nrow(df)) 
                         
                       }
                       
                     }
                     
                     
                     
                     #Lmidleg
                     print("Lmidleg X")
                     if(sum(is.na(df$Lmidleg.X))==nrow(df)){ #if the column only contains NAs
                       
                       df$LMidSteps.X<-rep(0,nrow(df))
                       df$LMidFreq.X<-rep(0,nrow(df))
                       
                     } else {
                       
                       df$Lmidleg.X<-repeat.before(df$Lmidleg.X)
                       a<-as.numeric(trackingjumps(df$Lmidleg.X,jumps=T,thresh=5,win=20,butter.low=T,cutoff=0.1))
                       
                       if(length(a)>=nrow(df)){
                         
                         df$LMidSteps.X<-a[1:nrow(df)] 
                         
                       } else {
                         
                         df$LMidSteps.X<-c(a,rep(0,(nrow(df)-length(a))))
                         
                       }
                       
                       df$LMidSteps.X<-peaks(time=df$Count,var=df$LMidSteps.X,thresh=0.25,halfwin=15) 
                       df$LMidSteps.X[df$LMidSteps.X!=1]<-0 
                       
                       if(sum(df$LMidSteps.X)>=2){
                         
                         df$LMidFreq.X<-rate.calc(df$LMidSteps.X,method="instantaneous") 
                         df$LMidFreq.X[df$LMidFreq.X>10]<-NA
                         df$LMidFreq.X<-repeat.before(df$LMidFreq.X)
                         df$LMidFreq.X<-df$LMidFreq.X/max(df$LMidFreq.X,na.rm=T)
                         
                       } else {
                         
                         df$LMidFreq.X<-rep(0,nrow(df)) 
                         
                       }
                       
                     }
                     
                     
                     
                     #Rmidleg
                     print("Rmidleg X")
                     if(sum(is.na(df$Rmidleg.X))==nrow(df)){ #if the column only contains NAs
                       
                       df$RMidSteps.X<-rep(0,nrow(df))
                       df$RMidFreq.X<-rep(0,nrow(df))
                       
                     } else {
                       
                       df$Rmidleg.X<-repeat.before(df$Rmidleg.X)
                       a<-as.numeric(trackingjumps(df$Rmidleg.X,jumps=T,thresh=5,win=20,butter.low=T,cutoff=0.1))
                       
                       if(length(a)>=nrow(df)){
                         
                         df$RMidSteps.X<-a[1:nrow(df)] 
                         
                       } else {
                         
                         df$RMidSteps.X<-c(a,rep(0,(nrow(df)-length(a))))
                         
                       }
                       
                       df$RMidSteps.X<-peaks(time=df$Count,var=df$RMidSteps.X,thresh=0.25,halfwin=15) 
                       df$RMidSteps.X[df$RMidSteps.X!=1]<-0 
                       
                       if(sum(df$RMidSteps)>=2){
                         
                         df$RMidFreq.X<-rate.calc(df$RMidSteps.X,method="instantaneous") 
                         df$RMidFreq.X[df$RMidFreq.X>10]<-NA
                         df$RMidFreq.X<-repeat.before(df$RMidFreq.X)
                         df$RMidFreq.X<-df$RMidFreq.X/max(df$RMidFreq.X,na.rm=T)
                         
                       } else {
                         
                         df$RMidFreq.X<-rep(0,nrow(df)) 
                         
                       }
                       
                     }
                     
                     
                     
                     
                     #heading angle
                     print("head")
                     
                     a<-ant.angle.up(plane="horizontal",antX=df$head.X,antY=df$head.Y,refX=df$tether.X,refY=df$tether.Y) # using the tether point as a reference
                     a<-as.numeric(trackingjumps(a,jumps=T,thresh=5,win=10,butter.low=T,cutoff=0.02)) #smooth
                     
                     if(length(a)>=nrow(df)){
                       
                       df$heading<-a[1:nrow(df)] 
                       
                     } else {
                       
                       df$heading<-c(a,rep(0,(nrow(df)-length(a))))
                       
                     }
                     
                     b<-df$head.X
                     c<-df$tether.X
                     b<-b-c
                     b<-as.numeric(trackingjumps(df$head.X,jumps=T,thresh=5,win=10,butter.low=T,cutoff=0.02))

                     if(length(b)>=nrow(df)){
                       
                       df$head.X.centersmooth<-b[1:nrow(df)] 
                       
                     } else {
                       
                       df$head.X.centersmooth<-c(a,rep(0,(nrow(df)-length(b))))
                       
                     }
                     
                     return(df)
                     
                   })


# Y (forward)
datlist.3D<-lapply(datlist.3D,
                   FUN=function(df){
                     
                     print(unique(df$Test))
                     
                     
                     #Lfrontleg
                     print("Lfrontleg Y")
                     
                     if(sum(is.na(df$Rfrontleg.Y))==nrow(df)){ #if the column only contains NAs
                       
                       df$RFrontSteps.Y<-rep(0,nrow(df))
                       df$RFrontFreq.Y<-rep(0,nrow(df))
                       
                     } else {
                       
                       df$Lfrontleg.Y<-repeat.before(df$Lfrontleg.Y)
                       a<-as.numeric(trackingjumps(df$Lfrontleg.Y,jumps=T,thresh=5,win=20,butter.low=T,cutoff=0.1))
                       
                       if(length(a)>=nrow(df)){
                         
                         df$LFrontSteps.Y<-a[1:nrow(df)]  #in some cases, the filtering adds frames. I don't understand why, but it's minor
                         
                       } else {
                         
                         df$LFrontSteps.Y<-c(a,rep(0,(nrow(df)-length(a))))
                         
                       }
                       
                       df$LFrontSteps.Y<-peaks(time=df$Count,var=df$LFrontSteps.Y,thresh=0.25,halfwin=15) 
                       df$LFrontSteps.Y[df$LFrontSteps.Y!=1]<-0 #only keep positive peaks
                       
                       if(sum(df$LFrontSteps.Y)>=2){ #if there are at least two peaks (can't calculate anything with just 0 or 1)
                         
                         df$LFrontFreq.Y<-rate.calc(df$LFrontSteps.Y,method="instantaneous") #freq per step
                         df$LFrontFreq.Y[df$LFrontFreq.Y>10]<-NA #remove frequency values above 10, because that is too fast for legs and is due to tracking errors
                         df$LFrontFreq.Y<-repeat.before(df$LFrontFreq.Y) #replace with previous frequency value
                         df$LFrontFreq.Y<-df$LFrontFreq.Y/max(df$LFrontFreq.Y,na.rm=T) #scaling
                         
                       } else {
                         
                         df$LFrontFreq.Y<-rep(0,nrow(df)) 
                         
                       }
                       
                     }
                     
                     
                     
                     #Rfrontleg
                     print("Rfrontleg Y")
                     if(sum(is.na(df$Rfrontleg.Y))==nrow(df)){ #if the column only contains NAs
                       
                       df$RFrontSteps.Y<-rep(0,nrow(df))
                       df$RFrontFreq.Y<-rep(0,nrow(df))
                       
                     } else {
                       
                       df$Rfrontleg.Y<-repeat.before(df$Rfrontleg.Y)
                       a<-as.numeric(trackingjumps(df$Rfrontleg.Y,jumps=T,thresh=5,win=20,butter.low=T,cutoff=0.1))
                       
                       if(length(a)>=nrow(df)){
                         
                         df$RFrontSteps.Y<-a[1:nrow(df)] 
                         
                       } else {
                         
                         df$RFrontSteps.Y<-c(a,rep(0,(nrow(df)-length(a))))
                         
                       }
                       
                       df$RFrontSteps.Y<-peaks(time=df$Count,var=df$RFrontSteps.Y,thresh=0.25,halfwin=15) 
                       df$RFrontSteps.Y[df$RFrontSteps.Y!=1]<-0 
                       
                       if(sum(df$RFrontSteps.Y)>=2){
                         
                         df$RFrontFreq.Y<-rate.calc(df$RFrontSteps.Y,method="instantaneous") 
                         df$RFrontFreq.Y[df$RFrontFreq.Y>10]<-NA
                         df$RFrontFreq.Y<-repeat.before(df$RFrontFreq.Y)
                         df$RFrontFreq.Y<-df$RFrontFreq.Y/max(df$RFrontFreq.Y,na.rm=T)
                         
                       } else {
                         
                         df$RFrontFreq.Y<-rep(0,nrow(df)) 
                         
                       }
                       
                     }
                     
                     
                     
                     #Lmidleg
                     print("Lmidleg Y")
                     if(sum(is.na(df$Lmidleg.Y))==nrow(df)){ #if the column only contains NAs
                       
                       df$LMidSteps.Y<-rep(0,nrow(df))
                       df$LMidFreq.Y<-rep(0,nrow(df))
                       
                     } else {
                       
                       df$Lmidleg.Y<-repeat.before(df$Lmidleg.Y)
                       a<-as.numeric(trackingjumps(df$Lmidleg.Y,jumps=T,thresh=5,win=20,butter.low=T,cutoff=0.1))
                       
                       if(length(a)>=nrow(df)){
                         
                         df$LMidSteps.Y<-a[1:nrow(df)] 
                         
                       } else {
                         
                         df$LMidSteps.Y<-c(a,rep(0,(nrow(df)-length(a))))
                         
                       }
                       
                       df$LMidSteps.Y<-peaks(time=df$Count,var=df$LMidSteps.Y,thresh=0.25,halfwin=15) 
                       df$LMidSteps.Y[df$LMidSteps.Y!=1]<-0 
                       
                       if(sum(df$LMidSteps.Y)>=2){
                         
                         df$LMidFreq.Y<-rate.calc(df$LMidSteps.Y,method="instantaneous") 
                         df$LMidFreq.Y[df$LMidFreq.Y>10]<-NA
                         df$LMidFreq.Y<-repeat.before(df$LMidFreq.Y)
                         df$LMidFreq.Y<-df$LMidFreq.Y/max(df$LMidFreq.Y,na.rm=T)
                         
                       } else {
                         
                         df$LMidFreq.Y<-rep(0,nrow(df)) 
                         
                       }
                       
                     }
                     
                     
                     
                     #Rmidleg
                     print("Rmidleg Y")
                     if(sum(is.na(df$Rmidleg.Y))==nrow(df)){ #if the column only contains NAs
                       
                       df$RMidSteps.Y<-rep(0,nrow(df))
                       df$RMidFreq.Y<-rep(0,nrow(df))
                       
                     } else {
                       
                       df$Rmidleg.Y<-repeat.before(df$Rmidleg.Y)
                       a<-as.numeric(trackingjumps(df$Rmidleg.Y,jumps=T,thresh=5,win=20,butter.low=T,cutoff=0.1))
                       
                       if(length(a)>=nrow(df)){
                         
                         df$RMidSteps.Y<-a[1:nrow(df)] 
                         
                       } else {
                         
                         df$RMidSteps.Y<-c(a,rep(0,(nrow(df)-length(a))))
                         
                       }
                       
                       df$RMidSteps.Y<-peaks(time=df$Count,var=df$RMidSteps.Y,thresh=0.25,halfwin=15) 
                       df$RMidSteps.Y[df$RMidSteps.Y!=1]<-0 
                       
                       if(sum(df$RMidSteps.Y)>=2){
                         
                         df$RMidFreq.Y<-rate.calc(df$RMidSteps.Y,method="instantaneous") 
                         df$RMidFreq.Y[df$RMidFreq.Y>10]<-NA
                         df$RMidFreq.Y<-repeat.before(df$RMidFreq.Y)
                         df$RMidFreq.Y<-df$RMidFreq.Y/max(df$RMidFreq.Y,na.rm=T)
                         
                       } else {
                         
                         df$RMidFreq.Y<-rep(0,nrow(df)) 
                         
                       }
                       
                     }
                     
                     return(df)
                     
                   })



# trial avg of the 4 front legs (most reliable tracking and it's enough to estimate walking behavior)
datlist.3D<-lapply(datlist.3D,
                   FUN=function(df){ 
                     
                     df2<-df[,c("LFrontFreq.Y","RFrontFreq.Y","LMidFreq.Y","RMidFreq.Y",
                                "LFrontFreq.X","RFrontFreq.X","LMidFreq.X","RMidFreq.X")]
                     
                     
                     Freqavg<-df2[,c("LFrontFreq.Y","RFrontFreq.Y","LMidFreq.Y","RMidFreq.Y")]
                     Freqavg<-apply(Freqavg,1,mean,na.rm=T)
                     Freqavg[is.nan(Freqavg)]<-0 # freq change at 0 if NaN
                     
                     df$LegsFreqavg<-Freqavg
                     
                     
                     return(df)
                     
                   })


data<-dplyr::bind_rows(datlist.3D)













#--------------------------------------------------------------------------------------------------------------------------
# Stepping
#--------------------------------------------------------------------------------------------------------------------------


# These trials had no walking activity at all and were excluded
non.walking<-c("T520","T522","T527","T530","T533","T564","T575","T617","T645","T662","T675","T676","T677","T678",
               "T679","T680","T683","T684","T686","T691","T700","T703","T705","T778","T780","T787","T821","T825",
               "T829","T842","T843","T866","T895","T899","T900","T903","T907","T926","T945","T946","T950","T951")






######################################
#   Stepping rate avg times-series   #
######################################


# pooled first 2 seconds of PB, PGb & PE 
stepping<-dplyr::filter(data,
                        Count>=-2 & Count<2,
                        Protocol %in% c("ProtocolB","ProtocolGb","ProtocolE"),
                        Test %!in% non.walking)


#centering
stepping<-split(stepping,f=stepping$Test)
stepping<-lapply(stepping,
                 FUN=function(df){ 
                   
                   psavg<-mean(df$LegsFreqavg[df$Count<0])
                   df$stepping<-df$LegsFreqavg-psavg 
                   
                   return(df)
                 })
stepping<-dplyr::bind_rows(stepping)



#avg
stepping.sem<-aggregate(stepping~Count+Odor,data=stepping,FUN="sd",na.rm=T)
stepping.n<-aggregate(stepping~Count+Odor,data=stepping,FUN="length")
stepping.sem$stepping<-stepping.sem$stepping/sqrt(stepping.n$stepping) #SEM

stepping<-aggregate(stepping~Count+Odor,data=stepping,FUN="mean",na.rm=T)
stepping$sem<-stepping.sem$stepping




#plot pooled times series
pdf(file=paste("./Plots/AvgSteppingRate_pooled_timeseries.pdf",sep=""),width=6,height=15,bg="white")

par(mfrow=c(3,1))
for(i in unique(stepping$Odor)){
  
  if(grepl("Blank",i)==TRUE) { next }
  
  dat<-dplyr::filter(stepping,Odor==i)
  
  plot(dat$Count,dat$stepping,
       type="l",lwd=2,col="black",xaxs="i",yaxs="i",
       ylim=c(-0.15,0.15),
       main=paste(i,"Stepping rate change"),ylab="Avg change [Hz]")
  
  polygon(x=c(dat$Count,rev(dat$Count)),
          y=c(dat$stepping+dat$sem,rev(dat$stepping-dat$sem)),
          col=rgb(0,0,0,alpha=100,names=NULL,maxColorValue=255),
          border=NA)
  
  abline(h=0,lty=2)
  abline(v=0)
  
}


dev.off()











###############################
#   Stepping rate modelling   #
###############################

# stepping<-data
stepping<-dplyr::filter(data, 
                        Test %!in% non.walking,
                        Protocol %in% c("ProtocolB","ProtocolGb","ProtocolE"),
                        Count>=-2,
)


stepping<-split(stepping,f=stepping$Test)
stepping<-lapply(stepping,
                 FUN=function(df){ 
                   
                   df$stepping<-df$LegsFreqavg-mean(df$LegsFreqavg[df$TestPhase=="PS"]) #center
                   
                   stepping.ps<-mean(df$stepping[df$TestPhase=="PS"],na.rm=T)
                   stepping.s<-mean(df$stepping[df$TestPhase=="S" & df$Count<2],na.rm=T)
                   
                   if(unique(unique(df$Protocol)=="ProtocolB")){
                     
                     stepping.psts<-mean(df$stepping[df$TestPhase=="PstS" & df$Count<8],na.rm=T)
                     
                   } else if(unique(df$Protocol)=="ProtocolGb"){
                     
                     stepping.psts<-mean(df$stepping[df$TestPhase=="PstS" & df$Count<6],na.rm=T)
                     
                   } else{
                     
                     stepping.psts<-NA
                     
                   }
                   
                   
                   df<-data.frame("Prot"=rep(unique(df$Protocol),3),
                                  "Odor"=rep(unique(df$Odor),3),
                                  "Test"=rep(unique(df$Test),3),
                                  "TestPhase"=c("PS","S","PstS"),
                                  "stepping"=c(stepping.ps,stepping.s,stepping.psts)
                   )
                   
                   return(df)
                 })
stepping<-dplyr::bind_rows(stepping)








### Modelling ###


modeldat<-dplyr::filter(stepping,
                        Odor !="Blank",
                        TestPhase %in% c("S","PstS"))


# retrieve animal ID
animal<-data[data$Test %in% modeldat$Test,c("Animal","Test")] #select correct animals based on test ids in dat
animal<-animal[!duplicated(animal),] #remove duplcated rows to have 1 iteration of Test-Animal pairs
animal<-rep(animal$Animal,each=2) #double each iteration because we have 2 windows per trial

modeldat$Animal<-animal #add animal ID to our model data


#view overall distribution
hist(modeldat$stepping)



set.seed(170822)
brmsmod<-brm(stepping ~ Odor + TestPhase + Odor:TestPhase + (1|Animal), 
             data=modeldat,
             family=student(),
             iter=10000, chains=4, cores=4,backend="cmdstanr", 
             control=list(adapt_delta=0.99,max_treedepth=20))





# pp check
pp<-pp_check(brmsmod,ndraws=200) 
pp


pp_check(brmsmod,type="stat",stat="mean",ndraws=200)
pp_check(brmsmod,type="stat",stat="max",ndraws=200)
pp_check(brmsmod,type="stat",stat="min",ndraws=200)



# convergence check
plot(brmsmod)

#quick effect plot
ce<-conditional_effects(brmsmod)
ce




### Check real data vs predicted data correlation

pred<-as.data.frame(predict(brmsmod)) #predicts data based on the model
dat<-na.omit(modeldat)
pred<-cbind(dat,pred)


# get means 
mean_reald<-aggregate(stepping~Odor,data=pred,"mean") 
mean_predd<-aggregate(Estimate~Odor,data=pred,"mean")
real_vs_pred_avg<-merge(mean_reald,mean_predd)

#plot correlation
plot(NA,xlim=c(-0.1,0.05),
     ylim=c(-0.2,0.2),
     ylab="pred",xlab="real") 


d<-dplyr::filter(pred,Odor=="Butanol")
points(x=rep(real_vs_pred_avg$stepping[1],nrow(d)),y=d$Estimate) # all data points for Butanol

d<-dplyr::filter(pred,Odor=="Linalool")
points(x=rep(real_vs_pred_avg$stepping[3],nrow(d)),y=d$Estimate) # all data points for Linalool

d<-dplyr::filter(pred,Odor=="Colony")
points(x=rep(real_vs_pred_avg$stepping[2],nrow(d)),y=d$Estimate) # all data points for Colony


points(Estimate~stepping,data=real_vs_pred_avg, # means
       pch=20,col="firebrick4",cex=2)


abline(a=0,b=1)










# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))






# Estimates, 95% CrI, and probabilities
ColS.quants<-round(quantile((simval[,1]+simval[,2]+simval[,4]+simval[,5]),probs=c(0.5,0.025,0.975)),2)
ColS.quants
print(paste("P(ColS)>0 = ",round(mean((simval[,1]+simval[,2]+simval[,4]+simval[,5])>0),2),sep=""))

ColPstS.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColPstS.quants
print(paste("P(ColPstS)>0 = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))


ButS.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
ButS.quants
print(paste("P(ButS)>0 = ",round(mean((simval[,1]+simval[,4])>0),2),sep=""))

ButPstS.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ButPstS.quants
print(paste("P(ButPstS)>0 = ",round(mean((simval[,1])>0),2),sep=""))


LinS.quants<-round(quantile((simval[,1]+simval[,3]+simval[,4]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
LinS.quants
print(paste("P(LinS)>0 = ",round(mean((simval[,1]+simval[,3]+simval[,4]+simval[,6])>0),2),sep=""))

LinPstS.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
LinPstS.quants
print(paste("P(LinPstS)>0 = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))








#effect plot
pdf(file=paste("./Plots/SteppingRateChange_pooled_2sWindows_EffectBarplot_allodors.pdf",sep=""),
    width=8,height=7,bg="white")


par(lwd=5)


m.fit<-matrix(ncol=3,nrow=2)
colnames(m.fit)<-c("Col","But","Lin")
rownames(m.fit)<-c("S","PstS")
m.fit[,1]<-c(ColS.quants[1],ColPstS.quants[1])
m.fit[,2]<-c(ButS.quants[1],ButPstS.quants[1])
m.fit[,3]<-c(LinS.quants[1],LinPstS.quants[1])

m.lwr<-matrix(ncol=3,nrow=2)
colnames(m.lwr)<-c("Col","But","Lin")
rownames(m.lwr)<-c("S","PstS")
m.lwr[,1]<-c(ColS.quants[2],ColPstS.quants[2])
m.lwr[,2]<-c(ButS.quants[2],ButPstS.quants[2])
m.lwr[,3]<-c(LinS.quants[2],LinPstS.quants[2])

m.upr<-matrix(ncol=3,nrow=2)
colnames(m.upr)<-c("Col","But","Lin")
rownames(m.upr)<-c("S","PstS")
m.upr[,1]<-c(ColS.quants[3],ColPstS.quants[3])
m.upr[,2]<-c(ButS.quants[3],ButPstS.quants[3])
m.upr[,3]<-c(LinS.quants[3],LinPstS.quants[3])


ylim=c(-0.15,0.15)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="[Hz]",main="Change in stepping rate")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.05,2),xpd=T,labels=rep(c("S","PstS"),3),cex=0.7)
text(x=barcenters[1,]+0.5,y=rep(par("usr")[3]-0.05,3),xpd=T,labels=rep(c("Col","But","Lin"),3),cex=1)



dev.off()
















#--------------------------------------------------------------------------------------------------------------------------
# Heading
#--------------------------------------------------------------------------------------------------------------------------



##########################################
#   Heading direction avg times-series   #
##########################################

# pooled first 2 seconds

heading<-dplyr::filter(data,
                       Count>-2 & Count<=2, #the first 2 seconds of the shifting odour are in the center as well.
                       Protocol %in% c("ProtocolB","ProtocolGb","ProtocolE")) #PB & PGb = center odour; PE = shifting odour


#absolute heading deviation for the center odor
heading$heading<-abs(heading$heading) 

#centering after abs() 
heading<-split(heading,f=heading$Test)
heading<-lapply(heading,
                FUN=function(df){
                  
                  psavg<-mean(df$heading[df$Count<0])
                  df$heading<-df$heading-psavg
                  
                  return(df)
                })
heading<-dplyr::bind_rows(heading)



#avg
heading.sem<-aggregate(heading~Count+Odor,data=heading,FUN="sd",na.rm=T)
heading.n<-aggregate(heading~Count+Odor,data=heading,FUN="length")
heading.sem$heading<-heading.sem$heading/sqrt(heading.n$heading) #SEM

heading<-aggregate(heading~Count+Odor,data=heading,FUN="mean",na.rm=T)

heading$sem<-heading.sem$heading






#plot pooled times series
pdf(file=paste("./Plots/AbsHeadingChange_pooledPBPGbPE_timeseries.pdf",sep=""),
    width=6,height=15,bg="white")


par(mfrow=c(3,1))
for(i in unique(heading$Odor)){
  
  if(grepl("Blank",i)==TRUE) { next }
  
  dat<-dplyr::filter(heading,Odor==i)
  
  plot(dat$Count,dat$heading,
       type="l",lwd=2,col="black",xaxs="i",yaxs="i",
       ylim=c(-1.5,2),
       xlim=c(-2,2),
       main=paste(i,"Abs heading change"),
ylab="heading change [°]")
  
  polygon(x=c(dat$Count,rev(dat$Count)),
          y=c(dat$heading+dat$sem,rev(dat$heading-dat$sem)),
          col=rgb(0,0,0,alpha=100,names=NULL,maxColorValue=255),
          border=NA)
  
  abline(h=0,lty=2)
  abline(v=0)
  
}


dev.off()









###################################
#   Heading direction modelling   #
###################################

# average in time bins

heading<-dplyr::filter(data,
                       Count>=-2, 
                       Protocol %in% c("ProtocolB","ProtocolGb","ProtocolE"))


heading<-split(heading,f=heading$Test)
heading<-lapply(heading,
                FUN=function(df){ 
                  
                  # print(unique(df$Test))
                  
                  df$heading<-abs(df$heading)
                  df$heading<-df$heading-mean(df$heading[df$Count<0]) #centering
                  
                  heading.ps<-mean(df$heading[df$TestPhase=="PS"],na.rm=T)
                  heading.s<-mean(df$heading[df$TestPhase=="S" & df$Count<2],na.rm=T)
                  
                  if(unique(df$Protocol)=="ProtocolB"){ #6s stim
                    
                    heading.psts<-mean(df$heading[df$TestPhase=="PstS" & df$Count<8],na.rm=T)
                    
                  } else if(unique(df$Protocol)=="ProtocolGb"){ #4s stim
                    
                    heading.psts<-mean(df$heading[df$TestPhase=="PstS" & df$Count<6],na.rm=T)
                    
                  } else{
                    
                    heading.psts<-NA  # the post-stimulus for the shifting experiment seems less comparable to the others, since it's far on the side now. so I ignore it
                    
                  }
                  
                  
                  df<-data.frame("Prot"=rep(unique(df$Protocol),3),
                                 "Odor"=rep(unique(df$Odor),3),
                                 "Test"=rep(unique(df$Test),3),
                                 "TestPhase"=c("PS","S","PstS"),
                                 "heading"=c(heading.ps,heading.s,heading.psts)
                  )
                  
                  
                  return(df)
                  
                })
heading<-dplyr::bind_rows(heading)








### Modelling ###

modeldat<-dplyr::filter(heading,
                   Odor !="Blank",
                   TestPhase %in% c("S","PstS"))


# retrieve animal ID
animal<-data[data$Test %in% modeldat$Test,c("Animal","Test")] #select correct animals based on test ids in dat
animal<-animal[!duplicated(animal),] #remove duplcated rows to have 1 iteration of Test-Animal pairs
animal<-rep(animal$Animal,each=2) #double each iteration because we have 2 windows per trial

modeldat$Animal<-animal #add animal ID to our model data


# look at overall distribution
hist(modeldat$heading)


# fit model
set.seed(170822)
brmsmod<-brm(heading ~ Odor + TestPhase + Odor:TestPhase + (1|Animal), 
             data=modeldat,
             family=student(),
             iter=10000, chains=4, cores=4,backend="cmdstanr", 
             control=list(adapt_delta=0.99,max_treedepth=20))





# pp check
pp<-pp_check(brmsmod,ndraws=200)
pp

pp_check(brmsmod,type="stat",stat="mean",ndraws=200)
pp_check(brmsmod,type="stat",stat="max",ndraws=200)
pp_check(brmsmod,type="stat",stat="min",ndraws=200)


#convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce


# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))


# Estimates, 95% CrI, and probabilities
ColS.quants<-round(quantile((simval[,1]+simval[,2]+simval[,4]+simval[,5]),probs=c(0.5,0.025,0.975)),2)
ColS.quants
print(paste("P(ColS)>0 = ",round(mean((simval[,1]+simval[,2]+simval[,4]+simval[,5])>0),2),sep=""))

ColPstS.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColPstS.quants
print(paste("P(ColPstS)>0 = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))


ButS.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
ButS.quants
print(paste("P(ButS)>0 = ",round(mean((simval[,1]+simval[,4])>0),2),sep=""))

ButPstS.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ButPstS.quants
print(paste("P(ButPstS)>0 = ",round(mean((simval[,1])>0),2),sep=""))


LinS.quants<-round(quantile((simval[,1]+simval[,3]+simval[,4]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
LinS.quants
print(paste("P(LinS)>0 = ",round(mean((simval[,1]+simval[,3]+simval[,4]+simval[,6])>0),2),sep=""))

LinPstS.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
LinPstS.quants
print(paste("P(LinPstS)>0 = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))





#effect plot
pdf(file=paste("./Plots/AbsHeadingChange_pooled_2sWindows_EffectBarplot_allodors.pdf",sep=""),
    width=8,height=7,bg="white")


par(lwd=5)


m.fit<-matrix(ncol=3,nrow=2)
colnames(m.fit)<-c("Col","But","Lin")
rownames(m.fit)<-c("S","PstS")
m.fit[,1]<-c(ColS.quants[1],ColPstS.quants[1])
m.fit[,2]<-c(ButS.quants[1],ButPstS.quants[1])
m.fit[,3]<-c(LinS.quants[1],LinPstS.quants[1])

m.lwr<-matrix(ncol=3,nrow=2)
colnames(m.lwr)<-c("Col","But","Lin")
rownames(m.lwr)<-c("S","PstS")
m.lwr[,1]<-c(ColS.quants[2],ColPstS.quants[2])
m.lwr[,2]<-c(ButS.quants[2],ButPstS.quants[2])
m.lwr[,3]<-c(LinS.quants[2],LinPstS.quants[2])

m.upr<-matrix(ncol=3,nrow=2)
colnames(m.upr)<-c("Col","But","Lin")
rownames(m.upr)<-c("S","PstS")
m.upr[,1]<-c(ColS.quants[3],ColPstS.quants[3])
m.upr[,2]<-c(ButS.quants[3],ButPstS.quants[3])
m.upr[,3]<-c(LinS.quants[3],LinPstS.quants[3])


ylim=c(-2,2)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="",main="Change in Abs heading")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.07,2),xpd=T,labels=rep(c("S","PstS"),3),cex=0.7)
text(x=barcenters[1,]+0.5,y=rep(par("usr")[3]-0.1,3),xpd=T,labels=rep(c("Col","But","Lin"),3),cex=1)



dev.off()















#--------------------------------------------------------------------------------------------------------------------------
# Presence density 
#--------------------------------------------------------------------------------------------------------------------------


#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)






########################
### Presence density ###
########################



# Parameters
#/////////////
Od="Linalool" # "Colony", "Butanol" or "Linalool"
view="back"  # view from the "back" (X & Z) or "top" (X & Y)
nbins=24   # number of bins in x and y 
#/////////////



# raw data 

dat<-dplyr::filter(AllData,
                   Protocol %in% c("ProtocolE","ProtocolB","ProtocolGb"),
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

dat$Ant.Z[dat$Ant.Z<0]<-0 #change everything below zero to zero, because 0 is the ground

# Normalization of coordinates for each trial & distribution metrics
datlist<-split(dat,f=dat$Test)  #separate the df into a list by trial id
datlist<-lapply(datlist,
                FUN=function(df){
                  
                  #centering coordinates around the head
                  df$Ant.X<-df$Ant.X-df$head.X
                  df$Ant.Y<-df$Ant.Y-df$head.Y
                  df$Ant.Z<-df$Ant.Z-df$head.Z
                  
                  # normalizing by maximum coordinates of all trials
                  df$Ant.X<-df$Ant.X/max(c(max(dat$Ant.X-dat$head.X,na.rm=T),abs(min(dat$Ant.X-dat$head.X,na.rm=T))))
                  df$Ant.Y<-df$Ant.Y/max(c(max(dat$Ant.Y-dat$head.Y,na.rm=T),abs(min(dat$Ant.Y-dat$head.Y,na.rm=T))))
                  df$Ant.Z<-df$Ant.Z/max(dat$Ant.Z,na.rm=T)
                  
                  return(df)
                  
                })
dat<-dplyr::bind_rows(datlist) #put df back together
dat<-dat[,-c(7:9)] #remove head coordinates






# Presence density matrices 

Presence.list<-list()
Hist.list<-list()
index=1
for(i in 1:2){
  
  phase<-c("PS","S")[i]
  subdat<-dplyr::filter(dat,TestPhase==phase)
  
  Hist.list[[index]]<-subdat #save for distributions below
  names(Hist.list)[index]<-phase
  
  x.bin<-seq(-1,1,length=nbins)  # bin borders between min and max coordinates (normalized), equal in size
  y.bin<-seq(-1,1,length=nbins)
  z.bin<-seq(min(dat$Ant.Z),1,length=nbins) # 0 at the head
  
  
  
  if(view=="top"){
    presence.df<-as.data.frame(table(findInterval(subdat$Ant.X,x.bin),findInterval(subdat$Ant.Y,y.bin)))
  } else{
    presence.df<-as.data.frame(table(findInterval(subdat$Ant.X,x.bin),findInterval(subdat$Ant.Z,z.bin)))
  }
  
  presence.df[,1]<-as.numeric(presence.df[,1])
  presence.df[,2]<-as.numeric(presence.df[,2])
  # Somehow, findInterval labels the intervals backwards? 
  # Need the biggest bin to always be max(nbins), because of how I normalise the coordinates.
  # it's only an issue for Y here
  if(view=="top"){
    presence.df[,2]<-presence.df[,2]+(nbins-max(presence.df[,2]))
  }
  
  spacemat<-matrix(data=0,nrow=nbins,ncol=nbins)
  spacemat[cbind(presence.df[,1],presence.df[,2])]<-presence.df[,3]
  
  #normalise by total observations
  spacemat<-spacemat/sum(spacemat)
  
  Presence.list[[index]]<-spacemat
  names(Presence.list)[index]<-phase
  index<-index+1
  
}










# plotting maps 


pdf(file=paste("./Plots/PresenceDensity_PooledS2s_",view,"view_",Od,"_ZeroHead.pdf",sep=""),
    width=7,height=4,bg="white",useDingbats=F)


cols<-viridis_pal()(20)
cols<-cols[-c(2:3)] #remove the second and third colors from the palette to have the low values stand out more
colfunc<-colorRampPalette(cols)



par(mfrow=c(1,2))


x=x.bin
if(view=="top"){
  y=y.bin
  ylab="y"
} else {
  y=z.bin
  ylab="z"
}



for(i in 1:length(Presence.list)){
  
  d<-Presence.list[[i]]

  image(x,y,d,col=colfunc(100),
        breaks=(seq(min(d),max(d),length.out=101)),
        xaxs="i",yaxs="i",xlab="x",ylab=ylab,
        sub=names(Presence.list)[i],cex.sub=1.5,font.sub=2)
  
}



dev.off()









# COLOR BAR
pdf(file=paste("./Plots/PresenceDensity_PooledS2s_",view,"view_",Od,"_ZeroHead_COLORBAR.pdf",sep=""),
    width=3,height=4,bg="white",useDingbats=F)

par(mfrow=c(1,2))

cols<-viridis_pal()(20)
cols<-cols[-c(2:3)] #remove the second and third colors from the palette to have the low values stand out more
colfunc<-colorRampPalette(cols)


for(i in 1:length(Presence.list)){
  
  d<-Presence.list[[i]]
  zmat<-matrix(ncol=100,nrow=1,data=seq(min(d),max(d),length.out=100))
  
  image(x=1,y=seq(min(d),max(d),length.out=100),zmat,col=colfunc(100),
        ylab="relative presence",xlab="",xaxt="n")
  
  if(i==1){
    title(main="PS")
  } else if(i==2){
    title(main="S")
  }
  
}

dev.off()







# plotting distributions 

if(view=="back"){
  
  pdf(file=paste("./Plots/ZcoordDensity_PooledS2s_",Od,".pdf",sep=""),width=5,height=5,bg="white",useDingbats=F)
  
  
  dps<-Hist.list$PS
  ds<-Hist.list$S
  
  adj=2
  
  
  dens<-density(dps$Ant.Z,adjust=adj)
  plot(dens$y,dens$x,type="l",lwd=5,
       col=FränziPurple,
       xlab="Density",
       ylab="Z",
       yaxs="i",xaxs="i",
       xlim=c(0,4))
  dens<-density(ds$Ant.Z,adjust=adj)
  lines(dens$y,dens$x,lwd=5,col=FränziOrange)
  

  dev.off()
  
  
  
  
} else if (view=="top"){
  
  
  pdf(file=paste("./Plots/XcoordDensity_PooledS2s_",Od,".pdf",sep=""),width=5,height=5,bg="white",useDingbats=F)
  
  
  dps<-Hist.list$PS
  ds<-Hist.list$S
  
  adj=2
  
  
  dens<-density(dps$Ant.X,adjust=adj)
  plot(dens$x,dens$y,type="l",lwd=5,
       col=FränziPurple,
       xlab="X",
       ylab="Density",
       yaxs="i",xaxs="i",
       ylim=c(0,1))
  dens<-density(ds$Ant.X,adjust=adj)
  lines(dens$x,dens$y,lwd=5,col=FränziOrange)
  
  
  dev.off()
  
  
}













######################################
### Modelling Distribution changes ###
######################################
# No need to run the code above for both a top view and a back view separately, 
# the code below uses the normalised coordinates data that includes both X and Z coordinates 
# and all is needed is the Hist.list


# X coordinates
### counting antenna presence in the center vs edges per trial ###

model.data<-dplyr::bind_rows(Hist.list)
datlist<-split(model.data,f=model.data$Test)  
datlist<-lapply(datlist,
                FUN=function(df){
                  
                  ps.ncenter<-sum(df$Ant.X[df$Count<0 & df$Count>=-2]>=(-0.5) & df$Ant.X[df$Count<0 & df$Count>=-2]<=0.5)
                  ps.nedge<-sum(df$Ant.X[df$Count<0 & df$Count>=-2]<(-0.5) | df$Ant.X[df$Count<0 & df$Count>=-2]>0.5)
                  
                  s.ncenter<-sum(df$Ant.X[df$Count<2 & df$Count>=0]>=(-0.5) & df$Ant.X[df$Count<2 & df$Count>=0]<=0.5)
                  s.nedge<-sum(df$Ant.X[df$Count<2 & df$Count>=0]<(-0.5) | df$Ant.X[df$Count<2 & df$Count>=0]>0.5)
                  
                  return(data.frame("Test"=unique(df$Test),
                                    "TestPhase"=rep(c("PS","S"),2),
                                    "location"=rep(c("center","edge"),each=2),
                                    "n"=c(ps.ncenter,s.ncenter,ps.nedge,s.nedge)))
                  
                })
model.data<-dplyr::bind_rows(datlist)

hist(model.data$n)



set.seed(171022)
brmsmod<-brm(n ~ location + TestPhase + location:TestPhase,
             data=model.data,
             family=gaussian(), #poisson isn't appropriate because observations are not independent
             iter=10000, chains=4, cores=4,
             backend="cmdstanr", control=list(adapt_delta=0.99,max_treedepth=20))



# pp check
pp<-pp_check(brmsmod,ndraws=200)
pp

#convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce 



# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))


# estimates and probabilities
PSedge.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),3)
PSedge.quants

Sedge.quants<-round(quantile((simval[,1]+simval[,2]+simval[,3]+simval[,4]),probs=c(0.5,0.025,0.975)),3)
Sedge.quants
print(paste("P(edge S>PS) = ",round(mean((simval[,1]+simval[,2]+simval[,3]+simval[,4])>(simval[,1]+simval[,2]),2)),sep=""))

PScenter.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),3)
PScenter.quants

Scenter.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),3)
Scenter.quants
print(paste("P(center S>PS) = ",round(mean((simval[,1]+simval[,3])>simval[,1],2)),sep=""))








# Z coordinates
### counting antenna presence above and below zero per trial ###

model.data<-dplyr::bind_rows(Hist.list)
datlist<-split(model.data,f=model.data$Test)  
datlist<-lapply(datlist,
                FUN=function(df){
                  
                  ps.nabove<-sum(df$Ant.Z[df$Count<0 & df$Count>=-2]>0)
                  ps.nbelow<-sum(df$Ant.Z[df$Count<0 & df$Count>=-2]<=0)
                  
                  s.nabove<-sum(df$Ant.Z[df$Count<2 & df$Count>=0]>0)
                  s.nbelow<-sum(df$Ant.Z[df$Count<2 & df$Count>=0]<=0)
                  
                  return(data.frame("Test"=unique(df$Test),
                                    "TestPhase"=rep(c("PS","S"),2),
                                    "location"=rep(c("above","below"),each=2),
                                    "n"=c(ps.nabove,s.nabove,ps.nbelow,s.nbelow)))
                  
                })
model.data<-dplyr::bind_rows(datlist)

hist(model.data$n)



set.seed(171022)
brmsmod<-brm(n ~ location + TestPhase + location:TestPhase,
             data=model.data,
             family=gaussian(), 
             iter=10000, chains=4, cores=4,
             backend="cmdstanr", control=list(adapt_delta=0.99,max_treedepth=20))



# pp check
pp<-pp_check(brmsmod,ndraws=200)
pp

#convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce 



# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))


# probs
PSbelow.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),3)
PSbelow.quants

Sbelow.quants<-round(quantile((simval[,1]+simval[,2]+simval[,3]+simval[,4]),probs=c(0.5,0.025,0.975)),3)
Sbelow.quants
print(paste("P(below S>PS) = ",round(mean((simval[,1]+simval[,2]+simval[,3]+simval[,4])>(simval[,1]+simval[,2]),2)),sep=""))

PSabove.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),3)
PSabove.quants

Sabove.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),3)
Sabove.quants
print(paste("P(above S>PS) = ",round(mean((simval[,1]+simval[,3])>simval[,1],2)),sep=""))


















#--------------------------------------------------------------------------------------------------------------------------
# Preference Index
#--------------------------------------------------------------------------------------------------------------------------


# Videos were tracked with BlobMaster3000 (https://github.com/YannickGuenzel/BlobMaster3000)
# and coordinates were pre-processed with code used in Günzel & McCollum et al. 2021. (https://github.com/Couzin-Fuchs-Lab/Guenzel-McCollum-et-al)
# From that, further analysis is done with the following code.




#load data set
ArenaData<-read.csv("../Data/ArenaData_ActiveSmelling.csv",sep=";",dec=".",header=T,check.names=F)



##########
### PI ###
##########

# Get relative time spent in location 

sheltersize=mean(ArenaData$shelter_r[ArenaData$Odor=="Lin"]) #shelter radius of Lin trials are 6.75cm
scaled1cm=sheltersize/6.75 #1cm in the scaled system

datlist<-split(ArenaData,f=ArenaData$Test)
datlist<-lapply(datlist,
                function(df){
                  
                  df$inOdor<-df$dist2odor<=df$shelter_r+scaled1cm*4 # 4cm buffer around shelter, T or F
                  df$inCtrl<-df$dist2control<=df$shelter_r+scaled1cm*4 
                  df$outsideOdor<-df$dist2odor>df$shelter_r+scaled1cm*4 & 
                    df$dist2center<=(1-(scaled1cm*3))
                  # animal outside shelter zone but not in edge zone (3cm) 
                  # (keeps instances where the shelter zone overlaps the edge zone)
                  df$outsideShelters<-df$dist2odor>df$shelter_r+scaled1cm*4 & 
                    df$dist2control>df$shelter_r+scaled1cm*4 &
                    df$dist2center<=(1-(scaled1cm*3))
                  
                  #proportion of time in locations
                  timeinOdor<-mean(df$inOdor)
                  timeinCtrl<-mean(df$inCtrl)
                  timeOutsideShelters<-mean(df$outsideShelters)
                  timeOutsideOdor<-mean(df$outsideOdor)
                  
                  #preference index = (TimeinOdor-TimeinCtrl) / (sum of both)
                  PI<-(sum(df$inOdor)-sum(df$inCtrl))/(sum(df$inOdor)+sum(df$inCtrl))
                  
                  
                  newdf<-data.frame("Test"=rep(unique(df$Test),4),
                                    "Odor"=rep(unique(df$Odor),4),
                                    "PropTime"=c(timeinOdor,timeinCtrl,timeOutsideShelters,timeOutsideOdor),
                                    "Location"=c("OdorShelter","CtrlShelter","NoShelter","OutsideOdor"),
                                    "PI"=rep(PI,4))
                  
                  
                  return(newdf)
                  
                })
PropTimeArena<-dplyr::bind_rows(datlist)







#############
### Model ###
#############


Model.dat<-PropTimeArena[,c("Odor","PI","Test")]
Model.dat<-Model.dat[!duplicated(Model.dat),]


set.seed(070922)
brmsmod<-brm(PI ~ Odor,data=Model.dat, 
             family=gaussian(),
             iter=10000, chains=4, cores=4, 
             backend="cmdstanr", control=list(adapt_delta=0.99,max_treedepth=20))



# pp check
pp<-pp_check(brmsmod,ndraws=200)
pp


#convergence check
plot(brmsmod)



#quick effect plot
ce<-conditional_effects(brmsmod)
ce



#draws
simval<-as.data.frame(as_draws_df(brmsmod))

#estimates
But<-round(quantile(simval[,1],probs=c(0.5,0.025,0.975)),2) 
Fec<-round(quantile(simval[,1]+simval[,2],probs=c(0.5,0.025,0.975)),2) 
Lin<-round(quantile(simval[,1]+simval[,3],probs=c(0.5,0.025,0.975)),2)


#probabilities
print(paste("P(But>0) = ",round(mean(simval[,1]>0),2),sep="")) 
print(paste("P(Fec>0) = ",round(mean((simval[,1]+simval[,2])>0),2),sep="")) 
print(paste("P(Lin>0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep="")) 



# effect plot

pdf(file=paste("./Plots/PI_effectplot.pdf",sep=""),
    width=5,height=9,bg="white")

par(lwd=5)

m.fit<-matrix(c(Fec[1],But[1],Lin[1]) )
rownames(m.fit)<-c("Fec","But","Lin")

m.lwr<-c(Fec[2],But[2],Lin[2]) 
m.upr<-c(Fec[3],But[3],Lin[3]) 


barcenters<-barplot(height=m.fit,plot=T,beside=T,
                    col=rgb(0,0,0,100/255),
                    ylim=c(-1,1),cex.axis=1,ylab="",main="PI")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr)
points(x=barcenters,y=m.fit,pch=19,cex=2)
axis(1,at=barcenters,
     labels=c("Fec","But","Lin"),lty=0)


dev.off()


















