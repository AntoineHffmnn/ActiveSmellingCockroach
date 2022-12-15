###########################################################################################################################
#                               Active smelling in the American cockroach: Fig. 4 & Fig. S2                               #
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



# Moving mean, median, sum, min, max or sd
# returns a df with center index of the bins & value
# if new.length != NULL, directly create a vector that can be added to a df instead of a separate output df
MovingFUN<-function(x,winsize,shift,FUN,new.length=NULL){
  
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

#for the overall range of the antennae
MovingRange<-function(x1,x2,winsize,shift,new.length=NULL){
  
  if(length(x1)!=length(x2)){
    stop("x1 and x2 have different lengths")
  }
  
  
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
  
  if(anyNA(x1)==TRUE){
    x1<-repeat.before(x1)
  }
  if(anyNA(x2)==TRUE){
    x2<-repeat.before(x2)
  }
  
  
  fun.out<-c()
  index<-c()
  i=1
  j=1
  while(i<(length(x1)-(winsize/2))){
    
    fun.out[j]<-abs(max(c(x1[i:(i+(winsize-1))],x2[i:(i+(winsize-1))]))-min(c(x1[i:(i+(winsize-1))],x2[i:(i+(winsize-1))])))  
    
    if(is.nan(fun.out[j])==T){ fun.out[j]<-NA }
    
    index[j]<-(i+(winsize/2))  
    i=i+shift
    j=j+1
    
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



# Euclidean distance between two vectors
euc.dist.vectors <- function(x1, x2){
  sqrt(sum((x1-x2)^2,na.rm=T))
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

data<-dplyr::filter(AllData,
                    Dataset=="3D")
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
# Stepping changes
#--------------------------------------------------------------------------------------------------------------------------



# These trials had no walking activity at all and were excluded
non.walking<-c("T520","T522","T527","T530","T533","T564","T575","T617","T645","T662","T675","T676","T677","T678",
               "T679","T680","T683","T684","T686","T691","T700","T703","T705","T778","T780","T787","T821","T825",
               "T829","T842","T843","T866","T895","T899","T900","T903","T907","T926","T945","T946","T950","T951")




######################################
#   Stepping rate avg times-series   #
######################################


stepping<-dplyr::filter(data,
                        Count>=-2,
                        Protocol == "ProtocolE",
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
pdf(file=paste("./Plots/AvgSteppingRate_MovingStream_timeseries.pdf",sep=""),width=6,height=15,bg="white")


par(mfrow=c(3,1))
for(i in unique(stepping$Odor)){
  
  if(grepl("Blank",i)==TRUE) { next }
  
  dat<-dplyr::filter(stepping,Odor==i)
  
  plot(dat$Count,dat$stepping,
       type="l",lwd=2,col="black",xaxs="i",yaxs="i",
       ylim=c(-0.2,0.35),
       xlim=c(-2,12),
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
                        Protocol == "ProtocolE",
                        Count>=-2,
)


stepping<-split(stepping,f=stepping$Test)
stepping<-lapply(stepping,
                 FUN=function(df){ 
                   
                   df$stepping<-df$LegsFreqavg-mean(df$LegsFreqavg[df$Count<0]) #center
                   
                   stepping.ps<-mean(df$stepping[df$Count<0],na.rm=T)
                   stepping.s1<-mean(df$stepping[df$Count>=0 & df$Count<2],na.rm=T) #static1
                   stepping.s2<-mean(df$stepping[df$Count>=2 & df$Count<5],na.rm=T) #shift1
                   stepping.s3<-mean(df$stepping[df$Count>=5 & df$Count<7],na.rm=T) #static2
                   stepping.s4<-mean(df$stepping[df$Count>=7 & df$Count<10],na.rm=T) #shift2
                   stepping.s5<-mean(df$stepping[df$Count>=10 & df$Count<12],na.rm=T) #static3
                   stepping.psts<-mean(df$stepping[df$Count>=12],na.rm=T)

                   df<-data.frame("Odor"=rep(unique(df$Odor),7),
                                  "Animal"=rep(unique(df$Animal),7),                                
                                  "Test"=rep(unique(df$Test),7),
                                  "TestPhase"=c("PS","S1","S2","S3","S4","S5","PstS"),
                                  "stepping"=c(stepping.ps,stepping.s1,stepping.s2,stepping.s3,
                                               stepping.s4,stepping.s5,stepping.psts)
                   )
                   
                   
                   return(df)
                 })
stepping<-dplyr::bind_rows(stepping)








### Modelling ###


modeldat<-dplyr::filter(stepping,
                        Odor !="Blank",
                        TestPhase %in% c("S1","S2","S3","S4","S5"))


#view overall distribution
hist(modeldat$stepping)



set.seed(170822)
brmsmod<-brm(stepping ~ Odor + TestPhase + Odor:TestPhase, 
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



# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))






# Estimates, 95% CrI, and probabilities
ColS1.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColS1.quants
print(paste("P(ColS1>0) = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))

ColS2.quants<-round(quantile((simval[,1]+simval[,2]+simval[,4]+simval[,8]),probs=c(0.5,0.025,0.975)),2)
ColS2.quants
print(paste("P(ColS2>0) = ",round(mean((simval[,1]+simval[,2]+simval[,4]+simval[,8])>0),2),sep=""))

ColS3.quants<-round(quantile((simval[,1]+simval[,2]+simval[,5]+simval[,10]),probs=c(0.5,0.025,0.975)),2)
ColS3.quants
print(paste("P(ColS3>0) = ",round(mean((simval[,1]+simval[,2]+simval[,5]+simval[,10])>0),2),sep=""))

ColS4.quants<-round(quantile((simval[,1]+simval[,2]+simval[,6]+simval[,12]),probs=c(0.5,0.025,0.975)),2)
ColS4.quants
print(paste("P(ColS4>0) = ",round(mean((simval[,1]+simval[,2]+simval[,6]+simval[,12])>0),2),sep=""))

ColS5.quants<-round(quantile((simval[,1]+simval[,2]+simval[,7]+simval[,14]),probs=c(0.5,0.025,0.975)),2)
ColS5.quants
print(paste("P(ColS5>0) = ",round(mean((simval[,1]+simval[,2]+simval[,7]+simval[,14])>0),2),sep=""))



ButS1.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ButS1.quants
print(paste("P(ButS1>0) = ",round(mean((simval[,1])>0),2),sep=""))

ButS2.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
ButS2.quants
print(paste("P(ButS2>0) = ",round(mean((simval[,1]+simval[,4])>0),2),sep=""))

ButS3.quants<-round(quantile((simval[,1]+simval[,5]),probs=c(0.5,0.025,0.975)),2)
ButS3.quants
print(paste("P(ButS3>0) = ",round(mean((simval[,1]+simval[,5])>0),2),sep=""))

ButS4.quants<-round(quantile((simval[,1]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
ButS4.quants
print(paste("P(ButS4>0) = ",round(mean((simval[,1]+simval[,6])>0),2),sep=""))

ButS5.quants<-round(quantile((simval[,1]+simval[,7]),probs=c(0.5,0.025,0.975)),2)
ButS5.quants
print(paste("P(ButS5>0) = ",round(mean((simval[,1]+simval[,7])>0),2),sep=""))



LinS1.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
LinS1.quants
print(paste("P(LinS1>0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))

LinS2.quants<-round(quantile((simval[,1]+simval[,3]+simval[,4]+simval[,9]),probs=c(0.5,0.025,0.975)),2)
LinS2.quants
print(paste("P(LinS2>0) = ",round(mean((simval[,1]+simval[,3]+simval[,4]+simval[,9])>0),2),sep=""))

LinS3.quants<-round(quantile((simval[,1]+simval[,3]+simval[,5]+simval[,11]),probs=c(0.5,0.025,0.975)),2)
LinS3.quants
print(paste("P(LinS3>0) = ",round(mean((simval[,1]+simval[,3]+simval[,5]+simval[,11])>0),2),sep=""))

LinS4.quants<-round(quantile((simval[,1]+simval[,3]+simval[,6]+simval[,13]),probs=c(0.5,0.025,0.975)),2)
LinS4.quants
print(paste("P(LinS4>0) = ",round(mean((simval[,1]+simval[,3]+simval[,6]+simval[,13])>0),2),sep=""))

LinS5.quants<-round(quantile((simval[,1]+simval[,3]+simval[,7]+simval[,15]),probs=c(0.5,0.025,0.975)),2)
LinS5.quants
print(paste("P(LinS5>0) = ",round(mean((simval[,1]+simval[,3]+simval[,7]+simval[,15])>0),2),sep=""))








#effect plot
pdf(file=paste("./Plots/SteppingRateChange_MovingStream_2sWindows_EffectBarplot_allodors.pdf",sep=""),
    width=8,height=7,bg="white")



par(lwd=5)

m.fit<-matrix(ncol=3,nrow=5)
colnames(m.fit)<-c("Col","But","Lin")
rownames(m.fit)<-c("S1","S2","S3","S4","S5")
m.fit[,1]<-c(ColS1.quants[1],ColS2.quants[1],ColS3.quants[1],ColS4.quants[1],ColS5.quants[1])
m.fit[,2]<-c(ButS1.quants[1],ButS2.quants[1],ButS3.quants[1],ButS4.quants[1],ButS5.quants[1])
m.fit[,3]<-c(LinS1.quants[1],LinS2.quants[1],LinS3.quants[1],LinS4.quants[1],LinS5.quants[1])

m.lwr<-matrix(ncol=3,nrow=5)
colnames(m.lwr)<-c("Col","But","Lin")
rownames(m.lwr)<-c("S1","S2","S3","S4","S5")
m.lwr[,1]<-c(ColS1.quants[2],ColS2.quants[2],ColS3.quants[2],ColS4.quants[2],ColS5.quants[2])
m.lwr[,2]<-c(ButS1.quants[2],ButS2.quants[2],ButS3.quants[2],ButS4.quants[2],ButS5.quants[2])
m.lwr[,3]<-c(LinS1.quants[2],LinS2.quants[2],LinS3.quants[2],LinS4.quants[2],LinS5.quants[2])

m.upr<-matrix(ncol=3,nrow=5)
colnames(m.upr)<-c("Col","But","Lin")
rownames(m.upr)<-c("S1","S2","S3","S4","S5")
m.upr[,1]<-c(ColS1.quants[3],ColS2.quants[3],ColS3.quants[3],ColS4.quants[3],ColS5.quants[3])
m.upr[,2]<-c(ButS1.quants[3],ButS2.quants[3],ButS3.quants[3],ButS4.quants[3],ButS5.quants[3])
m.upr[,3]<-c(LinS1.quants[3],LinS2.quants[3],LinS3.quants[3],LinS4.quants[3],LinS5.quants[3])


ylim=c(-0.25,0.35)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="[Hz]",main="Change in stepping rate - windows")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.01,3),xpd=T,labels=rep(c("S1","S2","S3","S4","S5"),3),cex=0.7)
text(x=barcenters[3,],y=rep(par("usr")[3]-0.025,3),xpd=T,labels=rep(c("Col","But","Lin"),3),cex=1)



dev.off()









#--------------------------------------------------------------------------------------------------------------------------
# Heading changes
#--------------------------------------------------------------------------------------------------------------------------




##########################################
#   Heading direction avg times-series   #
##########################################

# pooled first 2 seconds

heading<-dplyr::filter(data,
                       Count>-2,
                       Protocol == "ProtocolE") 


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
pdf(file=paste("./Plots/HeadingChange_MovingStream_timeseries.pdf",sep=""),
    width=6,height=15,bg="white")


par(mfrow=c(3,1))

for(i in unique(heading$Odor)){
  
  if(grepl("Blank",i)==TRUE) { next }
  
  dat<-dplyr::filter(heading,Odor==i)
  
  plot(dat$Count,dat$heading,
       type="l",lwd=2,col="black",xaxs="i",yaxs="i",
       ylim=c(-7,9),
       xlim=c(-2,14),
       main=paste(i,"Heading change"),
       ylab="[°]")
  
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


heading<-dplyr::filter(data, 
                        Protocol == "ProtocolE",
                        Count>=-2,
)


heading<-split(heading,f=heading$Test)
heading<-lapply(heading,
                 FUN=function(df){ 
                   
                   df$heading<-df$heading-mean(df$heading[df$Count<0]) #center
                   
                   heading.ps<-mean(df$heading[df$Count<0],na.rm=T)
                   heading.s1<-mean(df$heading[df$Count>=0 & df$Count<2],na.rm=T) #static1
                   heading.s2<-mean(df$heading[df$Count>=2 & df$Count<5],na.rm=T) #shift1
                   heading.s3<-mean(df$heading[df$Count>=5 & df$Count<7],na.rm=T) #static2
                   heading.s4<-mean(df$heading[df$Count>=7 & df$Count<10],na.rm=T) #shift2
                   heading.s5<-mean(df$heading[df$Count>=10 & df$Count<12],na.rm=T) #static3
                   heading.psts<-mean(df$heading[df$Count>=12],na.rm=T)
                   
                   df<-data.frame("Odor"=rep(unique(df$Odor),7),
                                  "Animal"=rep(unique(df$Animal),7),                                
                                  "Test"=rep(unique(df$Test),7),
                                  "TestPhase"=c("PS","S1","S2","S3","S4","S5","PstS"),
                                  "heading"=c(heading.ps,heading.s1,heading.s2,heading.s3,
                                              heading.s4,heading.s5,heading.psts)
                   )
                   
                   
                   return(df)
                 })
heading<-dplyr::bind_rows(heading)








### Modelling ###


modeldat<-dplyr::filter(heading,
                        Odor !="Blank",
                        TestPhase %in% c("S1","S2","S3","S4","S5"))


#view overall distribution
hist(modeldat$heading)



set.seed(170822)
brmsmod<-brm(heading ~ Odor + TestPhase + Odor:TestPhase, 
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



# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))




# Estimates, 95% CrI, and probabilities
ColS1.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColS1.quants
print(paste("P(ColS1>0) = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))

ColS2.quants<-round(quantile((simval[,1]+simval[,2]+simval[,4]+simval[,8]),probs=c(0.5,0.025,0.975)),2)
ColS2.quants
print(paste("P(ColS2>0) = ",round(mean((simval[,1]+simval[,2]+simval[,4]+simval[,8])>0),2),sep=""))

ColS3.quants<-round(quantile((simval[,1]+simval[,2]+simval[,5]+simval[,10]),probs=c(0.5,0.025,0.975)),2)
ColS3.quants
print(paste("P(ColS3>0) = ",round(mean((simval[,1]+simval[,2]+simval[,5]+simval[,10])>0),2),sep=""))

ColS4.quants<-round(quantile((simval[,1]+simval[,2]+simval[,6]+simval[,12]),probs=c(0.5,0.025,0.975)),2)
ColS4.quants
print(paste("P(ColS4>0) = ",round(mean((simval[,1]+simval[,2]+simval[,6]+simval[,12])>0),2),sep=""))

ColS5.quants<-round(quantile((simval[,1]+simval[,2]+simval[,7]+simval[,14]),probs=c(0.5,0.025,0.975)),2)
ColS5.quants
print(paste("P(ColS5>0) = ",round(mean((simval[,1]+simval[,2]+simval[,7]+simval[,14])>0),2),sep=""))



ButS1.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ButS1.quants
print(paste("P(ButS1>0) = ",round(mean((simval[,1])>0),2),sep=""))

ButS2.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
ButS2.quants
print(paste("P(ButS2>0) = ",round(mean((simval[,1]+simval[,4])>0),2),sep=""))

ButS3.quants<-round(quantile((simval[,1]+simval[,5]),probs=c(0.5,0.025,0.975)),2)
ButS3.quants
print(paste("P(ButS3>0) = ",round(mean((simval[,1]+simval[,5])>0),2),sep=""))

ButS4.quants<-round(quantile((simval[,1]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
ButS4.quants
print(paste("P(ButS4>0) = ",round(mean((simval[,1]+simval[,6])>0),2),sep=""))

ButS5.quants<-round(quantile((simval[,1]+simval[,7]),probs=c(0.5,0.025,0.975)),2)
ButS5.quants
print(paste("P(ButS5>0) = ",round(mean((simval[,1]+simval[,7])>0),2),sep=""))



LinS1.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
LinS1.quants
print(paste("P(LinS1>0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))

LinS2.quants<-round(quantile((simval[,1]+simval[,3]+simval[,4]+simval[,9]),probs=c(0.5,0.025,0.975)),2)
LinS2.quants
print(paste("P(LinS2>0) = ",round(mean((simval[,1]+simval[,3]+simval[,4]+simval[,9])>0),2),sep=""))

LinS3.quants<-round(quantile((simval[,1]+simval[,3]+simval[,5]+simval[,11]),probs=c(0.5,0.025,0.975)),2)
LinS3.quants
print(paste("P(LinS3>0) = ",round(mean((simval[,1]+simval[,3]+simval[,5]+simval[,11])>0),2),sep=""))

LinS4.quants<-round(quantile((simval[,1]+simval[,3]+simval[,6]+simval[,13]),probs=c(0.5,0.025,0.975)),2)
LinS4.quants
print(paste("P(LinS4>0) = ",round(mean((simval[,1]+simval[,3]+simval[,6]+simval[,13])>0),2),sep=""))

LinS5.quants<-round(quantile((simval[,1]+simval[,3]+simval[,7]+simval[,15]),probs=c(0.5,0.025,0.975)),2)
LinS5.quants
print(paste("P(LinS5>0) = ",round(mean((simval[,1]+simval[,3]+simval[,7]+simval[,15])>0),2),sep=""))








#effect plot
pdf(file=paste("./Plots/HeadingChange_MovingStream_2sWindows_EffectBarplot_allodors.pdf",sep=""),
    width=8,height=7,bg="white")



par(lwd=5)

m.fit<-matrix(ncol=3,nrow=5)
colnames(m.fit)<-c("Col","But","Lin")
rownames(m.fit)<-c("S1","S2","S3","S4","S5")
m.fit[,1]<-c(ColS1.quants[1],ColS2.quants[1],ColS3.quants[1],ColS4.quants[1],ColS5.quants[1])
m.fit[,2]<-c(ButS1.quants[1],ButS2.quants[1],ButS3.quants[1],ButS4.quants[1],ButS5.quants[1])
m.fit[,3]<-c(LinS1.quants[1],LinS2.quants[1],LinS3.quants[1],LinS4.quants[1],LinS5.quants[1])

m.lwr<-matrix(ncol=3,nrow=5)
colnames(m.lwr)<-c("Col","But","Lin")
rownames(m.lwr)<-c("S1","S2","S3","S4","S5")
m.lwr[,1]<-c(ColS1.quants[2],ColS2.quants[2],ColS3.quants[2],ColS4.quants[2],ColS5.quants[2])
m.lwr[,2]<-c(ButS1.quants[2],ButS2.quants[2],ButS3.quants[2],ButS4.quants[2],ButS5.quants[2])
m.lwr[,3]<-c(LinS1.quants[2],LinS2.quants[2],LinS3.quants[2],LinS4.quants[2],LinS5.quants[2])

m.upr<-matrix(ncol=3,nrow=5)
colnames(m.upr)<-c("Col","But","Lin")
rownames(m.upr)<-c("S1","S2","S3","S4","S5")
m.upr[,1]<-c(ColS1.quants[3],ColS2.quants[3],ColS3.quants[3],ColS4.quants[3],ColS5.quants[3])
m.upr[,2]<-c(ButS1.quants[3],ButS2.quants[3],ButS3.quants[3],ButS4.quants[3],ButS5.quants[3])
m.upr[,3]<-c(LinS1.quants[3],LinS2.quants[3],LinS3.quants[3],LinS4.quants[3],LinS5.quants[3])


ylim=c(-6,9)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="[Hz]",main="Change in heading - windows")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.01,3),xpd=T,labels=rep(c("S1","S2","S3","S4","S5"),3),cex=0.7)
text(x=barcenters[3,],y=rep(par("usr")[3]-1,3),xpd=T,labels=rep(c("Col","But","Lin"),3),cex=1)



dev.off()










#--------------------------------------------------------------------------------------------------------------------------
# Antenna - Stream Correlation
#--------------------------------------------------------------------------------------------------------------------------


#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)



# Run for each Odor


###############
# Correlation #
###############

# Parameters
#Stream center coordinates:
PEstream=c(rep(0,500+200),seq(0,30,length.out=300),rep(30,200),seq(30,60,length.out=300),rep(60,200+500)) # 500 PS, 1200 S, 500 PstS 
#/////////////////
Od="Colony"
win=100 #sliding window
sh=1 #shift size 
set.seed(210720) #for the jitter 
#////////////////


jit<-jitter(rep(0,length(PEstream)),factor=0.001) #adding a bit of noise to not have sd=0 and be able to calculate rho
stream<-PEstream


### CONTROL DATA
ctrldat<-dplyr::filter(AllData,
                       Protocol=="ProtocolE",
                       Odor=="Blank",
                       Count>-5)
ctrldat<-ctrldat[,c("Animal","Test","Count","TestPhase","Lant.X","Rant.X","head.X")]
ctrldat$stream<-stream


# Calculating Rho on each trial
datlist<-split(ctrldat,f=ctrldat$Test)  #separate the df into a list by trial id
correl.list<-lapply(datlist,
                    FUN=function(df){
                      
                      test<-unique(df$Test)
                      animal<-unique(df$Animal)
                      
                      # Normalizing the coordinates the same way as in Fig2
                      #centering coordinates around the head
                      df$Lant.X<-df$Lant.X-df$head.X
                      df$Rant.X<-df$Rant.X-df$head.X
                      
                      #normalizing stream coordinates 
                      df$stream<-df$stream/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                      df$stream<-df$stream+jit #add jitter
                      
                      # normalizing by maximum coordinates within trial
                      df$Lant.X.norm<-df$Lant.X/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                      df$Rant.X.norm<-df$Rant.X/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                      
                      
                      
                      
                      #moving mean - Right
                      df$Rant.X.Mean<-MovingFUN(df$Rant.X.norm,win,sh,new.length=nrow(df),FUN="mean")
                      df<-df[is.na(df$Rant.X.Mean)==F,] #cutting the parts with NAs for cor() to work (start and end)
                      #moving mean - Left
                      df$Lant.X.Mean<-MovingFUN(df$Lant.X.norm,win,sh,new.length=nrow(df),FUN="mean") 
                      df<-df[is.na(df$Lant.X.Mean)==F,]#cutting the parts with NAs for cor() to work (start and end)
                      
                      #PS
                      r<-cor(df[df$TestPhase=="PS","Rant.X.Mean"],df[df$TestPhase=="PS","stream"],method="pearson") #Right
                      l<-cor(df[df$TestPhase=="PS","Lant.X.Mean"],df[df$TestPhase=="PS","stream"],method="pearson") #Left
                      rho.PS<-data.frame("Animal"=animal,
                                         "Odor"="ctrl",
                                         "Test"=test,
                                         "TestPhase"="PS",
                                         "R"=r,
                                         "L"=l)
                      #S (all)
                      r<-cor(df[df$TestPhase=="S","Rant.X.Mean"],
                             df[df$TestPhase=="S","stream"],method="pearson")
                      l<-cor(df[df$TestPhase=="S","Lant.X.Mean"],
                             df[df$TestPhase=="S","stream"],method="pearson")
                      rho.S<-data.frame("Animal"=animal,
                                        "Odor"="ctrl",
                                        "Test"=test,
                                        "TestPhase"="S",
                                        "R"=r,
                                        "L"=l)
                      #PstS
                      r<-cor(df[df$TestPhase=="PstS","Rant.X.Mean"],df[df$TestPhase=="PstS","stream"],method="pearson")
                      l<-cor(df[df$TestPhase=="PstS","Lant.X.Mean"],df[df$TestPhase=="PstS","stream"],method="pearson")
                      rho.PstS<-data.frame("Animal"=animal,
                                           "Odor"="ctrl",
                                           "Test"=test,
                                           "TestPhase"="PstS",
                                           "R"=r,
                                           "L"=l)
                      
                      rho<-rbind(rho.PS,rho.S,rho.PstS)

                      
                      return(rho)
                      
                    })
correl.ctrl<-dplyr::bind_rows(correl.list) 



### ODOR DATA
dat<-dplyr::filter(AllData,
                   Protocol=="ProtocolE",
                   Odor==Od,
                   Count>-5)
dat<-dat[,c("Animal","Test","Count","TestPhase","Lant.X","Rant.X","head.X")]
dat$stream<-stream 

# Calculating Rho on each trial
datlist<-split(dat,f=dat$Test)  #separate the df into a list by trial id
correl.list<-lapply(datlist,
                    FUN=function(df){
                      
                      test<-unique(df$Test)
                      animal<-unique(df$Animal)
                      
                      # Normalizing the coordinates the same way as in Fig2
                      #centering coordinates around the head
                      df$Lant.X<-df$Lant.X-df$head.X
                      df$Rant.X<-df$Rant.X-df$head.X
                      
                      #normalizing stream coordinates
                      df$stream<-df$stream/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                      df$stream<-df$stream+jit #add jitter
                      
                      # normalizing by maximum coordinates within trial
                      df$Lant.X.norm<-df$Lant.X/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                      df$Rant.X.norm<-df$Rant.X/max(c(abs(df$Lant.X),abs(abs(df$Rant.X)),na.rm=T))
                      
                      
                      
                      
                      
                      #moving mean - Right
                      df$Rant.X.Mean<-MovingFUN(df$Rant.X.norm,win,sh,new.length=nrow(df),FUN="mean") 
                      df<-df[is.na(df$Rant.X.Mean)==F,] #cutting the parts with NAs for cor() to work (start and end)
                      #moving mean - Left
                      df$Lant.X.Mean<-MovingFUN(df$Lant.X.norm,win,sh,new.length=nrow(df),FUN="mean") 
                      df<-df[is.na(df$Lant.X.Mean)==F,]#cutting the parts with NAs for cor() to work (start and end)
                      
                      #PS
                      r<-cor(df[df$TestPhase=="PS","Rant.X.Mean"],df[df$TestPhase=="PS","stream"],method="pearson") #Right
                      l<-cor(df[df$TestPhase=="PS","Lant.X.Mean"],df[df$TestPhase=="PS","stream"],method="pearson") #Left
                      rho.PS<-data.frame("Animal"=animal,
                                         "Odor"=Od,
                                         "Test"=test,
                                         "TestPhase"="PS",
                                         "R"=r,
                                         "L"=l)
                      #S (all)
                      r<-cor(df[df$TestPhase=="S","Rant.X.Mean"],
                             df[df$TestPhase=="S","stream"],method="pearson")
                      l<-cor(df[df$TestPhase=="S","Lant.X.Mean"],
                             df[df$TestPhase=="S","stream"],method="pearson")
                      rho.S<-data.frame("Animal"=animal,
                                        "Odor"=Od,
                                        "Test"=test,
                                        "TestPhase"="S",
                                        "R"=r,
                                        "L"=l)
                      #PstS
                      r<-cor(df[df$TestPhase=="PstS","Rant.X.Mean"],df[df$TestPhase=="PstS","stream"],method="pearson")
                      l<-cor(df[df$TestPhase=="PstS","Lant.X.Mean"],df[df$TestPhase=="PstS","stream"],method="pearson")
                      rho.PstS<-data.frame("Animal"=animal,
                                           "Odor"=Od,
                                           "Test"=test,
                                           "TestPhase"="PstS",
                                           "R"=r,
                                           "L"=l)
                      
                      rho<-rbind(rho.PS,rho.S,rho.PstS)
                      
                      
                      return(rho)
                      
                    })
correl.odor<-dplyr::bind_rows(correl.list) 



correl.PE<-rbind(correl.ctrl,correl.odor) 








#################
### Modelling ###
#################


ModelData<-correl.PE

ModelData$Odor[ModelData$Odor!="ctrl"]<-"odor"

ModelData$group<-paste(ModelData$Odor,ModelData$TestPhase,sep="_")



# fit
set.seed(180321)
brmsmod<-brm(R ~ group + (1|Animal),
             # L ~ group + (1|Animal),
             data=ModelData,
             family=student(),
             iter=10000, chains=4, cores=4, backend="cmdstanr", 
             control=list(adapt_delta=0.99,max_treedepth=20))



# pp check
pp<-pp_check(brmsmod,ndraws=200) + xlim(-2,2)
pp

#convergence
plot(brmsmod)

#quick effect plot
ce<-conditional_effects(brmsmod)
ce


# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))


# Estimates, 95% CrI, and probabilities
odorS<-round(quantile((simval[,1]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
odorS
print(paste("P(odor S > 0) = ",round(mean((simval[,1]+simval[,6])>0),2),sep=""))

ctrlS<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
ctrlS
print(paste("P(ctrl S > 0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))

# odor > ctrl?
print(paste("P(odor > ctrl) = ",round(mean((simval[,1]+simval[,6])>(simval[,1]+simval[,3])),2),sep=""))





# effect plot
# values are specified manually from the results of the code runs

pdf(file=paste("./Plots/FIG3/pdfs/Pearson_MovingStream_EffectBarplot_allodors.pdf",sep=""),
    width=5,height=6,bg="white",useDingbats=F)


m.fit<-matrix(ncol=2,nrow=4)
colnames(m.fit)<-c("L","R")
rownames(m.fit)<-c("ctrlS","ColS","ButS","LinS")
m.fit[,1]<-c(-0.47,0.5,-0.07,-0.26) #L
m.fit[,2]<-c(0.21,0.41,0.06,0.23) #R

m.lwr<-matrix(ncol=2,nrow=4)
colnames(m.lwr)<-c("L","R")
rownames(m.lwr)<-c("ctrlS","ColS","ButS","LinS")
m.lwr[,1]<-c(-0.59,-0.37,-0.26,-0.32) 
m.lwr[,2]<-c(-0.28,0.35,-0.02,-0.07) 

m.upr<-matrix(ncol=2,nrow=4)
colnames(m.upr)<-c("L","R")
rownames(m.upr)<-c("ctrlS","ColS","ButS","LinS")
m.upr[,1]<-c(0.01,0.62,0.05,-0.19) 
m.upr[,2]<-c(0.34,0.47,0.39,0.52) 



ylim=c(-0.7,0.7)


par(lwd=5)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    col=c(rgb(0,0,0,100,maxColorValue=255),rep(rgb(firebrick3[1],firebrick3[2],firebrick3[3],100,maxColorValue=255),3),
                          rgb(0,0,0,100,maxColorValue=255),rep(rgb(dodgerblue3[1],dodgerblue3[2],dodgerblue3[3],100,maxColorValue=255),3)),
                    border=c("black",rep("firebrick3",3),
                             "black",rep("dodgerblue3",3)),
                    xaxt="n",ylim=ylim,
                    cex.axis=1,ylab="mean rho",main="stream-sweepcenter correlation")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col=c("black",rep("firebrick3",3),
                                                             "black",rep("dodgerblue3",3)))
text(x=barcenters,y=rep(-0.68,8),labels=rep(c("Ctrl","Col","But","Lin"),2),xpd=T,srt=45)
par(lwd=1)



dev.off()












#--------------------------------------------------------------------------------------------------------------------------
# Distance to the moving stream
#--------------------------------------------------------------------------------------------------------------------------



#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)




###########################
# Time-series and windows #
###########################


# Stream center coordinates:
PEstream=c(rep(0,400+200),seq(0,30,length.out=300),rep(30,200),seq(30,60,length.out=300),rep(60,200+500)) # 500 PS, 1200 S, 500 PstS 


data<-dplyr::filter(AllData,
                    Protocol == "ProtocolE",
                    Count>-4)

#stream position
stream<-PEstream


# Calculating the distance to the stream for each frame on each trial
### For the time-series
datlist<-split(data,f=data$Animal) 
datlist<-lapply(datlist,
                FUN=function(df){
                  
                  
                  ### CTRL ###
                  
                  dctrl<-dplyr::filter(df,Odor=="Blank")
                  
                  # Normalizing the coordinates the same way as in Fig2
                  #centering coordinates around the head
                  dctrl$Lant.X<-dctrl$Lant.X-dctrl$head.X
                  dctrl$Rant.X<-dctrl$Rant.X-dctrl$head.X
                  # normalizing stream 
                  stream.ctrl<-stream/max(c(abs(dctrl$Lant.X),abs(dctrl$Rant.X),na.rm=T))
                  # normalizing by maximum coordinates within trial
                  dctrl$Lant.X.norm<-dctrl$Lant.X/max(c(abs(dctrl$Lant.X),abs(dctrl$Rant.X),na.rm=T))
                  dctrl$Rant.X.norm<-dctrl$Rant.X/max(c(abs(dctrl$Lant.X),abs(dctrl$Rant.X),na.rm=T))
                  
                  #Right
                  dctrl$Rdist<-abs(stream.ctrl-dctrl$Rant.X.norm)
                  dctrl$Rdist<-dctrl$Rdist-mean(dctrl[dctrl$Count<=0,"Rdist"])
                  #Left
                  dctrl$Ldist<-abs(stream.ctrl-dctrl$Lant.X.norm)
                  dctrl$Ldist<-dctrl$Ldist-mean(dctrl[dctrl$Count<=0,"Ldist"])
                  
                  
                  
                  ### ODOR ###
                  
                  dod<-dplyr::filter(df,Odor == "Colony") #time series for colony only, can be changed
                  
                  # Normalizing the coordinates the same way as in Fig2
                  #centering coordinates around the head
                  dod$Lant.X<-dod$Lant.X-dod$head.X
                  dod$Rant.X<-dod$Rant.X-dod$head.X
                  # normalizing stream 
                  stream.odor<-stream/max(c(abs(dod$Lant.X),abs(dod$Rant.X),na.rm=T))
                  # normalizing by maximum coordinates within trial
                  dod$Lant.X.norm<-dod$Lant.X/max(c(abs(dod$Lant.X),abs(dod$Rant.X),na.rm=T))
                  dod$Rant.X.norm<-dod$Rant.X/max(c(abs(dod$Lant.X),abs(dod$Rant.X),na.rm=T))
                  
                  #Right
                  dod$Rdist<-abs(stream.odor-dod$Rant.X.norm)
                  dod$Rdist<-dod$Rdist-mean(dod[dod$Count<=0,"Rdist"])
                  #Left
                  dod$Ldist<-abs(stream.odor-dod$Lant.X.norm)
                  dod$Ldist<-dod$Ldist-mean(dod[dod$Count<=0,"Ldist"])
                  
                  
                  
                  #ctrl correction
                  dod$Rdist.correct<-dod$Rdist-dctrl$Rdist
                  dod$Ldist.correct<-dod$Ldist-dctrl$Ldist
                  
                  
                  return(dod)
                  
                })
data.ts<-dplyr::bind_rows(datlist) #put df back together



### For static P1, P2, P3 and each antenna
datlist<-split(data,f=data$Animal) 
datlist<-lapply(datlist,
                FUN=function(df){
                  
                  
                  ### CTRL ###
                  
                  d<-dplyr::filter(df,Odor=="Blank")
                  
                  # Normalizing the coordinates the same way as in Fig2
                  #centering coordinates around the head
                  d$Lant.X<-d$Lant.X-d$head.X
                  d$Rant.X<-d$Rant.X-d$head.X
                  # normalizing stream 
                  stream.ctrl<-stream/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                  # normalizing by maximum coordinates within trial
                  d$Lant.X.norm<-d$Lant.X/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                  d$Rant.X.norm<-d$Rant.X/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                  
                  #Right
                  d$Rdist<-abs(stream.ctrl-d$Rant.X.norm)
                  #Left
                  d$Ldist<-abs(stream.ctrl-d$Lant.X.norm)
                  
                  # mean per static
                  Rdist.PS1.ctrl<-mean(d[d$Count<=-1,"Rdist"])
                  Rdist.PS2.ctrl<-mean(d[d$Count>-1 & d$Count<=2,"Rdist"])
                  Rdist.PS3.ctrl<-mean(d[d$Count>-2 & d$Count<=0,"Rdist"])
                  Rdist.static1.ctrl<-mean(d[d$Count>0 & d$Count<=2,"Rdist"])
                  Rdist.static2.ctrl<-mean(d[d$Count>5 & d$Count<=7,"Rdist"])
                  Rdist.static3.ctrl<-mean(d[d$Count>10 & d$Count<=12,"Rdist"])
                  
                  Ldist.PS1.ctrl<-mean(d[d$Count<=-1,"Ldist"])
                  Ldist.PS2.ctrl<-mean(d[d$Count>-1 & d$Count<=2,"Ldist"])
                  Ldist.PS3.ctrl<-mean(d[d$Count>-2 & d$Count<=0,"Ldist"])
                  Ldist.static1.ctrl<-mean(d[d$Count>0 & d$Count<=2,"Ldist"])
                  Ldist.static2.ctrl<-mean(d[d$Count>5 & d$Count<=7,"Ldist"])
                  Ldist.static3.ctrl<-mean(d[d$Count>10 & d$Count<=12,"Ldist"])
                  
                  # delta
                  Rdist.baseline.ctrl.delta<-Rdist.PS2.ctrl-Rdist.PS1.ctrl #baseline variation before the shift
                  Rdist.static1.ctrl.delta<-Rdist.static1.ctrl-Rdist.PS3.ctrl #delta with only 2s PS
                  Rdist.static2.ctrl.delta<-Rdist.static2.ctrl-Rdist.static1.ctrl
                  Rdist.static3.ctrl.delta<-Rdist.static3.ctrl-Rdist.static1.ctrl
                  
                  Ldist.baseline.ctrl.delta<-Ldist.PS2.ctrl-Ldist.PS1.ctrl #baseline variation before the shift
                  Ldist.static1.ctrl.delta<-Ldist.static1.ctrl-Ldist.PS3.ctrl #delta with only 2s PS
                  Ldist.static2.ctrl.delta<-Ldist.static2.ctrl-Ldist.static1.ctrl
                  Ldist.static3.ctrl.delta<-Ldist.static3.ctrl-Ldist.static1.ctrl
                  
                  
                  
                  
                  
                  ### ODOR ###
                  
                  dflist<-split(df,f=df$Odor)
                  dflist<-lapply(dflist,
                                 FUN=function(d){
                                   
                                   
                                   # Normalizing the coordinates the same way as in Fig2
                                   #centering coordinates around the head
                                   d$Lant.X<-d$Lant.X-d$head.X
                                   d$Rant.X<-d$Rant.X-d$head.X
                                   # normalizing stream 
                                   stream.odor<-stream/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                                   # normalizing by maximum coordinates within trial
                                   d$Lant.X.norm<-d$Lant.X/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                                   d$Rant.X.norm<-d$Rant.X/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                                   
                                   #Right
                                   d$Rdist<-abs(stream.odor-d$Rant.X.norm)
                                   #Left
                                   d$Ldist<-abs(stream.odor-d$Lant.X.norm)
                                   
                                   # mean per static
                                   
                                   Rdist.PS1<-mean(d[d$Count<=-2,"Rdist"])
                                   Rdist.PS2<-mean(d[d$Count>-2 & d$Count<=0,"Rdist"])
                                   Rdist.static1<-mean(d[d$Count>0 & d$Count<=2,"Rdist"])
                                   Rdist.static2<-mean(d[d$Count>5 & d$Count<=7,"Rdist"])
                                   Rdist.static3<-mean(d[d$Count>10 & d$Count<=12,"Rdist"])
                                   
                                   Ldist.PS1<-mean(d[d$Count<=-2,"Ldist"])
                                   Ldist.PS2<-mean(d[d$Count>-2 & d$Count<=0,"Ldist"]) 
                                   Ldist.static1<-mean(d[d$Count>0 & d$Count<=2,"Ldist"])
                                   Ldist.static2<-mean(d[d$Count>5 & d$Count<=7,"Ldist"])
                                   Ldist.static3<-mean(d[d$Count>10 & d$Count<=12,"Ldist"])
                                   
                                   # delta
                                   Rdist.PS.delta<-Rdist.PS2-Rdist.PS1
                                   Rdist.static1.delta<-Rdist.static1-Rdist.PS2
                                   Rdist.static2.delta<-Rdist.static2-Rdist.PS2
                                   Rdist.static3.delta<-Rdist.static3-Rdist.PS2
                                   
                                   Ldist.PS.delta<-Ldist.PS2-Ldist.PS1
                                   Ldist.static1.delta<-Ldist.static1-Ldist.PS2
                                   Ldist.static2.delta<-Ldist.static2-Ldist.PS2
                                   Ldist.static3.delta<-Ldist.static3-Ldist.PS2
                                   
                                   
                                   # control correction
                                   Rdist.static1.delta.corrected<-Rdist.static1.delta-Rdist.baseline.ctrl.delta
                                   Rdist.static2.delta.corrected<-Rdist.static2.delta-Rdist.static2.ctrl.delta
                                   Rdist.static3.delta.corrected<-Rdist.static3.delta-Rdist.static3.ctrl.delta
                                   
                                   Ldist.static1.delta.corrected<-Ldist.static1.delta-Ldist.baseline.ctrl.delta
                                   Ldist.static2.delta.corrected<-Ldist.static2.delta-Ldist.static2.ctrl.delta
                                   Ldist.static3.delta.corrected<-Ldist.static3.delta-Ldist.static3.ctrl.delta
                                   
                                   
                                   # avg of antennae and all windows after correction
                                   dist.change.corrected.avg<-mean(c(Ldist.static1.delta.corrected,
                                                                     Ldist.static2.delta.corrected,
                                                                     Ldist.static3.delta.corrected,
                                                                     Rdist.static1.delta.corrected,
                                                                     Rdist.static2.delta.corrected,
                                                                     Rdist.static3.delta.corrected))
                                   
                                   
                                   
                                   return(data.frame("Animal"=rep(unique(d$Animal),6),
                                                     "Test"=rep(unique(d$Test),6),
                                                     "Odor"=rep(unique(d$Odor),6),
                                                     "Ant"=rep(c("L","R"),each=3),
                                                     "window"=rep(c("static1","static2","static3"),2),
                                                     "dist.change.corrected"=c(Ldist.static1.delta.corrected,
                                                                               Ldist.static2.delta.corrected,
                                                                               Ldist.static3.delta.corrected,
                                                                               Rdist.static1.delta.corrected,
                                                                               Rdist.static2.delta.corrected,
                                                                               Rdist.static3.delta.corrected),
                                                     "dist.change.corrected.avg"=rep(dist.change.corrected.avg,6)
                                   )
                                   )
                                   
                                   
                                   
                                 })
                  df<-dplyr::bind_rows(dflist) #put df back together
                  
                  
                  
                  return(df)
                  
                  
                  
                })
data<-dplyr::bind_rows(datlist) #put df back together

data.model<-data







### Plotting the time-series 


#avg
n<-length(unique(data.ts$Test))

Lsem<-aggregate(Ldist.correct~Count,data=data.ts,FUN="sd",na.rm=T)
Lsem$Ldist.correct<-Lsem$Ldist.correct/sqrt(n)
Lsem$smooth<-MovingFUN(Lsem$Ldist.correct,100,1,FUN="mean",new.length=nrow(Lsem))
Lsem$smooth<-repeat.before(Lsem$smooth)

Rsem<-aggregate(Rdist.correct~Count,data=data.ts,FUN="sd",na.rm=T)
Rsem$Rdist.correct<-Rsem$Rdist.correct/sqrt(n)
Rsem$smooth<-MovingFUN(Rsem$Rdist.correct,100,1,FUN="mean",new.length=nrow(Rsem))
Rsem$smooth<-repeat.before(Rsem$smooth)

data.ts.l<-aggregate(Ldist.correct~Count,data=data.ts,FUN="mean",na.rm=T)
data.ts.l$smooth<-MovingFUN(data.ts.l$Ldist.correct,100,1,FUN="mean",new.length=nrow(data.ts.l))
data.ts.l$smooth<-repeat.before(data.ts.l$smooth)

data.ts.r<-aggregate(Rdist.correct~Count,data=data.ts,FUN="mean",na.rm=T)
data.ts.r$smooth<-MovingFUN(data.ts.r$Rdist.correct,100,1,FUN="mean",new.length=nrow(data.ts.r))
data.ts.r$smooth<-repeat.before(data.ts.r$smooth)





pdf(file=paste("./Plots/FIG3/pdfs/StreamDist_withintrialctrlcorrected_MovingStream_",Od,".pdf",sep=""),
    width=10,height=8,bg="white",useDingbats=F)


ylim<-c(-0.4,0.2)
par(mar=c(5,5,4,2)+0.1)

plot(data.ts.r$Count,data.ts.r$smooth-200,type="l",col="dodgerblue3",
     main=paste(Od,"- Average corrected distance of sweep center to stream center"),
     ylab="Corrected distance [mm]",xlab="Time from odor onset [s]",
     ylim=ylim,
     xlim=c(-2,14),
     cex.main=1,font.main=2,cex.axis=1,font.axis=2,cex.lab=1,font.lab=2)
abline(h=0,lty=2)
#L
polygon(x=c(data.ts.l$Count,rev(data.ts.l$Count)),
        y=c(data.ts.l$smooth+Lsem$smooth,rev(data.ts.l$smooth-Lsem$smooth)),  #sem
        col=rgb(firebrick3[1],firebrick3[2],firebrick3[3],alpha=70,names=NULL,maxColorValue=255),border=NA)
lines(data.ts.l$Count,data.ts.l$smooth,col="firebrick3",lwd=4)
#R
polygon(x=c(data.ts.r$Count,rev(data.ts.r$Count)),
        y=c(data.ts.r$smooth+Rsem$smooth,rev(data.ts.r$smooth-Rsem$smooth)),  #sem
        col=rgb(dodgerblue3[1],dodgerblue3[2],dodgerblue3[3],alpha=70,names=NULL,maxColorValue=255),border=NA)
lines(data.ts.r$Count,data.ts.r$smooth,col="dodgerblue3",lwd=4)


abline(v=c(0,12))



dev.off()


#####









#############
# Modelling #
#############

# data
data.model.col<-dplyr::filter(data.model,Odor=="Colony")
hist(data.model.col$dist.change.corrected)
# for the avg overall of an odor
data.model.odoravg<-data.model
data.model.odoravg<-data.model.odoravg[!duplicated(data.model.odoravg), ] #remove duplicates
hist(data.model.odoravg$dist.change.corrected.avg)



# Fit
set.seed(090421)
brmsmod<-brm(# dist.change.corrected ~ Ant + window + Ant:window + (1|Animal), #for Colony
             #  data=data.model.col, #for Colony
             dist.change.corrected.avg ~ Odor + (1|Animal), #for the overall avg per odor
             data=data.model.odoravg, #for for the overall avg per odor
             family=student(),
             iter=10000, chains=4, cores=4, 
             backend="cmdstanr", control=list(adapt_delta=0.99,max_treedepth=20))




# pp check
pp<-pp_check(brmsmod,ndraws=200) + xlim(-2,2)
pp


#convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce




# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))



# Estimates, 95% CrI, and probabilities

# for Colony left & right antennae
ColP1L.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ColP1L.quants
print(paste("P(ColP1L>0) = ",round(mean((simval[,1])>0),2),sep=""))

ColP2L.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
ColP2L.quants
print(paste("P(ColP2L>0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))

ColP3L.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
ColP3L.quants
print(paste("P(ColP3L>0) = ",round(mean((simval[,1]+simval[,4])>0),2),sep=""))

ColP1R.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColP1R.quants
print(paste("P(ColP1R>0) = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))

ColP2R.quants<-round(quantile((simval[,1]+simval[,2]+simval[,3]+simval[,5]),probs=c(0.5,0.025,0.975)),2)
ColP2R.quants
print(paste("P(ColP2R>0) = ",round(mean((simval[,1]+simval[,2]+simval[,3]+simval[,5])>0),2),sep=""))

ColP3R.quants<-round(quantile((simval[,1]+simval[,2]+simval[,4]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
ColP3R.quants
print(paste("P(ColP3R>0) = ",round(mean((simval[,1]+simval[,2]+simval[,4]+simval[,6])>0),2),sep=""))



# for overall avg per odor
Col.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
Col.quants
print(paste("P(Col>0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))

But.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
But.quants
print(paste("P(But>0) = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))

Lin.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
Lin.quants
print(paste("P(Lin>0) = ",round(mean((simval[,1]+simval[,4])>0),2),sep=""))








# Plotting effect plots

# for colony
pdf(file=paste("./Plots/DistanceChanges_ctrlcorrectionwithinAnimal_MovingStream_Colony_EffectBarplot.pdf",sep="_"),
  width=5,height=7,bg="white",useDingbats=F)


par(lwd=5)

m.fit<-matrix(ncol=2,nrow=3)
colnames(m.fit)<-c("L","R")
rownames(m.fit)<-c("P1","P2","P3")
m.fit[,1]<-c(ColP1L.quants[1],ColP2L.quants[1],ColP3L.quants[1])
m.fit[,2]<-c(ColP1R.quants[1],ColP2R.quants[1],ColP3R.quants[1])

m.lwr<-matrix(ncol=2,nrow=3)
colnames(m.lwr)<-c("L","R")
rownames(m.lwr)<-c("P1","P2","P3")
m.lwr[,1]<-c(ColP1L.quants[2],ColP2L.quants[2],ColP3L.quants[2])
m.lwr[,2]<-c(ColP1R.quants[2],ColP2R.quants[2],ColP3R.quants[2])

m.upr<-matrix(ncol=2,nrow=3)
colnames(m.upr)<-c("L","R")
rownames(m.upr)<-c("P1","P2","P3")
m.upr[,1]<-c(ColP1L.quants[3],ColP2L.quants[3],ColP3L.quants[3])
m.upr[,2]<-c(ColP1R.quants[3],ColP2R.quants[3],ColP3R.quants[3])



ylim=c(-0.5,0.3)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="",main="Corrected change in Distance - Colony")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.01,2),xpd=T,labels=rep(c("P1","P2","P3"),2),cex=0.7)
text(x=barcenters[2,c(1,2)],y=rep(-0.6,2),xpd=T,labels=c("L","R"),cex=0.7)


dev.off()






# avg per odor
pdf(file=paste("./Plots/DistanceChanges_ctrlcorrectionwithinAnimal_MovingStream_AvgPerOdor_EffectBarplot.pdf",sep="_"),
    width=4,height=8,bg="white",useDingbats=F)


par(lwd=5)

m.fit<-matrix(ncol=1,nrow=3)
rownames(m.fit)<-c("Col","But","Lin")
m.fit[,1]<-c(Col.quants[1],But.quants[1],Lin.quants[1])

m.lwr<-matrix(ncol=1,nrow=3)
rownames(m.lwr)<-c("Col","But","Lin")
m.lwr[,1]<-c(Col.quants[2],But.quants[2],Lin.quants[2])

m.upr<-matrix(ncol=1,nrow=3)
rownames(m.upr)<-c("Col","But","Lin")
m.upr[,1]<-c(Col.quants[3],But.quants[3],Lin.quants[3])



ylim=c(-0.3,0.1)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="",main="Corrected change in Distance
- overall avg")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.01,2),xpd=T,labels=rep(c("Col","But","Lin"),2),cex=0.7)


dev.off()















#--------------------------------------------------------------------------------------------------------------------------
# Antennal Range
#--------------------------------------------------------------------------------------------------------------------------



#load data set
AllData<-read.csv("../Data/AllData_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)




###########################
# Time-series and windows #
###########################

data<-dplyr::filter(AllData,
                    Protocol == "ProtocolE",
                    Count>-4)


# Calculating the range for each frame on each trial
### For the time-series
datlist<-split(data,f=data$Animal) 
datlist<-lapply(datlist,
                FUN=function(df){
                  
                  
                  
                  ### CTRL ###
                  
                  dctrl<-dplyr::filter(df,Odor=="Blank")
                  
                  # Normalizing the coordinates the same way as in Fig2
                  #centering coordinates around the head
                  dctrl$Lant.X<-dctrl$Lant.X-dctrl$head.X
                  dctrl$Rant.X<-dctrl$Rant.X-dctrl$head.X
                  # normalizing by maximum coordinates within trial
                  dctrl$Lant.X.norm<-dctrl$Lant.X/max(c(abs(dctrl$Lant.X),abs(dctrl$Rant.X),na.rm=T))
                  dctrl$Rant.X.norm<-dctrl$Rant.X/max(c(abs(dctrl$Lant.X),abs(dctrl$Rant.X),na.rm=T))
                  
                  #Right
                  dctrl$Rant.X.Max<-MovingFUN(dctrl$Rant.X.norm,win=100,sh=1,FUN="max",new.length=nrow(dctrl)) 
                  dctrl$Rant.X.Min<-MovingFUN(dctrl$Rant.X.norm,win=100,sh=1,FUN="min",new.length=nrow(dctrl))
                  dctrl$Rrange<-dctrl$Rant.X.Max-dctrl$Rant.X.Min
                  dctrl$Rrange<-dctrl$Rrange-mean(dctrl$Rrange[dctrl$Count<=0],na.rm=T) #center
                  #Left
                  dctrl$Lant.X.Max<-MovingFUN(dctrl$Lant.X.norm,win=100,sh=1,FUN="max",new.length=nrow(dctrl)) 
                  dctrl$Lant.X.Min<-MovingFUN(dctrl$Lant.X.norm,win=100,sh=1,FUN="min",new.length=nrow(dctrl))
                  dctrl$Lrange<-dctrl$Lant.X.Max-dctrl$Lant.X.Min
                  dctrl$Lrange<-dctrl$Lrange-mean(dctrl$Lrange[dctrl$Count<=0],na.rm=T) #center
                  
                  #overall range
                  dctrl$overallrange<-MovingRange(dctrl$Lant.X.norm,dctrl$Rant.X.norm,win=100,sh=1,new.length=nrow(dctrl))
                  dctrl$overallrange<-dctrl$overallrange-mean(dctrl$overallrange[dctrl$Count<=0],na.rm=T) #center
                  
                  
                  
                  
                  ### ODOR ###
                  
                  dod<-dplyr::filter(df,Odor=="Colony")
                  
                  # Normalizing the coordinates the same way as in Fig2
                  #centering coordinates around the head
                  dod$Lant.X<-dod$Lant.X-dod$head.X
                  dod$Rant.X<-dod$Rant.X-dod$head.X
                  # normalizing by maximum coordinates within trial
                  dod$Lant.X.norm<-dod$Lant.X/max(c(abs(dod$Lant.X),abs(dod$Rant.X),na.rm=T))
                  dod$Rant.X.norm<-dod$Rant.X/max(c(abs(dod$Lant.X),abs(dod$Rant.X),na.rm=T))
                  
                  #Right
                  dod$Rant.X.Max<-MovingFUN(dod$Rant.X.norm,win=100,sh=1,FUN="max",new.length=nrow(dod)) 
                  dod$Rant.X.Min<-MovingFUN(dod$Rant.X.norm,win=100,sh=1,FUN="min",new.length=nrow(dod))
                  dod$Rrange<-dod$Rant.X.Max-dod$Rant.X.Min
                  dod$Rrange<-dod$Rrange-mean(dod$Rrange[dod$Count<=0],na.rm=T) #center
                  #Left
                  dod$Lant.X.Max<-MovingFUN(dod$Lant.X.norm,win=100,sh=1,FUN="max",new.length=nrow(dod)) 
                  dod$Lant.X.Min<-MovingFUN(dod$Lant.X.norm,win=100,sh=1,FUN="min",new.length=nrow(dod))
                  dod$Lrange<-dod$Lant.X.Max-dod$Lant.X.Min
                  dod$Lrange<-dod$Lrange-mean(dod$Lrange[dod$Count<=0],na.rm=T) #center
                  
                  #overall range
                  dod$overallrange<-MovingRange(dod$Lant.X.norm,dod$Rant.X.norm,win=100,sh=1,new.length=nrow(dod))
                  dod$overallrange<-dod$overallrange-mean(dod$overallrange[dod$Count<=0],na.rm=T) #center
                  
                  #ctrl correction
                  dod$Rrange.correct<-dod$Rrange-dctrl$Rrange
                  dod$Lrange.correct<-dod$Lrange-dctrl$Lrange
                  dod$overallrange.correct<-dod$overallrange-dctrl$overallrange
                  
                  
                  
                  return(dod)
                  
                })
data.ts<-dplyr::bind_rows(datlist) #put df back together



### For static P1, P2, P3 and each antenna
datlist<-split(data,f=data$Animal) 
datlist<-lapply(datlist,
                FUN=function(df){
                  
                  
                  ### CTRL ###
                  
                  d<-dplyr::filter(df,Odor=="Blank")
                  
                  # Normalizing the coordinates the same way as in Fig2
                  #centering coordinates around the head
                  d$Lant.X<-d$Lant.X-d$head.X
                  d$Rant.X<-d$Rant.X-d$head.X
                  # normalizing by maximum coordinates within trial
                  d$Lant.X.norm<-d$Lant.X/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                  d$Rant.X.norm<-d$Rant.X/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                  
                  #Right
                  d$Rant.X.Max<-MovingFUN(d$Rant.X.norm,win=100,sh=1,FUN="max",new.length=nrow(d)) #maximum per 50fr window
                  d$Rant.X.Min<-MovingFUN(d$Rant.X.norm,win=100,sh=1,FUN="min",new.length=nrow(d))
                  d$Rrange<-d$Rant.X.Max-d$Rant.X.Min
                  #Left
                  d$Lant.X.Max<-MovingFUN(d$Lant.X.norm,win=100,sh=1,FUN="max",new.length=nrow(d)) #maximum per 50fr window
                  d$Lant.X.Min<-MovingFUN(d$Lant.X.norm,win=100,sh=1,FUN="min",new.length=nrow(d))
                  d$Lrange<-d$Lant.X.Max-d$Lant.X.Min
                  #overall range
                  d$overallrange<-MovingRange(d$Lant.X.norm,d$Rant.X.norm,win=100,sh=1,new.length=nrow(d))
                  d$overallrange<-d$overallrange-mean(d$overallrange[d$Count<=0],na.rm=T) #center
                  
                  
                  # mean per static
                  Rrange.PS1.ctrl<-mean(d[d$Count<=-1,"Rrange"],na.rm=T)
                  Rrange.PS2.ctrl<-mean(d[d$Count>-1 & d$Count<=2,"Rrange"],na.rm=T)
                  Rrange.PS3.ctrl<-mean(d[d$Count>-2 & d$Count<=0,"Rrange"],na.rm=T)
                  Rrange.static1.ctrl<-mean(d[d$Count>0 & d$Count<=2,"Rrange"],na.rm=T)
                  Rrange.static2.ctrl<-mean(d[d$Count>5 & d$Count<=7,"Rrange"],na.rm=T)
                  Rrange.static3.ctrl<-mean(d[d$Count>10 & d$Count<=12,"Rrange"],na.rm=T)
                  
                  Lrange.PS1.ctrl<-mean(d[d$Count<=-1,"Lrange"],na.rm=T)
                  Lrange.PS2.ctrl<-mean(d[d$Count>-1 & d$Count<=2,"Lrange"],na.rm=T)
                  Lrange.PS3.ctrl<-mean(d[d$Count>-2 & d$Count<=0,"Lrange"],na.rm=T)
                  Lrange.static1.ctrl<-mean(d[d$Count>0 & d$Count<=2,"Lrange"],na.rm=T)
                  Lrange.static2.ctrl<-mean(d[d$Count>5 & d$Count<=7,"Lrange"],na.rm=T)
                  Lrange.static3.ctrl<-mean(d[d$Count>10 & d$Count<=12,"Lrange"],na.rm=T)
                  
                  Overallrange.PS1.ctrl<-mean(d[d$Count<=-1,"overallrange"],na.rm=T)
                  Overallrange.PS2.ctrl<-mean(d[d$Count>-1 & d$Count<=2,"overallrange"],na.rm=T)
                  Overallrange.PS3.ctrl<-mean(d[d$Count>-2 & d$Count<=0,"overallrange"],na.rm=T)
                  Overallrange.static1.ctrl<-mean(d[d$Count>0 & d$Count<=2,"overallrange"],na.rm=T)
                  Overallrange.static2.ctrl<-mean(d[d$Count>5 & d$Count<=7,"overallrange"],na.rm=T)
                  Overallrange.static3.ctrl<-mean(d[d$Count>10 & d$Count<=12,"overallrange"],na.rm=T)
                  
                  # delta
                  Rrange.baseline.ctrl.delta<-Rrange.PS2.ctrl-Rrange.PS1.ctrl #baseline variation before the shift
                  Rrange.static1.ctrl.delta<-Rrange.static1.ctrl-Rrange.PS3.ctrl #delta with only 2s PS
                  Rrange.static2.ctrl.delta<-Rrange.static2.ctrl-Rrange.static1.ctrl
                  Rrange.static3.ctrl.delta<-Rrange.static3.ctrl-Rrange.static1.ctrl
                  
                  Lrange.baseline.ctrl.delta<-Lrange.PS2.ctrl-Lrange.PS1.ctrl #baseline variation before the shift
                  Lrange.static1.ctrl.delta<-Lrange.static1.ctrl-Lrange.PS3.ctrl #delta with only 2s PS
                  Lrange.static2.ctrl.delta<-Lrange.static2.ctrl-Lrange.static1.ctrl
                  Lrange.static3.ctrl.delta<-Lrange.static3.ctrl-Lrange.static1.ctrl
                  
                  Overallrange.baseline.ctrl.delta<-Overallrange.PS2.ctrl-Overallrange.PS1.ctrl #baseline variation before the shift
                  Overallrange.static1.ctrl.delta<-Overallrange.static1.ctrl-Overallrange.PS3.ctrl #delta with only 2s PS
                  Overallrange.static2.ctrl.delta<-Overallrange.static2.ctrl-Overallrange.static1.ctrl
                  Overallrange.static3.ctrl.delta<-Overallrange.static3.ctrl-Overallrange.static1.ctrl
                  
                  
                  
                  
                  ### ODOR ###
                  
                  dflist<-split(df,f=df$Odor)
                  dflist<-lapply(dflist,
                                 FUN=function(d){
                                   
                                   # Normalizing the coordinates the same way as in Fig2
                                   #centering coordinates around the head
                                   d$Lant.X<-d$Lant.X-d$head.X
                                   d$Rant.X<-d$Rant.X-d$head.X
                                   # normalizing by maximum coordinates within trial
                                   d$Lant.X.norm<-d$Lant.X/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                                   d$Rant.X.norm<-d$Rant.X/max(c(abs(d$Lant.X),abs(d$Rant.X),na.rm=T))
                                   
                                   #Right
                                   d$Rant.X.Max<-MovingFUN(d$Rant.X.norm,win=100,sh=1,FUN="max",new.length=nrow(d)) #maximum per 50fr window
                                   d$Rant.X.Min<-MovingFUN(d$Rant.X.norm,win=100,sh=1,FUN="min",new.length=nrow(d))
                                   d$Rrange<-d$Rant.X.Max-d$Rant.X.Min
                                   #Left
                                   d$Lant.X.Max<-MovingFUN(d$Lant.X.norm,win=100,sh=1,FUN="max",new.length=nrow(d)) #maximum per 50fr window
                                   d$Lant.X.Min<-MovingFUN(d$Lant.X.norm,win=100,sh=1,FUN="min",new.length=nrow(d))
                                   d$Lrange<-d$Lant.X.Max-d$Lant.X.Min
                                   #overall range
                                   d$overallrange<-MovingRange(d$Lant.X.norm,d$Rant.X.norm,win=100,sh=1,new.length=nrow(d))
                                   d$overallrange<-d$overallrange-mean(d$overallrange[d$Count<=0],na.rm=T) #center
                                   
                                   # mean per static
                                   Rrange.PS1<-mean(d[d$Count<=-2,"Rrange"],na.rm=T)
                                   Rrange.PS2<-mean(d[d$Count>-2 & d$Count<=0,"Rrange"],na.rm=T)
                                   Rrange.static1<-mean(d[d$Count>0 & d$Count<=2,"Rrange"],na.rm=T)
                                   Rrange.static2<-mean(d[d$Count>5 & d$Count<=7,"Rrange"],na.rm=T)
                                   Rrange.static3<-mean(d[d$Count>10 & d$Count<=12,"Rrange"],na.rm=T)
                                   
                                   Lrange.PS1<-mean(d[d$Count<=-2,"Lrange"],na.rm=T)
                                   Lrange.PS2<-mean(d[d$Count>-2 & d$Count<=0,"Lrange"],na.rm=T)
                                   Lrange.static1<-mean(d[d$Count>0 & d$Count<=2,"Lrange"],na.rm=T)
                                   Lrange.static2<-mean(d[d$Count>5 & d$Count<=7,"Lrange"],na.rm=T)
                                   Lrange.static3<-mean(d[d$Count>10 & d$Count<=12,"Lrange"],na.rm=T)
                                   
                                   Overallrange.PS1<-mean(d[d$Count<=-1,"overallrange"],na.rm=T)
                                   Overallrange.PS2<-mean(d[d$Count>-1 & d$Count<=2,"overallrange"],na.rm=T)
                                   Overallrange.PS3<-mean(d[d$Count>-2 & d$Count<=0,"overallrange"],na.rm=T)
                                   Overallrange.static1<-mean(d[d$Count>0 & d$Count<=2,"overallrange"],na.rm=T)
                                   Overallrange.static2<-mean(d[d$Count>5 & d$Count<=7,"overallrange"],na.rm=T)
                                   Overallrange.static3<-mean(d[d$Count>10 & d$Count<=12,"overallrange"],na.rm=T)
                                   
                                   # delta
                                   Rrange.PS.delta<-Rrange.PS2-Rrange.PS1
                                   Rrange.static1.delta<-Rrange.static1-Rrange.PS2
                                   Rrange.static2.delta<-Rrange.static2-Rrange.PS2
                                   Rrange.static3.delta<-Rrange.static3-Rrange.PS2
                                   
                                   Lrange.PS.delta<-Lrange.PS2-Lrange.PS1
                                   Lrange.static1.delta<-Lrange.static1-Lrange.PS2
                                   Lrange.static2.delta<-Lrange.static2-Lrange.PS2
                                   Lrange.static3.delta<-Lrange.static3-Lrange.PS2
                                   
                                   Overallrange.PS.delta<-Overallrange.PS2-Overallrange.PS1
                                   Overallrange.static1.delta<-Overallrange.static1-Overallrange.PS2
                                   Overallrange.static2.delta<-Overallrange.static2-Overallrange.PS2
                                   Overallrange.static3.delta<-Overallrange.static3-Overallrange.PS2
                                   
                                   
                                   # control correction
                                   Rrange.static1.delta.corrected<-Rrange.static1.delta-Rrange.baseline.ctrl.delta
                                   Rrange.static2.delta.corrected<-Rrange.static2.delta-Rrange.static2.ctrl.delta
                                   Rrange.static3.delta.corrected<-Rrange.static3.delta-Rrange.static3.ctrl.delta
                                   
                                   Lrange.static1.delta.corrected<-Lrange.static1.delta-Lrange.baseline.ctrl.delta
                                   Lrange.static2.delta.corrected<-Lrange.static2.delta-Lrange.static2.ctrl.delta
                                   Lrange.static3.delta.corrected<-Lrange.static3.delta-Lrange.static3.ctrl.delta
                                   
                                   Overallrange.static1.delta.corrected<-Overallrange.static1.delta-Overallrange.baseline.ctrl.delta
                                   Overallrange.static2.delta.corrected<-Overallrange.static2.delta-Overallrange.static2.ctrl.delta
                                   Overallrange.static3.delta.corrected<-Overallrange.static3.delta-Overallrange.static3.ctrl.delta
                                   
                                   
                                   
                                   # avg of all windows after correction
                                   overallrange.change.corrected.avg<-mean(c(Overallrange.static1.delta.corrected,
                                                                             Overallrange.static2.delta.corrected,
                                                                             Overallrange.static3.delta.corrected))
                                   
                                   
                                   
                                   
                                   return(data.frame("Animal"=rep(unique(d$Animal),6),
                                                     "Test"=rep(unique(d$Test),6),
                                                     "Odor"=rep(unique(d$Odor),6),
                                                     "Ant"=rep(c("L","R"),each=3),
                                                     "window"=rep(c("static1","static2","static3"),2),
                                                     "range.change.corrected"=c(Lrange.static1.delta.corrected,
                                                                                Lrange.static2.delta.corrected,
                                                                                Lrange.static3.delta.corrected,
                                                                                Rrange.static1.delta.corrected,
                                                                                Rrange.static2.delta.corrected,
                                                                                Rrange.static3.delta.corrected),
                                                     "overallrange.change"=c(Overallrange.static1.delta,
                                                                             Overallrange.static2.delta,
                                                                             Overallrange.static3.delta,
                                                                             Overallrange.static1.delta,
                                                                             Overallrange.static2.delta,
                                                                             Overallrange.static3.delta),
                                                     "overallrange.change.corrected.avg"=rep(overallrange.change.corrected.avg,6)
                                   )
                                   )
                                   
                                   
                                   
                                 })
                  df<-dplyr::bind_rows(dflist) #put df back together
                  
                  
                  
                  return(df)
                  
                  
                  
                  
                })
data<-dplyr::bind_rows(datlist) #put df back together

data.model<-data








### Plotting the time-series 


#data
n<-length(unique(data.ts$Test))

Lsem<-aggregate(Lrange.correct~Count,data=data.ts,FUN="sd",na.rm=T)
Lsem$Lrange.correct<-Lsem$Lrange.correct/sqrt(n)
Lsem$smooth<-MovingFUN(Lsem$Lrange.correct,100,1,FUN="mean",new.length=nrow(Lsem))
Lsem$smooth<-repeat.before(Lsem$smooth)

Rsem<-aggregate(Rrange.correct~Count,data=data.ts,FUN="sd",na.rm=T)
Rsem$Rrange.correct<-Rsem$Rrange.correct/sqrt(n)
Rsem$smooth<-MovingFUN(Rsem$Rrange.correct,100,1,FUN="mean",new.length=nrow(Rsem))
Rsem$smooth<-repeat.before(Rsem$smooth)

data.ts.l<-aggregate(Lrange.correct~Count,data=data.ts,FUN="mean",na.rm=T)
data.ts.l$smooth<-MovingFUN(data.ts.l$Lrange.correct,100,1,FUN="mean",new.length=nrow(data.ts.l))
data.ts.l$smooth<-repeat.before(data.ts.l$smooth)

data.ts.r<-aggregate(Rrange.correct~Count,data=data.ts,FUN="mean",na.rm=T)
data.ts.r$smooth<-MovingFUN(data.ts.r$Rrange.correct,100,1,FUN="mean",new.length=nrow(data.ts.r))
data.ts.r$smooth<-repeat.before(data.ts.r$smooth)




pdf(file=paste("./Plots/AntRange_withintrialctrlcorrected_MovingStream_Colony.pdf",sep=""),
    width=10,height=8,bg="white",useDingbats=F)


ylim<-c(-0.3,0.3)


par(mar=c(5,5,4,2)+0.1)

plot(data.ts.r$Count,data.ts.r$Rrange.correct-200,type="l",col="dodgerblue3",
     main=paste(Od[1],"- Average corrected antennal range"),
     ylab="Corrected range [mm]",xlab="Time from odor onset [s]",
     ylim=ylim,
     xlim=c(-2,14),
     xaxs="i",yaxs="i",
     cex.main=1,font.main=2,cex.axis=1,font.axis=2,cex.lab=1,font.lab=2)
abline(h=0,lty=2)
#L
polygon(x=c(data.ts.l$Count,rev(data.ts.l$Count)),
        y=c(data.ts.l$Lrange.correct+Lsem$Lrange.correct,rev(data.ts.l$Lrange.correct-Lsem$Lrange.correct)),  #sem
        col=rgb(firebrick3[1],firebrick3[2],firebrick3[3],alpha=70,names=NULL,maxColorValue=255),border=NA)
lines(data.ts.l$Count,data.ts.l$Lrange.correct,col="firebrick3",lwd=4)
#R
polygon(x=c(data.ts.r$Count,rev(data.ts.r$Count)),
        y=c(data.ts.r$Rrange.correct+Rsem$Rrange.correct,rev(data.ts.r$Rrange.correct-Rsem$Rrange.correct)),  #sem
        col=rgb(dodgerblue3[1],dodgerblue3[2],dodgerblue3[3],alpha=70,names=NULL,maxColorValue=255),border=NA)
lines(data.ts.r$Count,data.ts.r$Rrange.correct,col="dodgerblue3",lwd=4)


abline(v=c(0,12))



dev.off()



#####









#############
# Modelling #
#############

# data
data.model.col<-dplyr::filter(data.model,Odor=="Colony")
hist(data.model.col$range.change.corrected)
# for the overall range per position for colony
data.model.overallrange<-dplyr::filter(data.model,Odor=="Colony",Ant=="L")
data.model.overallrange<-data.model.overallrange[!duplicated(data.model.overallrange), ] #remove duplicates
hist(data.model.overallrange$overallrange.change)
# for the avg overall of an odor
data.model.odoravg<-dplyr::filter(data.model,Ant=="L")
data.model.odoravg<-data.model.odoravg[,c("Animal","Odor","overallrange.change","overallrange.change.corrected.avg")]
data.model.odoravg<-data.model.odoravg[!duplicated(data.model.odoravg), ] #remove duplicates
hist(data.model.odoravg$overallrange.change.corrected.avg)



# Fit
set.seed(090421)
brmsmod<-brm( range.change.corrected ~ Ant + window + Ant:window + (1|Animal), #for Colony antennal range
              data=data.model.col, #for Colony antennal range
              
             # overallrange.change ~ window + (1|Animal), #for Colony overall range
             # data=data.model.overallrange, #for Colony overall range
             
             # overallrange.change.corrected.avg ~ Odor + (1|Animal), #for the overall range avg per odor
             # data=data.model.odoravg, #for for the overall avg per odor
             
             family=student(),
             iter=10000, chains=4, cores=4, 
             backend="cmdstanr", control=list(adapt_delta=0.99,max_treedepth=20))




# pp check
pp<-pp_check(brmsmod,ndraws=200) + xlim(-2,2)
pp


#convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce



# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))



# Estimates, 95% CrI, and probabilities

# for Colony left & right antennae ranges
#L
ColP1L.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ColP1L.quants
print(paste("P(ColP1L>0) = ",round(mean((simval[,1])>0),2),sep=""))

ColP2L.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
ColP2L.quants
print(paste("P(ColP2L>0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))

ColP3L.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
ColP3L.quants
print(paste("P(ColP3L>0) = ",round(mean((simval[,1]+simval[,4])>0),2),sep=""))

#R
ColP1R.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColP1R.quants
print(paste("P(ColP1R>0) = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))

ColP2R.quants<-round(quantile((simval[,1]+simval[,2]+simval[,3]+simval[,5]),probs=c(0.5,0.025,0.975)),2)
ColP2R.quants
print(paste("P(ColP2R>0) = ",round(mean((simval[,1]+simval[,2]+simval[,3]+simval[,5])>0),2),sep=""))

ColP3R.quants<-round(quantile((simval[,1]+simval[,2]+simval[,4]+simval[,6]),probs=c(0.5,0.025,0.975)),2)
ColP3R.quants
print(paste("P(ColP3R>0) = ",round(mean((simval[,1]+simval[,2]+simval[,4]+simval[,6])>0),2),sep=""))




# for Colony overall range
ColP1.quants<-round(quantile((simval[,1]),probs=c(0.5,0.025,0.975)),2)
ColP1.quants
print(paste("P(ColP1>0) = ",round(mean((simval[,1])>0),2),sep=""))

ColP2.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
ColP2.quants
print(paste("P(ColP2>0) = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))

ColP3.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
ColP3.quants
print(paste("P(ColP3>0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))


# for overall average per odor
Col.quants<-round(quantile((simval[,1]+simval[,3]),probs=c(0.5,0.025,0.975)),2)
Col.quants
print(paste("P(Col>0) = ",round(mean((simval[,1]+simval[,3])>0),2),sep=""))

But.quants<-round(quantile((simval[,1]+simval[,2]),probs=c(0.5,0.025,0.975)),2)
But.quants
print(paste("P(But>0) = ",round(mean((simval[,1]+simval[,2])>0),2),sep=""))

Lin.quants<-round(quantile((simval[,1]+simval[,4]),probs=c(0.5,0.025,0.975)),2)
Lin.quants
print(paste("P(Lin>0) = ",round(mean((simval[,1]+simval[,4])>0),2),sep=""))







# Plotting effect plots

pdf(file=paste("./Plots/RangeChanges_ctrlcorrectionwithinAnimal_MovingStream_EffectBarplot.pdf",sep="_"),
    width=8,height=7,bg="white",useDingbats=F)


par(mfrow=c(1,3))


par(lwd=5)


#L&R
m.fit<-matrix(ncol=2,nrow=3)
colnames(m.fit)<-c("L","R")
rownames(m.fit)<-c("P1","P2","P3")
m.fit[,1]<-c(ColP1L.quants[1],ColP2L.quants[1],ColP3L.quants[1])
m.fit[,2]<-c(ColP1R.quants[1],ColP2R.quants[1],ColP3R.quants[1])

m.lwr<-matrix(ncol=2,nrow=3)
colnames(m.lwr)<-c("L","R")
rownames(m.lwr)<-c("P1","P2","P3")
m.lwr[,1]<-c(ColP1L.quants[2],ColP2L.quants[2],ColP3L.quants[2])
m.lwr[,2]<-c(ColP1R.quants[2],ColP2R.quants[2],ColP3R.quants[2])

m.upr<-matrix(ncol=2,nrow=3)
colnames(m.upr)<-c("L","R")
rownames(m.upr)<-c("P1","P2","P3")
m.upr[,1]<-c(ColP1L.quants[3],ColP2L.quants[3],ColP3L.quants[3])
m.upr[,2]<-c(ColP1R.quants[3],ColP2R.quants[3],ColP3R.quants[3])



ylim=c(-0.25,0.35)
barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="",main="Corrected change 
in Range - Colony")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.01,2),xpd=T,labels=rep(c("P1","P2","P3"),2),cex=0.7)
text(x=barcenters[2,c(1,2)],y=rep(-0.6,2),xpd=T,labels=c("L","R"),cex=0.7)



#avg per position
m.fit<-matrix(ncol=1,nrow=3)
rownames(m.fit)<-c("P1","P2","P3")
m.fit[,1]<-c(ColP1.quants[1],ColP2.quants[1],ColP3.quants[1])

m.lwr<-matrix(ncol=1,nrow=3)
rownames(m.lwr)<-c("P1","P2","P3")
m.lwr[,1]<-c(ColP1.quants[2],ColP2.quants[2],ColP3.quants[2])

m.upr<-matrix(ncol=1,nrow=3)
rownames(m.upr)<-c("P1","P2","P3")
m.upr[,1]<-c(ColP1.quants[3],ColP2.quants[3],ColP3.quants[3])


barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="",main="Change in overall range
per position - Colony")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.01,2),xpd=T,labels=rep(c("P1","P2","P3"),2),cex=0.7)



#avg per odor
m.fit<-matrix(ncol=1,nrow=3)
rownames(m.fit)<-c("Col","But","Lin")
m.fit[,1]<-c(Col.quants[1],But.quants[1],Lin.quants[1])

m.lwr<-matrix(ncol=1,nrow=3)
rownames(m.lwr)<-c("Col","But","Lin")
m.lwr[,1]<-c(Col.quants[2],But.quants[2],Lin.quants[2])

m.upr<-matrix(ncol=1,nrow=3)
rownames(m.upr)<-c("Col","But","Lin")
m.upr[,1]<-c(Col.quants[3],But.quants[3],Lin.quants[3])


barcenters<-barplot(height=m.fit,beside=T,plot=T,
                    xaxt="n",ylim=ylim,
                    cex.axis=2,ylab="",main="Avg Change in overall range
per odor")
segments(x0=barcenters,y0=m.lwr,x1=barcenters,y1=m.upr,col="black")
text(x=barcenters,y=rep(par("usr")[3]-0.01,2),xpd=T,labels=rep(c("Col","But","Lin"),2),cex=0.7)



dev.off()












