###########################################################################################################################
#                                   Active smelling in the American cockroach: Fig. 6                                     #
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
require("brms")
require("cmdstanr")
require("bayesplot")
require("ggpubr")
require("fitdistrplus")
require("imager")
require("abind")


#colors
gray60<-as.vector(col2rgb(col =c("gray60"),alpha=F ))
gray40<-as.vector(col2rgb(col =c("gray40"),alpha=F ))




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
      
      index[j]<-(i)  
      i=i+shift
      j=j+1
      
    }
    
  } else if(method=="back"){
    
    i=winsize
    j=1
    while(i<=length(x)){
      
      if(FUN=="mean"){
        fun.out[j]<-mean(x[(i-(winsize-1)):i])
      } else if(FUN=="sum"){
        fun.out[j]<-sum(x[(i-(winsize-1)):i])
      } else if(FUN=="min"){
        fun.out[j]<-min(x[(i-(winsize-1)):i])
      } else if(FUN=="max"){
        fun.out[j]<-max(x[(i-(winsize-1)):i])
      } else if(FUN=="median"){
        fun.out[j]<-median(x[(i-(winsize-1)):i])
      } else if(FUN=="sd"){
        fun.out[j]<-sd(x[(i-(winsize-1)):i])
      }
      
      if(is.nan(fun.out[j])==T){ fun.out[j]<-NA }
      
      index[j]<-(i)
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
    
    m<-mean(var[(d[i]-halfwin):(d[i]+halfwin)],na.rm=T) # surrounding average (+- halfwin)
    
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
    
    if(sum(peakvect)==0){ #if there are no peaks 
      
      return(0)
      
    } else if(sum(peakvect)==1){
      
      return(0)
      
    } else {
      
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
}













#--------------------------------------------------------------------------------------------------------------------------------
# Smoke Distribution: Horizontal sweep - selected frames
#--------------------------------------------------------------------------------------------------------------------------------




#####################################
### Import and pre-process frames ###
#####################################

png.files<-list.files(path = "../Data/SmokeData_selectedframes/horizontal/selected_png/sequence")  

horizontal_pnglist<-list()
for (i in 1:length(png.files)){  
  
  png.name<-png.files[i]
  
  png.image<-imager::load.image(paste("../Data/SmokeData_selectedframes/horizontal/selected_png/sequence/",png.name,sep="")) 
  
  
  #extract pixel values
  frame<-png.image
  frame<-imager::grayscale(frame)
  frame<-as.data.frame(frame)
  #create a matrix with pixel values
  frame.mat<-matrix(data=NA,nrow=max(frame$x),ncol=max(frame$y))
  frame.mat[cbind(frame$x,rev(frame$y))]<-frame$value
  #crop
  frame.mat<-frame.mat[800:1500,500:825]
  
  
  horizontal_pnglist[[i]]<-frame.mat
  names(horizontal_pnglist)[i]<-paste("png",png.name,sep="_")
  
}



### Mosaic ###

png(filename="./Plots/horizontal_mosaic_raw.png",width=12,height=4,units="in",bg="white",res=300)


par(mfrow=c(3,3),
    mar=c(0,0.2,1,0.2))

colfunc<-colorRampPalette(c("black","white"))

for(i in 1:length(horizontal_pnglist)){
  
  frame<-horizontal_pnglist[[i]]
  name<-names(horizontal_pnglist)[i]
  
  #plot
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(20),
        axes=F,main=name)
  
  #scale (105px = 10mm)
  segments(x0=602.5,x1=655,y0=50,y1=50,
           lwd=3,lend="square",col="white")
  segments(x0=655,x1=655,y0=50,y1=102.5,
           lwd=3,lend="square",col="white")
  
  
}


dev.off()












#############################
#### Frames Thresholding ####
#############################


### Getting frame thresholding ###

horizontal_framethresh<-list()
for(i in 1:length(horizontal_pnglist)){
  
  current<-horizontal_pnglist[[i]]
  
  #thresholding
  frame.quants<-quantile(current,probs=0.75)
  current[current<=frame.quants[1]]<-0
  current[current>frame.quants[1]]<-1
  
  #save
  horizontal_framethresh[[i]]<-current
  
}





### Plot all frames ###

png(filename="./Plots/horizontal_frames_thresh75.png",width=12,height=6,units="in",bg="white",res=300)


par(mfrow=c(3,3),
    mar=c(0,0.2,1,0.2))

colfunc<-colorRampPalette(c("black","white"))

for(i in 1:length(horizontal_framethresh)){
  
  frame<-horizontal_framethresh[[i]]
  
  #plot
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(20),
        axes=F,main=i)
  
}


dev.off()





### Save processed frames as png ###

for(i in 1:length(horizontal_framethresh)){
  
  frame<-horizontal_framethresh[[i]]
  
  #flip around
  frame<-apply(frame,2,rev)
  frame<-apply(frame,1,rev)
  frame<-t(frame)
  frame<-apply(frame,2,rev)
  frame.cimg<-imager::as.cimg(frame)
  
  #make frame counter for ffmpeg
  if(i<10){
    n<-paste("00",i,sep="")
  } else if(i>=10 & i<100){
    n<-paste("0",i,sep="")
  } else if(i>=100 & i<1000){
    n<-i
  }
  
  #save
  imager::save.image(frame.cimg,
                     paste("../Data/SmokeData_selectedframes/horizontal/processed_png/thresholded75/unmasked/",n,".png",sep=""))
  
}












############################
#### Pixel distribution ####
############################
#from the processed, thresholded frames with erased animal in photoshop



### Import selected processed masked frames ###

processed.files<-list.files(path = "../Data/SmokeData_selectedframes/horizontal/processed_png/thresholded75/masked")

horizontal_processedlist<-list()
for (i in 1:length(processed.files)){  
  
  name<-processed.files[i]
  
  processed.image<-imager::load.image(paste("../Data/SmokeData_selectedframes/horizontal/processed_png/thresholded75/masked/",name,sep=""))
  
  
  #extract pixel values
  frame<-processed.image
  frame.mat<-frame[,,,1]
  frame.mat<-frame.mat[,rev(1:ncol(frame.mat))]#flip vertically
  # make sure there are only 1s and 0s
  frame.mat[frame.mat!=0]<-1
  
  horizontal_processedlist[[i]]<-frame.mat
  names(horizontal_processedlist)[i]<-name
  
}





### plot frames ###

png(filename="./Plots/horizontal_thresh75_masked_ROI.png",width=12,height=3,units="in",bg="white",res=300)

par(mfrow=c(1,3))

for(i in 1:length(horizontal_processedlist)){
  
  frame<-horizontal_processedlist[[i]]
  
  #plot
  colfunc<-colorRampPalette(c("black","white"))
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(2))
  
  abline(v=c(360,600),col="red") #ROI
  
}

dev.off()








### Distribution within the ROI ###
# sum of X for each Y value = number of white pixels


ROI<-c(360,600) #limits


horizontal_pixeldistr.list<-list()
for(i in 1:length(horizontal_processedlist)){
  
  frame<-horizontal_processedlist[[i]]
  frame<-frame[ROI[1]:ROI[2],] #subset with ROI
  
  pixel.distr<-apply(frame,2,sum) #margin=2 is columns = Y axis in this case
  pixel.distr<-rev(pixel.distr) #Y axis is flipped in the matrix, I want Y=0 on top
  
  horizontal_pixeldistr.list[[i]]<-pixel.distr
  names(horizontal_pixeldistr.list)[i]<-names(horizontal_processedlist)[i]
  
}






# plot distribution

pdf(file="./Plots/horizontal_distr_widths.pdf",width=10,height=5,bg="white",useDingbats=F)


par(mfrow=c(1,2))

frame.distr<-horizontal_pixeldistr.list[[1]]
plot(y=frame.distr,x=(1:length(frame.distr))/105, #105 is 10mm
     type="l",col="gray",lwd=2,
     main="Width of the smoke - horizontal",
     ylab="number of white pixels",xlab="Y coordinate [mm]")

frame.distr<-horizontal_pixeldistr.list[[2]]
lines(y=frame.distr,x=(1:length(frame.distr))/105,col="gray40",lwd=2)

frame.distr<-horizontal_pixeldistr.list[[3]]
lines(y=frame.distr,x=(1:length(frame.distr))/105,col="black",lwd=2)


# Smoke Width #

widths<-rep(NA,3)

frame1.distr<-horizontal_pixeldistr.list[[1]]
frame1.width<-length(frame1.distr[frame1.distr>0]) #number of Y rows with white pixels in them
frame.width<-frame1.width/frame1.width #Normalise by the starting width
widths[1]<-frame.width

frame.distr<-horizontal_pixeldistr.list[[2]]
frame.width<-length(frame.distr[frame.distr>0]) 
frame.width<-frame.width/frame1.width
widths[2]<-frame.width

frame.distr<-horizontal_pixeldistr.list[[3]]
frame.width<-length(frame.distr[frame.distr>0]) #number of Y rows with white pixels in them
frame.width<-frame.width/frame1.width
widths[3]<-frame.width


#plot widhts
par(lwd=2)
barcenters<-barplot(height=widths,beside=T,plot=T,width=0.5,
                    ylim=c(0,2),xaxt="n",las=1,
                    col=c("gray","gray40","black"),border="black",
                    ylab="",main="Normalised smoke width - horizontal sweep")
axis(1,at=barcenters,labels=c("before","inward","outward"),lty=0,cex.axis=1)


dev.off()

















#--------------------------------------------------------------------------------------------------------------------------------
# Smoke Distribution: Vertical sweep - selected frames
#--------------------------------------------------------------------------------------------------------------------------------



#####################################
### Import and pre-process frames ###
#####################################

png.files<-list.files(path = "../Data/SmokeData_selectedframes/vertical/selected_png/sequence")  

vertical_pnglist<-list()
for (i in 1:length(png.files)){  
  
  png.name<-png.files[i]
  
  png.image<-imager::load.image(paste("../Data/SmokeData_selectedframes/vertical/selected_png/sequence/",png.name,sep="")) 
  
  
  #extract pixel values
  frame<-png.image
  frame<-imager::grayscale(frame)
  frame<-as.data.frame(frame)
  #create a matrix with pixel values
  frame.mat<-matrix(data=NA,nrow=max(frame$x),ncol=max(frame$y))
  frame.mat[cbind(frame$x,rev(frame$y))]<-frame$value
  #crop
  frame.mat<-frame.mat[800:1500,500:825]
  
  
  vertical_pnglist[[i]]<-frame.mat
  names(vertical_pnglist)[i]<-paste("png",png.name,sep="_")
  
}







### Mosaic ###

png(filename="./Plots/vertical_mosaic_raw.png",width=12,height=4,units="in",bg="white",res=300)


par(mfrow=c(3,3),
    mar=c(0,0.2,1,0.2))

colfunc<-colorRampPalette(c("black","white"))

for(i in 1:length(vertical_pnglist)){
  
  frame<-vertical_pnglist[[i]]
  name<-names(vertical_pnglist)[i]
  
  #plot
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(20),
        axes=F,main=name)
  
  #scale (105px = 10mm)
  segments(x0=602.5,x1=655,y0=50,y1=50,
           lwd=3,lend="square",col="white")
  segments(x0=655,x1=655,y0=50,y1=102.5,
           lwd=3,lend="square",col="white")
  
  
}


dev.off()









################################
#### Raw Frames Thresholded ####
################################


### thresholding ###

vertical_framethresh<-list()
for(i in 1:length(vertical_pnglist)){
  
  current<-vertical_pnglist[[i]]
  
  #thresholding
  frame.quants<-quantile(current,probs=0.75)
  current[current<=frame.quants[1]]<-0
  current[current>frame.quants[1]]<-1
  
  #save
  vertical_framethresh[[i]]<-current
  
}





### Plot all frames ###

png(filename=".Plots/vertical_frames_thresh75.png",width=10,height=4,units="in",bg="white",res=300)


par(mfrow=c(3,3),
    mar=c(0,0.2,1,0.2))

colfunc<-colorRampPalette(c("black","white"))

for(i in 1:length(vertical_framethresh)){
  
  frame<-vertical_framethresh[[i]]
  
  #plot
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(20),
        axes=F,main=i)
  
}


dev.off()





### Save processed frames as png ###

for(i in 1:length(vertical_framethresh)){
  
  frame<-vertical_framethresh[[i]]
  
  #flip around
  frame<-apply(frame,2,rev)
  frame<-apply(frame,1,rev)
  frame<-t(frame)
  frame<-apply(frame,2,rev)
  frame.cimg<-imager::as.cimg(frame)
  
  #make frame counter fo ffmpeg
  if(i<10){
    n<-paste("00",i,sep="")
  } else if(i>=10 & i<100){
    n<-paste("0",i,sep="")
  } else if(i>=100 & i<1000){
    n<-i
  }
  
  #save
  imager::save.image(frame.cimg,
                     paste("../Data/SmokeData_selectedframes/vertical/processed_png/thresholded75/unmasked/",n,".png",sep=""))

}











############################
#### Pixel distribution ####
############################
#from the processed, thresholded frames with erased animal in photoshop


### Import processed frames ###

processed.files<-list.files(path = "../Data/SmokeData_selectedframes/vertical/processed_png/thresholded75/masked")

vertical_processedlist<-list()
for (i in 1:length(processed.files)){  
  
  name<-processed.files[i]
  
  processed.image<-imager::load.image(paste("../Data/SmokeData_selectedframes/vertical/processed_png/thresholded75/masked/",name,sep=""))
  
  
  #extract pixel values
  frame<-processed.image
  frame.mat<-frame[,,,1]
  frame.mat<-frame.mat[,rev(1:ncol(frame.mat))]#flip vertically
  # make sure there are only 1s and 0s
  frame.mat[frame.mat!=0]<-1
  
  vertical_processedlist[[i]]<-frame.mat
  names(vertical_processedlist)[i]<-name
  
}





### plot frame ###

png(filename=".Plots/vertical_thresh75_masked_ROI.png",width=12,height=3,units="in",bg="white",res=300)

par(mfrow=c(1,3))

for(i in 1:length(vertical_processedlist)){
  
  frame<-vertical_processedlist[[i]]
  
  #plot
  colfunc<-colorRampPalette(c("black","white"))
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(2))
  
  abline(v=c(360,600),col="red") #ROI
  
}

dev.off()







### Distribution within the ROI ###
# sum of X for each Y value = number of white pixels


ROI<-c(360,600) #limits


vertical_pixeldistr.list<-list()
for(i in 1:length(vertical_processedlist)){
  
  frame<-vertical_processedlist[[i]]
  frame<-frame[ROI[1]:ROI[2],] #subset with ROI
  
  pixel.distr<-apply(frame,2,sum) #margin=2 is columns = Y axis in this case
  pixel.distr<-rev(pixel.distr) #Y axis is flipped in the matrix, I want Y=0 on top
  
  vertical_pixeldistr.list[[i]]<-pixel.distr
  names(vertical_pixeldistr.list)[i]<-names(vertical_processedlist)[i]
  
}




# plot distribution

pdf(file=".Plots/vertical_distr_widths.pdf",width=10,height=5,bg="white",useDingbats=F)


par(mfrow=c(1,2))

frame.distr<-vertical_pixeldistr.list[[1]]
plot(y=frame.distr,x=(1:length(frame.distr))/105, #105 is 10mm,
     type="l",col="gray",lwd=2,
     # xlim=c(0,250),
     main="Width of the smoke - vertical",
     ylab="number of white pixels",xlab="Y coordinate")

frame.distr<-vertical_pixeldistr.list[[2]]
lines(y=frame.distr,x=(1:length(frame.distr))/105,col="gray40",lwd=2)

frame.distr<-vertical_pixeldistr.list[[3]]
lines(y=frame.distr,x=(1:length(frame.distr))/105,col="black",lwd=2)


# Smoke Width #

widths<-rep(NA,3)

frame1.distr<-vertical_pixeldistr.list[[1]]
frame1.width<-length(frame1.distr[frame1.distr>0]) #number of Y rows with white pixels in them
frame.width<-frame1.width/frame1.width #Normalise by the starting width
widths[1]<-frame.width

frame.distr<-vertical_pixeldistr.list[[2]]
frame.width<-length(frame.distr[frame.distr>0]) 
frame.width<-frame.width/frame1.width
widths[2]<-frame.width

frame.distr<-vertical_pixeldistr.list[[3]]
frame.width<-length(frame.distr[frame.distr>0]) #number of Y rows with white pixels in them
frame.width<-frame.width/frame1.width
widths[3]<-frame.width


#plot widhts
par(lwd=2)
barcenters<-barplot(height=widths,beside=T,plot=T,width=0.5,
                    ylim=c(0,2.5),xaxt="n",las=1,
                    col=c("gray","gray40","black"),border="black",
                    ylab="",main="Normalised smoke width - vertical sweep")
axis(1,at=barcenters,labels=c("before","down-sweep","after"),lty=0,cex.axis=1)



dev.off()




















#--------------------------------------------------------------------------------------------------------------------------------
# Smoke Distribution: Control - dead animal 
#--------------------------------------------------------------------------------------------------------------------------------




#####################################
### Import and pre-process frames ###
#####################################

png.files<-list.files(path = "../Data/SmokeData_selectedframes/ctrl/selected_png/sequence") 

ctrl_pnglist<-list()
for (i in 1:length(png.files)){  
  
  png.name<-png.files[i]
  
  png.image<-imager::load.image(paste("../Data/SmokeData_selectedframes/ctrl/selected_png/sequence/",png.name,sep=""))
  
  
  #extract pixel values
  frame<-png.image
  frame<-imager::grayscale(frame)
  frame<-as.data.frame(frame)
  #create a matrix with pixel values
  frame.mat<-matrix(data=NA,nrow=max(frame$x),ncol=max(frame$y))
  frame.mat[cbind(frame$x,rev(frame$y))]<-frame$value
  #crop
  frame.mat<-frame.mat[800:1500,500:825]
  
  
  ctrl_pnglist[[i]]<-frame.mat
  names(ctrl_pnglist)[i]<-paste("png",png.name,sep="_")
  
}








### Mosaic ###

png(filename="./Plots/ctrl_mosaic_raw.png",width=12,height=4,units="in",bg="white",res=300)


par(mfrow=c(3,3),
    mar=c(0,0.2,1,0.2))

colfunc<-colorRampPalette(c("black","white"))

for(i in 1:length(ctrl_pnglist)){
  
  frame<-ctrl_pnglist[[i]]
  name<-names(ctrl_pnglist)[i]
  
  #plot
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(20),
        axes=F,main=name)
  
  #scale (105px = 10mm)
  segments(x0=602.5,x1=655,y0=50,y1=50,
           lwd=3,lend="square",col="white")
  segments(x0=655,x1=655,y0=50,y1=102.5,
           lwd=3,lend="square",col="white")
  
  
}


dev.off()















################################
#### Raw Frames Thresholded ####
################################



### thresholding ###

ctrl_framethresh<-list()
for(i in 1:length(ctrl_pnglist)){
  
  current<-ctrl_pnglist[[i]]
  
  #thresholding
  frame.quants<-quantile(current,probs=0.75)
  current[current<=frame.quants[1]]<-0
  current[current>frame.quants[1]]<-1
  
  #save
  ctrl_framethresh[[i]]<-current
  
}





### Plot all frames ###

png(filename="./Plot/ctrl_frames_thresh75.png",width=10,height=4,units="in",bg="white",res=300)


par(mfrow=c(3,3),
    mar=c(0,0.2,1,0.2))

colfunc<-colorRampPalette(c("black","white"))

for(i in 1:length(ctrl_framethresh)){
  
  frame<-ctrl_framethresh[[i]]
  
  #plot
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(20),
        axes=F,main=i)
  
  
  # abline(v=c(360,600),col="red") #ROI
  
}


dev.off()





### Save processed frames as png ###

for(i in 1:length(ctrl_framethresh)){
  
  frame<-ctrl_framethresh[[i]]
  
  #flip around
  frame<-apply(frame,2,rev)
  frame<-apply(frame,1,rev)
  frame<-t(frame)
  frame<-apply(frame,2,rev)
  frame.cimg<-imager::as.cimg(frame)
  
  #make frame counter fo ffmpeg
  if(i<10){
    n<-paste("00",i,sep="")
  } else if(i>=10 & i<100){
    n<-paste("0",i,sep="")
  } else if(i>=100 & i<1000){
    n<-i
  }
  
  #save
  imager::save.image(frame.cimg,
                     paste("../Data/SmokeData_selectedframes/ctrl/processed_png/thresholded75/unmasked/",n,".png",sep=""))
  
}











############################
#### Pixel distribution ####
############################
#from the processed, thresholded frames with erased animal in photoshop



### Import processed frames ###

processed.files<-list.files(path = "../Data/SmokeData_selectedframes/ctrl/processed_png/thresholded75/masked")

ctrl_processedlist<-list()
for (i in 1:length(processed.files)){  
  
  name<-processed.files[i]
  
  processed.image<-imager::load.image(paste("../Data/SmokeData_selectedframes/ctrl/processed_png/thresholded75/masked/",name,sep=""))
  
  
  #extract pixel values
  frame<-processed.image
  frame.mat<-frame[,,,1]
  frame.mat<-frame.mat[,rev(1:ncol(frame.mat))]#flip vertically
  # make sure there are only 1s and 0s
  frame.mat[frame.mat!=0]<-1
  
  ctrl_processedlist[[i]]<-frame.mat
  names(ctrl_processedlist)[i]<-name
  
}





### plot frame ###

png(filename="ctrl_thresh75_masked_ROI.png",width=12,height=3,units="in",bg="white",res=300)

par(mfrow=c(1,3))

for(i in 1:length(ctrl_processedlist)){
  
  frame<-ctrl_processedlist[[i]]
  
  #plot
  colfunc<-colorRampPalette(c("black","white"))
  image(1:nrow(frame),1:ncol(frame),frame,col=colfunc(2))
  
  abline(v=c(360,600),col="red") #ROI
  
}

dev.off()







### Distribution within the ROI ###
# sum of X for each Y value = number of white pixels


ROI<-c(360,600) #limits


ctrl_pixeldistr.list<-list()
for(i in 1:length(ctrl_processedlist)){
  
  frame<-ctrl_processedlist[[i]]
  frame<-frame[ROI[1]:ROI[2],] #subset with ROI
  
  pixel.distr<-apply(frame,2,sum) #margin=2 is columns = Y axis in this case
  pixel.distr<-rev(pixel.distr) #Y axis is flipped in the matrix, I want Y=0 on top
  
  ctrl_pixeldistr.list[[i]]<-pixel.distr
  names(ctrl_pixeldistr.list)[i]<-names(ctrl_processedlist)[i]
  
}




# plot distribution

pdf(file="./Plots/ctrl_distr_widths.pdf",width=10,height=5,bg="white",useDingbats=F)


par(mfrow=c(1,2))

frame.distr<-ctrl_pixeldistr.list[[1]]
plot(y=frame.distr,x=(1:length(frame.distr))/105,
     type="l",col="gray",lwd=2,
     # xlim=c(0,250),
     main="Width of the smoke - ctrl",
     ylab="number of white pixels",xlab="Y coordinate")

frame.distr<-ctrl_pixeldistr.list[[2]]
lines(y=frame.distr,x=(1:length(frame.distr))/105,col="gray40",lwd=2)

frame.distr<-ctrl_pixeldistr.list[[3]]
lines(y=frame.distr,x=(1:length(frame.distr))/105,col="black",lwd=2)


# Smoke Width #

widths<-rep(NA,3)

frame1.distr<-ctrl_pixeldistr.list[[1]]
frame1.width<-length(frame1.distr[frame1.distr>0]) #number of Y rows with white pixels in them
frame.width<-frame1.width/frame1.width #Normalise by the starting width
widths[1]<-frame.width

frame.distr<-ctrl_pixeldistr.list[[2]]
frame.width<-length(frame.distr[frame.distr>0]) 
frame.width<-frame.width/frame1.width
widths[2]<-frame.width

frame.distr<-ctrl_pixeldistr.list[[3]]
frame.width<-length(frame.distr[frame.distr>0]) #number of Y rows with white pixels in them
frame.width<-frame.width/frame1.width
widths[3]<-frame.width


#plot widhts
par(lwd=2)
barcenters<-barplot(height=widths,beside=T,plot=T,width=0.5,
                    ylim=c(0,2),xaxt="n",las=1,
                    col=c("gray","gray40","black"),border="black",
                    ylab="",main="Normalised smoke width - dead")
axis(1,at=barcenters,labels=c("before","down-sweep","after"),lty=0,cex.axis=1)



dev.off()























#--------------------------------------------------------------------------------------------------------------------------------
# PID analysis
#--------------------------------------------------------------------------------------------------------------------------------
# creating the data for the analysis


# load data set with the synchronised PID recordings and antennal movements
Data<-read.csv("../Data/PIDdata_activesmelling.csv",sep=";",dec=".",header=T,check.names=F)




# Calculate combined angular speed

#need to remove NAs
H.angle<-Data$H.angle
H.angle[is.na(H.angle)]<-0
V.angle<-Data$V.angle
V.angle[is.na(V.angle)]<-0

Data$ResultingAngle<-sqrt((H.angle^2)+(V.angle^2)) 

datlist<-split(Data,f=Data$Test)
datlist<-lapply(datlist,
                FUN=function(df){ 
                  
                  df$ResultingAngularSpeed<-abs(c(NA,diff(df$ResultingAngle)))*100 #speed °/s
                  df$ResultingAngularSpeed<-repeat.before(df$ResultingAngularSpeed)
                  
                  return(df)
                  
                })
Data<-dplyr::bind_rows(datlist)





# moving average, then choose window around a local max
datlist<-split(Data,f=Data$Test)
datlist<-lapply(datlist,
                FUN=function(df){ 
                  
                  df$RollingMeanSpeed<-MovingFUN(df$ResultingAngularSpeed,method="centered",winsize=100,shift=1,
                                                 FUN="mean",new.length=nrow(df))
                  #smoothing further
                  df$RollingMeanSpeed.2<-MovingFUN(df$RollingMeanSpeed,method="centered",winsize=100,shift=1,
                                                   FUN="mean",new.length=nrow(df))
                  
                  return(df)
                  
                })
Data<-dplyr::bind_rows(datlist)



# plot
lapply(datlist,
       FUN=function(df){
         
         # pdf(file=paste(unique(df$Test),"rawspeed&PID.pdf",sep=""),width=10,height=7.5,bg="white",useDingbats=F)
         
         
         par(mar=c(5,4,4,4)+0.3)
         
         plot(df$Count*0.01,df$pidvoltage,col="forestgreen",type="l",lwd=2,
              ylim=c(0,max(df$pidvoltage)+0.1),ylab="PID [V]",xlab="Time [s]",
              xaxs="i",yaxs="i",
              main=unique(df$Test))
         
         par(new=TRUE)
         speed<-df$ResultingAngularSpeed*0.01
         plot(df$Count*0.01,speed,type="l",lwd=2,ylim=c(0,4),col="gray",
              axes=FALSE,bty="n",xlab="",ylab="",xaxs="i",yaxs="i",)
         axis(side=4,at=c(0,50,100,150,200,250)*0.01,labels=c(0,50,100,150,200,250))
         mtext("Angular speed [°/s]",side=4,line=3)
         
         
         # just overlaying the PID trace again for it to be on top of the speed
         par(new=TRUE)
         plot(df$Count*0.01,df$pidvoltage,col="forestgreen",type="l",lwd=2,
              ylim=c(0,max(df$pidvoltage)+0.1),
              axes=FALSE,bty="n",xlab="",ylab="",xaxs="i",yaxs="i",)
         
         
         # dev.off()
         
       })





# select max windows based on smoothed mean speed
winsize=200 

# pick the max windows of each trial
windowlist<-split(Data,f=Data$Test)
windowlist<-lapply(windowlist,
               FUN=function(df){
                 
                 a<-which(df$RollingMeanSpeed.2==max(df$RollingMeanSpeed.2,na.rm=T) )
                 df$maxspeed<-0
                 df$maxspeed[a]<-1
                 
                 #select a winsize+1 frame window around that frame
                 if(sum(df$maxspeed)!=0 & unique(df$mvt)=="free"){
                   
                   a<-which(df$maxspeed==1)
                   
                   if((a-winsize/2)<0){ #if the peak is too close to the starting edge to go to a-winsize/2
                     
                     df<-df[1:(winsize+1),] 
                     
                   } else if((a+winsize/2)>nrow(df)){ #if the peak is too close to the end edge to go to a+winsize/2
                     
                     df<-df[(nrow(df)-winsize):nrow(df),] 
                     
                   } else { # if the peak is far enough from the edges
                     
                     df<-df[(a-winsize/2):(a+winsize/2),] #select a winsize+1 frame window around that frame
                     
                   }
                   
                   return(df)
                   
                   
                 } else { 
                   
                   return(NULL)
                   
                 }
                 
               })
windowlist<-windowlist[-which(sapply(windowlist, is.null))] #remove empty slots




# Choose windows in ctrl data
ctrllist<-split(Data,f=Data$Test)
#first window in ctrl
ctrllist1<-lapply(ctrllist,
                  FUN=function(df){
                    
                    if(unique(df$mvt)=="fixed"){
                      
                      a<-200 # just take a point in the trial
                      df<-df[(a-winsize/2):(a+winsize/2),] #select a winsize1 frame window around that frame
                      
                      return(df)
                      
                    } else { 
                      
                      return(NULL)
                      
                    }
                    
                  })
ctrllist1<-ctrllist1[-which(sapply(ctrllist1, is.null))] #remove empty slots
names(ctrllist1)<-paste(names(ctrllist1),"1",sep="_")

#second window later in ctrl
ctrllist2<-lapply(ctrllist,
                  FUN=function(df){
                    
                    if(unique(df$mvt)=="fixed"){
                      
                      a<-600 # just take a point in the trial
                      df<-df[(a-winsize/2):(a+winsize/2),] #select a winsize+1 frame window around that frame
                      
                      return(df)
                      
                    } else { 
                      
                      return(NULL)
                      
                    }
                    
                  })
ctrllist2<-ctrllist2[-which(sapply(ctrllist2, is.null))] #remove empty slots
names(ctrllist2)<-paste(names(ctrllist2),"2",sep="_")

#merge
ctrllist<-c(ctrllist1,ctrllist2)



# MERGE ALL
fulllist<-c(windowlist,ctrllist)









#--------------------------------------------------------------------------------------------------------------------------------
# PID Peaks analysis
#--------------------------------------------------------------------------------------------------------------------------------


# NEEDS "fulllist" created above


### Peak rate on first derivative ###
Peaks.list<-lapply(fulllist,
                        FUN=function(df){ 
                          
                          print(unique(df$Test))
                          
                          #PID peaks
                          pidpeaks<-peaks(time=(1:length(df$pidvoltage)),var=df$pidvoltage,
                                          halfwin=7,thresh=0.03) #halfwin=5, thresh=0.03
                          #PID values at peaks
                          peakvalues<-df$pidvoltage[which(pidpeaks!=0)]
                          #only keep consecutively altering peaks (ignore consecutive positives or negatives)
                          pidpeaks.nonzero<-pidpeaks[pidpeaks!=0]
                          nondoubles<-which(diff(pidpeaks.nonzero)!=0)
                          peakvalues<-peakvalues[nondoubles]
                          
                          #Mean contrasts
                          peak.contrasts<-diff(peakvalues)
                          peak.contrasts<-peak.contrasts[peak.contrasts>=0] #positive contrasts only
                          meancontrast<-mean((peak.contrasts),na.rm=T)
                          if(is.nan(meancontrast)==T){ meancontrast<-0 } #no peaks, no contrast
                          
                          
                          #rate of PID +-peaks
                          peaks.nondoubles<-which(df$pidvoltage %in% peakvalues)
                          newpeaks<-rep(0,nrow(df))
                          newpeaks[peaks.nondoubles]<-1
                          newpeaksrate<-rate.calc(peakvect=newpeaks,method="instantaneous")
                          meannewpeaksrate<-mean(newpeaksrate,na.rm=T) #mean in the window
                          
                          
                          #avg angular speed in that window
                          avgspeed<-mean(df$RollingMeanSpeed,na.rm=T)
                          
                          
                          
                          #new df
                          metrics<-data.frame("window"=paste(unique(na.omit(df$Test)),min(df$Count),sep="_"),
                                              "meanpidpeakrate.all"=meannewpeaksrate,
                                              "meancontrast"=meancontrast,
                                              "meanantspeed"=avgspeed
                          )
                          
                          return(metrics)
                          
                          
                        }
)
Peaks.df<-dplyr::bind_rows(Peaks.list)





plot(Peaks.df$meanantspeed,Peaks.df$meancontrast)
plot(Peaks.df$meanantspeed,Peaks.df$meanpidpeakrate.all)










###########################
### Modelling Peak rate ###
###########################


ModelData<-Peaks.df

hist(ModelData$meanpidpeakrate.all)


# Fit
set.seed(13062022)
brmsmod<-brm(meanpidpeakrate.all ~ meanantspeed,
             data=ModelData, 
             family=gaussian(),
             iter=10000, chains=4, cores=4, backend="cmdstanr",
             control=list(adapt_delta=0.99,max_treedepth=20))

summary(brmsmod)


# pp check
pp<-pp_check(brmsmod,ndraws=200) 
pp

# convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce



# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))



#saving parameters for the plot
peakrate_estimates<-ce$meanantspeed$estimate
peakrate_speeds<-ce$meanantspeed$meanantspeed
peakrate_upper<-ce$meanantspeed$upper__
peakrate_lower<-ce$meanantspeed$lower__









##########################
### Modelling Contrast ###
##########################


ModelData<-Peaks.df

hist(ModelData$meancontrast)


# Fit
set.seed(13062022)
brmsmod<-brm(meancontrast ~ meanantspeed,
             data=ModelData, 
             family=gaussian(),
             iter=10000, chains=4, cores=4, backend="cmdstanr",
             control=list(adapt_delta=0.99,max_treedepth=20))

summary(brmsmod)



# pp check
pp<-pp_check(brmsmod,ndraws=200) 
pp

# convergence check
plot(brmsmod)


#quick effect plot
ce<-conditional_effects(brmsmod)
ce



# MCMC samples
simval<-as.data.frame(as_draws_df(brmsmod))



#saving parameters for the plot
contrast_estimates<-ce$meanantspeed$estimate
contrast_speeds<-ce$meanantspeed$meanantspeed
contrast_upper<-ce$meanantspeed$upper__
contrast_lower<-ce$meanantspeed$lower__











### effect plots

pdf(file=paste("./Plots/PIDContrast&Peakrate_vs_antspeed_effectplots.pdf",sep=""),
  width=10,height=6,bg="white",useDingbats=F)


par(mfrow=c(1,2))


#peak rate
plot(ModelData$meanantspeed,rep(NA,nrow(ModelData)),
     ylim=c(0,14), 
     yaxs="i",xaxs="i",
     ylab="Mean PID peak rate [Hz]")

#credible interval
polygon(x=c(peakrate_speeds,rev(peakrate_speeds)),
  y=c(peakrate_upper,rev(peakrate_lower)),
  col=rgb(gray40[1],gray40[2],gray40[3],alpha=100,names=NULL,maxColorValue=255),border=NA)
#points
points(ModelData$meanantspeed,ModelData$meanpidpeakrate.all,
       pch=21,cex=2,lwd=2,
       bg=rgb(gray60[1],gray60[2],gray60[3],alpha=0.6*255,maxColorValue=255),)
lines(peakrate_speeds,peakrate_estimates,lwd=2)





# contrast
plot(ModelData$meanantspeed,rep(NA,nrow(ModelData)),
     ylim=c(0,0.85), #meancontrast
     yaxs="i",xaxs="i",
     ylab="Mean pos contrasts [V]")

#credible interval
polygon(x=c(contrast_speeds,rev(contrast_speeds)),
        y=c(contrast_upper,rev(contrast_lower)),
        col=rgb(gray40[1],gray40[2],gray40[3],alpha=100,names=NULL,maxColorValue=255),border=NA)
#points
points(ModelData$meanantspeed,ModelData$meancontrast,
       pch=21,cex=2,lwd=2,
       bg=rgb(gray60[1],gray60[2],gray60[3],alpha=0.6*255,maxColorValue=255))
lines(contrast_speeds,contrast_estimates,lwd=2)



dev.off()









