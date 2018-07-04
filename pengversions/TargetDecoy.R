###################
#this script is used for Target-Decoy Search and untargeted chemical analysis
###Hui, 20170921

####################for drinking water, cpp funtion, the oxygen number could be 1.2*carbon+3, cutoff e5, S/N>5, change is.BrCl
library(xcms)
library(MassSpecWavelet)
library(Rcpp)
library(RcppArmadillo)
library(isopat)
data(iso_list)
library(ggplot2)
setwd("C:/procedure/R work/Target_Decoy")
source("TargetDecoyFun.r")

##########identify lockmass#############
setwd("C:/procedure/R work/Target_Decoy/data")
msfiles<-list.files()
mzwin<-5###2.5ppm for mz cutoff
timewin<-0.5###30 sec for rt cutoff
timewin2<-2####60 sec for rt cutoff, since library was established for a long time
xset<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=10,nSlaves=1,polarity="negative")##peak width, the min and max range of chromatographic peaks in seconds
result<-findlock(xset,2000,0.002)##xset, intensity threshold, mzstep
setwd("C:/procedure/R work/Target_Decoy")
write.table(result, file="lockmass.csv", sep = ',',row.names=FALSE,col.names=c("mz","minintensity","sampleID"))


#################plot LockMass######################
setwd("C:/procedure/R work/Target_Decoy/CalData")
msfiles<-list.files()
testMass<-151.0395*(1-0*10^(-6))
polarity<--1##if neg -1, if pos 1
LockMass.POS<-c(81.070425,93.070425,105.070425,139.11229,151.042199,171.138505,413.26696,445.120583)##Lock Mass for positive
xrawdata<-xcmsRaw(msfiles[1],includeMSn=TRUE)
setwd("C:/procedure/R work/Target_Decoy")
LockMass.NEG<-read.table("lockmass.csv",header=TRUE,sep=',')
LockMass.NEG<-LockMass.NEG$Lock
ppmshift<-plotlock(xrawdata,testMass,2)
ppm.matrix<-matrix(rep(0,3*nrow(ppmshift)*ncol(ppmshift)),ncol=3,nrow=nrow(ppmshift)*ncol(ppmshift))
for (i in 1:nrow(ppmshift)){
     for (j in 1:ncol(ppmshift)){
     ppm.matrix[ncol(ppmshift)*(i-1)+j,1]<-LockMass.POS[i]##m/z of lockmass
     ppm.matrix[ncol(ppmshift)*(i-1)+j,2]<-ppmshift[i,j]##ppm shift
     ppm.matrix[ncol(ppmshift)*(i-1)+j,3]<-i##groupid}
}}
colnames(ppm.matrix)<-c('mz','ppm','group')
ppm.matrix<-data.frame(ppm.matrix)
shift3<-ppm.matrix$ppm
ppm.matrix$scan<-rep(1:ncol(ppmshift),nrow(ppmshift))
ppm.matrix$scan<-ppm.matrix$scan*30/length(ppm.matrix$scan)
dataplot<-ggplot(ppm.matrix,aes(x=scan,y=ppm,color=factor(group)))+
geom_point(size=3)
dataplot+theme(axis.text.x = element_text(size=10,color='black'),
          axis.title.x=element_text(face="bold",size=12,color='black'),
          axis.text.y = element_text(size=10,color='black'),
          axis.title.y=element_text(face="bold",size=12,color='black'),
          plot.title=element_text(size=14,face='bold'))+
          scale_y_continuous(name='ppm',limits=c(-10,5))+
          xlab('scan')+
          ylab('ppm')
          
test<-ppm.matrix[3:500,]
model<-loess(ppm~scan,test,span=0.05)
plot(test$scan,test$ppm,size=1)
lines(test$scan,model$fit,col='red')          


####################Mass Calibration###############
#setwd(file.path(maindir, "housedust", "rawdata"))
setwd("C:/Users/skutarna/Dropbox/MSc/Mass Spec/Target Decoy/housedust/rawdata")
msfiles<-list.files()
off.set<--10###bigshift to narrow down
for (i in 1:length(msfiles)){
    #setwd(file.path(maindir, "housedust", "rawdata"))
    setwd("C:/Users/skutarna/Dropbox/MSc/Mass Spec/Target Decoy/housedust/rawdata")
    xrawdata<-xcmsRaw(msfiles[i],includeMSn=TRUE)
    xrawdata@env$mz<-xrawdata@env$mz*(1-off.set*10^(-6))
    if (polarity==1){
    xrawdata<-MassCal(xrawdata,LockMass.POS,'POS',10)
    xrawdata<-MassCal(xrawdata,LockMass.POS,'POS',5)}
    if (polarity==-1){
    xrawdata<-MassCal(xrawdata,LockMass.NEG,'NEG',10)###two step calibration
    xrawdata<-MassCal(xrawdata,LockMass.NEG,'NEG',5)}
    #setwd(file.path(maindir, "housedust", "caldata"))
    setwd("C:/Users/skutarna/Dropbox/MSc/Mass Spec/Target Decoy/housedust/caldata")
    write.mzdata(xrawdata,msfiles[i])####save the calibrated data to new files
}

###################Figure#############
test<-data.frame(mz=table[,1],ppm=table[,2]*10^6)                   
ggplot(test,aes(mz,ppm))+geom_point(color='red')+geom_smooth()+
theme(axis.text.x = element_text(size=10,color='black'),
panel.background = element_rect(fill=NA), 
axis.line.x=element_line(colour="black"),
axis.line.y=element_line(colour="black"),
          axis.title.x=element_text(face="bold",size=12,color='black'),
          axis.text.y = element_text(size=10,color='black'),
          axis.title.y=element_text(face="bold",size=12,color='black'),
          plot.title=element_text(size=14,face='bold'))+
          xlab('m/z')+
          ylab('ppm')  

######smoothing the data to average the m/z values##########
setwd("C:/procedure/R work/Target_Decoy/CalData")
msfiles<-list.files()
mzwin<-5###2.5ppm for mz cutoff
timewin<-0.5###30 sec for rt cutoff
timewin2<-2####60 sec for rt cutoff, since library was established for a long time
xset<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,20),snthresh=10,nSlaves=1)##peak width, the min and max range of chromatographic peaks in seconds
newpeaks<-xset@peaks##########smooth the m/z values
minmz<-min(newpeaks[,1])
maxmz<-max(newpeaks[,1])
minrt<-min(newpeaks[,4])
maxrt<-max(newpeaks[,4])

#################function to smooth the three scanning points to get the m/z######
##########m/z=mz1*W1+mz2*W2+mz3*W3, where Wi=1*Intensityi/sum(intensity)
for (i in 1:length(msfiles)){
     index<-which(newpeaks[,11]==i)
     xraw<-xcmsRaw(msfiles[i],includeMSn=TRUE)
for (j in 1:length(index)){
     mz.value<-newpeaks[index[j],1]
     mzmin<-max(minmz,mz.value-mz.value*2*10^(-6))
     mzmax<-min(maxmz,mz.value+mz.value*2*10^(-6))
     rt.value<-newpeaks[index[j],4]
     rtmin<-max(minrt,rt.value-30)
     rtmax<-min(maxrt,rt.value+30)
     peak<-rawEIC(xraw,mzrange=cbind(mzmin,mzmax),rtrange=cbind(rtmin,rtmax))
     if (max(peak$intensity)==0){next}
     index.peak<-which.max(peak$intensity)
     indexpeak<-peak$scan[index.peak[1]]
     index.min<-max(1,indexpeak-1)
     index.max<-min(max(peak$scan),indexpeak+1)
     scanNum.min<-c(xraw@scanindex[index.min],xraw@scanindex[index.min+1])
     id.min<-(scanNum.min[1]+1):scanNum.min[2]
     idmin<-which(abs(xraw@env$mz[id.min]-mz.value)<2*10^(-6)*mz.value)
     if (length(idmin)<1){mz.min<-mz.value}else{
     idmin<-idmin[1]
     mz.min<-xraw@env$mz[id.min[idmin]]}
     scanNum.max<-c(xraw@scanindex[index.max],xraw@scanindex[index.max+1])
     id.max<-(scanNum.max[1]+1):scanNum.max[2]
     idmax<-which(abs(xraw@env$mz[id.max]-mz.value)<2*10^(-6)*mz.value)
     if (length(idmax)<1){mz.max<-mz.value}else{
     idmax<-idmax[1]
     mz.max<-xraw@env$mz[id.max[idmax]]}
     mz.cal<-(mz.value*peak$intensity[index.peak]+mz.max*peak$intensity[max(1,index.peak-1)]+mz.min*peak$intensity[min(length(peak$scan),index.peak+1)])/(peak$intensity[max(1,index.peak-1)]+peak$intensity[index.peak]+peak$intensity[min(length(peak$scan),index.peak+1)])
     newpeaks[index[j],1]<-mz.cal
}}
setwd("C:/procedure/R work/Target_Decoy")
write.table(newpeaks, file="peaks.csv", sep = ',',row.names=FALSE)





############read the library############
setwd("C:/procedure/R work/Target_Decoy/library")
Library<-read.table("library.csv",header=TRUE,sep=',')
TargetDatabase<-TargetConstruct(Library)
ImAdducts<-c(4.002603,7.016005,9.012183,11.009305,19.992439,26.981541,27.976928,39.962383,39.962591,44.955914,47.947947)##Only those ions with lesser than 50 m/z,implausible adducts
DecoyDatabase<-DecoyConstruct(TargetDatabase,ImAdducts)

############peaks detection##############################
setwd("C:/procedure/R work/Target_Decoy/CalData")
xset.all<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,10),snthresh=10,nSlaves=min(6,length(msfiles)))
xset_group<-group(xset.all,bw=30,minsamp=1,minfrac=0.1,mzwid=0.001)
len<-length(xset_group@groupidx)#group number
len2<-length(msfiles)#data files
Peak.init<-array(rep(0,len*(len2+3)),dim=c(len,(len2+3)))##columns are m/z, rt,sampleID 
for (i in 1:len){
    temp<-unlist(xset_group@groupidx[i])
    len3<-length(temp)
    for (j in 1:len3){
         index1<-xset_group@peaks[temp[j],11]
         Peak.init[i,1]<-xset_group@peaks[temp[j],1]##mz
         Peak.init[i,2]<-xset_group@peaks[temp[j],4]/60##rt
         Peak.init[i,3]<-xset_group@peaks[temp[j],11]###window id
         Peak.init[i,index1+3]<-max(xset_group@peaks[temp[j],9],Peak.init[i,index1+3])##intensity,select the maximum one if there are multiple peaks for one sample
         }}

for (i in 1:nrow(Peak.init)){####replace the sample id as the one with maximal abundance
    sample.id<-Peak.init[i,4:ncol(Peak.init)]
    Peak.init[i,3]<-which.max(sample.id)
}         

######################retention time alignment#################
setwd("C:/procedure/R work/Target_Decoy/library")
Reference<-read.table("Reference.csv",header=TRUE,sep=',')
ReferDatabase<-TargetConstruct(Reference)
setwd("C:/procedure/R work/Target_Decoy/CalData")
refermatch<-DatabaseMatch(Peak.init,ReferDatabase,2*10^(-6),msfiles)
index.refer<-which(refermatch[,10]>0)
refermatch<-refermatch[index.refer,]
for (i in 1:nrow(refermatch)){
     refermatch[i,3]<-ReferDatabase[refermatch[i,6],1] 
}

for (i in 1:nrow(refermatch)){
    reg<-lm(refermatch[,2]~refermatch[,3])
    coeff<-reg$coefficients
    res<-refermatch[,3]*coeff[2]+coeff[1]-refermatch[,2]
    index<-which.max(abs(res))##find the maximal residual, and delete it
    refermatch<-refermatch[-index,]
    if (coeff[2]>0.99||nrow(refermatch)<2){break}}
plot(refermatch[,3]*coeff[2]+coeff[1],refermatch[,2])
TargetDatabase[,1]<-TargetDatabase[,1]*coeff[2]+coeff[1]###correct the retention time of the database
DecoyDatabase[,1]<-DecoyDatabase[,1]*coeff[2]+coeff[1]    


############Database Match#########################
Target.score<-DatabaseMatch(Peak.init,TargetDatabase,2*10^(-6),msfiles)
Decoy.score<-DatabaseMatch(Peak.init,DecoyDatabase,2*10^(-6),msfiles)
index<-which(Target.score[,10]>1.5)
index1<-which(Decoy.score[,10]>1.5)
plot(Target.score[index,10])
points(Decoy.score[index1,10],col="red")

##########FDR calculation######################
score.save<-NULL
for (k in 1:100){##repeat the decoy database for 100 times, for random data
    print(paste('iteration...',k,sep=''))
    DecoyDatabase<-DecoyConstruct(TargetDatabase,ImAdducts)
    Decoy.score<-DatabaseMatch(Peak.init,DecoyDatabase,2*10^(-6),msfiles)
    for (j in 1:99){
        score.fdr<-0.03*j
        index<-which(Target.score[,10]>score.fdr)
        index1<-which(Decoy.score[,10]>score.fdr)
        if (length(index1)/length(index)<0.05){
        score.save<-c(score.save,score.fdr)
        break}
    }}
index<-which(Target.score[,10]>mean(score.save))
FinalID<-Target.score[index,]

########delete the duplicate information, the isotopic peaks##############
score.save<-NULL
ratio<-NULL
for (j in 1:99){
    score.fdr<-0.03*j
    score.save<-c(score.save,score.fdr)
    index<-which(Target.score[,10]>score.fdr)
    index1<-which(Decoy.score[,10]>score.fdr)
    ratio<-c(ratio,1-length(index1)/length(index))
}
plot(score.save,ratio)