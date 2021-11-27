###
# set input variables
###

source("DifferentiationPredictor_functions.R")

# name the input files
limTrainFile<- "GSE16997_series_entrez.txt"
the337File<- "UNC337arraydata_imputedCollapsed.txt"
mmsFile <- "14unsortedSamples_noReplicates_forPaper.txt"
mouseFile<- "Jason-122-Mouse.txt"

# name the output file for the figure
outputFile<- "Figure4.pdf"
# name the output file for the model
limModelFile<- "differentiationCentroids_LimDWD.txt"


####
# start script
####
# dwd training with Lim
###
x<-readarray(limTrainFile,hr=4)
x$xd<-standardize(medianCtr(x$xd))

# set pL as the origin
plmeds<-apply(x$xd[,x$classes$type=="pL"],1,median,na.rm=T)
x.plcent<-x$xd-plmeds

# relabel MaSC/pL as 1/2
mascVspL<-as.vector(x$classes$type)
mascVspL[mascVspL=="MaSC"]<-2
mascVspL[mascVspL=="pL"]<-1
mascVspL[!(mascVspL=="1" | mascVspL=="2")]<-NA
mascVspL.x<-x.plcent[,!is.na(mascVspL)]
mascVspL<-mascVspL[!is.na(mascVspL)]
# get the DWD loading
mascVspL.dwd<-dwd(mascVspL.x,mascVspL)

# relabel pL/mL as 1/2
mLVspL<-as.vector(x$classes$type)
mLVspL[mLVspL=="mL"]<-2
mLVspL[mLVspL=="pL"]<-1
mLVspL[!(mLVspL=="1" | mLVspL=="2")]<-NA
mLVspL.x<-x.plcent[,!is.na(mLVspL)]
mLVspL<-mLVspL[!is.na(mLVspL)]
# get the DWD loading
mLVspL.dwd<-dwd(mLVspL.x,mLVspL)

# format and write the two DWD vectors
trainedVectors<-cbind(as.vector(t(mascVspL.dwd)),as.vector(t(mLVspL.dwd)))
rownames(trainedVectors)<-rownames(x$xd)
colnames(trainedVectors)<-c("MS","mL")
write.table(trainedVectors,limModelFile,sep="\t",col.names=NA)

####
#  read in data and prepare for predictions
####

# assumed that data sets are annoted with common (Entrez) IDs
the337<-readarray(the337File,hr=5)
mms<-readarray(mmsFile,hr=2)
mouse<-readarray(mouseFile,hr=4)
mouse$xd<-medianCtr(mouse$xd)
y25<-readarray(limModelFile,hr=1)

# open the file for plotting
pdf(outputFile,width=10)
layout(matrix(c(1,2,3,4,4,4),nrow=2,byrow=T),widths=c(1,0.5,1,1,1,1),heights=c(1,1.3))

####
# prediction in MMS
####
diffScore<-assignDiffScore.dwd(mms,y25)

par(las=1,cex.axis=1.5,cex.lab=1.5,mar=c(5,7,2,0))
ylims<-range(-0.3,0.33)
stripchart(diffScore~rep("MMS",dim(mms$xd)[2]),main="A",pch="",xaxt="n",ylim=ylims,xlab="")
title(ylab="Differentiation Score",line=5)
axis(2,at=seq(-0.3,0.3,.15),labels=seq(-0.3,0.3,.15))
text(c(-0.17,-0.17,-0.17),c(.3,0,-.28),labels=c("Mature Luminal","Luminal Progenitor","MaSC"),cex=1.2)

segments(-0.17,0.05,-0.17,0.22,lwd=2)
segments(-0.16,0.19,-0.17,0.22,lwd=2)
segments(-0.18,0.19,-0.17,0.22,lwd=2)

segments(-0.17,-0.05,-0.17,-0.22,lwd=2)
segments(-0.16,-0.08,-0.17,-0.05,lwd=2)
segments(-0.18,-0.08,-0.17,-0.05,lwd=2)

par(las=1,cex.axis=1.5,cex.lab=1.5,mar=c(5,0,2,0))
boxplot(diffScore~rep("MMS",dim(mms$xd)[2]),outline=F,border=8,main="B. MMS",ylim=ylims,yaxt="n")
axis(1,at=1,label="MMS")
mmscols<-as.vector(mms$classes$subtype)
mmscols[mmscols=="CL"]<-"yellow"
mmscols[mmscols=="not"]<-"black"
xs<-jitter(rep(1,dim(mms$xd)[2]),5)
points(xs,diffScore,pch=3,cex=1.5,col=1)

###
# UNC 337
###
the337$classes$subtype<-as.vector(the337$classes$subtype)
the337$classes$subtype[the337$classes$subtype=="LumA"]<-"LA"
the337$classes$subtype[the337$classes$subtype=="LumB"]<-"LB"
the337$classes$subtype[the337$classes$subtype=="Her2"]<-"H2"
the337$classes$subtype[the337$classes$subtype=="Basal"]<-"BL"
the337$classes$subtype[the337$classes$subtype=="Claudin"]<-"CL"
the337$classes$subtype[the337$classes$subtype=="Normal Breast"]<-"NBL"
the337$classes$subtype<-factor(the337$classes$subtype,levels=c("LA","LB","H2","BL","CL","NBL"))
diffScore<-assignDiffScore.dwd(the337,y25)

par(las=1,cex.axis=1.5,cex.lab=1.5,mar=c(5,0,2,2))
boxplot(diffScore~the337$classes$subtype,main="C. UNC 337",col=c("blue","lightblue","hotpink","red","yellow","green"),
												ylim=ylims,ylab="Differentiation Score",outline=F,yaxt="n")
stripchart(diffScore~the337$classes$subtype,method="jitter",vertical=T,add=T,pch=3,cex=.5)
segments(4,0.1,4,0.15,lwd=2)
segments(5,0.1,5,0.15,lwd=2)
segments(4,0.15,5,0.15,lwd=2)
text(4.5,0.175,labels=c("*"),cex=2.5)
t.test(diffScore[the337$classes$subtype=="BL"],diffScore[the337$classes$subtype=="CL"])

####
# prediction in JH mouse
####
diffScore<-assignDiffScore.dwd(mouse,y25)
mouse$classes[,3]<-factor(mouse$classes[,3],levels=c("MMTV-Neu/PyMT","WAPT121/WAPTag","WAP-Myc","WAP-Int3","p53null/p53het IR","TgC3(1)-Tag","Normal Breast","DMBA/Wnt1","BRCA1/p53/Wnt1","Claudin-low"))

par(las=2,cex.axis=1.5,cex.lab=1.5,mar=c(13,7,2,2))
boxplot(diffScore~mouse$classes[,3],ylab="",main="D. Mouse tumors",outline=F,border=8,col=c(rep("white",9),"yellow"),yaxt="n")
stripchart(diffScore~mouse$classes[,3],add=T,pch=3,vertical=T,method="jitter",cex=1.25)
title(ylab="Differentiation Score",line=5)
axis(2,at=seq(-0.2,0.2,.05),labels=seq(-0.2,0.2,.05))

dev.off()
