# require(AnalyzeFMRI)

### Create new output folder (add _1, _2 etc if output folder exists) ###
newOutputFolder <- function(main,folder) {
  origwd <- getwd()
  setwd(paste(origwd,main,sep="/"))
  outputFolder <- file.path(getwd(), folder)
  checkFolder <- dir(pattern = paste(folder,"_[0-9]",sep=""))
  lastNum <- 0 
  if(length(checkFolder)>0)  {
    lastFolder <- checkFolder[length(checkFolder)]
    lastNum <- as.numeric(unlist(strsplit(lastFolder,"[_]"))[length(unlist(strsplit(lastFolder,"[_]")))])
  }
  outputFolder <- paste(outputFolder,"_", lastNum+1, sep="")
  dir.create(outputFolder, showWarnings = FALSE)
  setwd(origwd)
  cat(paste("Created output folder:",outputFolder,"\n\n"))
  return(outputFolder)
}

### Print start/finish time ###
printTime <- function(type) {
  if(type=="start") {
    startTime <<- Sys.time()
    cat(paste("\n### Started at ",format(startTime, "%a %b %d %X %Y")," ###\n\n",sep=""))
  }
  if(type=="end") {
    td <- as.numeric(floor(difftime(Sys.time(), startTime, units=c("secs"))))
    if(td/3600>1){th=floor(td/3600);td=td-(th*3600)}else{th=0};if(td/60>1){tm=floor(td/60);td=td-(tm*60)}else{tm=0}
    cat(paste("\n\n### Finished at ",format(Sys.time(), "%a %b %d %X %Y")," - Runtime: ",th,"h",tm,"m",td,"s ###\n",sep=""))
  }
}

#print progress function
lastPrint <- 10
printProgress <- function(current,total) {
  percent <- round((current/total)*10)
  if(percent!=lastPrint&percent>0) { cat(percent); lastPrint <<- percent }
}

### Analytic method of thresholding weights (Gaonkar & Davatzikos, 2013) ###
thresholdWeights<- function(DATA,w,pthresh) {
  threshold_z <- qnorm(1-pthresh/2) #p value >> z-score
  K <- DATA%*%t(DATA) 
  J <- matrix(1,nrow(DATA))
  KK <- ginv(K)
  Z <- KK+(KK%*%J%*%ginv(-t(J)%*%KK%*%J)%*%t(J)%*%KK)
  C <- t(DATA)%*%Z
  SD <- apply(C,1,function(x)sqrt(sum(x*x))) #get SD for each row
  SDt <- SD*threshold_z #threshold SD
  wt <- ifelse(abs(w)>SDt,w,0) #threshold weight vector
  return(wt)
}

### Script to write weights to a nifti file ###
writeMap <- function(w,output) {
  MAP <- MASK
  MAP[which(MASK!=0,arr.ind=TRUE)] <- w
  #f.write.analyze(SVM_MAP,output,size="float",pixdim=c(2,2,2)) #write as analyze - not good!
  header <- f.complete.hdr.nifti.list.create(file="map.nii", dim=c(3,91,109,91,1,1,1,1), datatype=16,
                                             pixdim=c(-1,2,2,2,0,0,0,0), descrip="Created with R",
                                             qform.code=2, sform.code=2, xyzt.units="\n",
                                             qoffset.x=90, qoffset.y=-126, qoffset.z=-72,
                                             srow.x=c(-2,0,0,90), srow.y=c(0,2,0,-126), srow.z=c(0,0,2,-72))
  f.write.nifti(MAP,file=output,size="float",L=header,nii=T) #write as nifti file
}

### Script to print weights to an image for visualization ###
printMap <- function(views,sl,transp,w,mri,output) {
  require(raster)
  
  #create map
  map <- MASK
  map[which(MASK!=0,arr.ind=TRUE)] <- w
  
  #define image colors
  palmap <- colorRampPalette(c("cyan", "blue", "black", "red", "yellow"), bias=1) #create ramp palette
  colmap <- palmap(21) #sample 21 colors from palette
  colmap <- paste(colmap,substring(rgb(0,0,transp,maxColorValue=100),6,7),sep="") #define alpha as HEX (ex: #FF00FF80)
  colmap[11] <- "#00000000" #make central color transparent
  palmri <- colorRampPalette(c("black", "white"), bias=1)
  colmri <- palmri(255)
  
  #define volume margins (based on map)
  if(any(c(dim(map)[1]!=91,dim(map)[2]!=109,dim(map)[3]!=91))){
    stop("Image dimensions do not match predefined values in script. Please make the necessary changes to the script!")}
  marx <- c(11,82); mary <- c(11,100); marz <- c(4,77)
  
  #setup output device
  png(paste(output,".png",sep=""),height=91*2*length(views)*sl[1],width=109*2*sl[2])
  tiles=c(length(views)*sl[1],sl[2]) #calculate ammount of tiles in image
  par(mfrow=tiles,mar=c(1,1,1,1),oma=c(1,1,1,1))
  
  #define ploting function
  plotmap <- function(view,slices,xmx,ymx) {
    for(s in slices) {
      if(view=="x") { m<-mri[s,,,]; z<-map[s,,,] #get slice
                      m<-t(m[,ncol(m):1]); z<-t(z[,ncol(z):1]) } #flip image
      if(view=="y") { m<-mri[,s,,]; z<-map[,s,,] #get slice
                      m<-t(m[,ncol(m):1]); z<-t(z[,ncol(z):1]) } #flip image
      if(view=="z") { m<-mri[,,s,]; z<-map[,,s,] } #get slice
      z[1,1]<-min(map); z[1,2]<-max(map) #center color pallette
      z<-ifelse(z==0,NaN,z) #make voxels with value=0 -> NaN so they are transparent
      rz<-raster(z,xmn=1,xmx=xmx,ymn=1,ymx=ymx); rm<-raster(m,xmn=1,xmx=xmx,ymn=1,ymx=ymx)
      plot(rm,col=colmri,axes=F,legend=F,box=F); plot(rz,col=colmap,axes=F,add=T,legend=F,box=F)
    }
  }
  
  #calculate slices
  ns<-sl[1]*sl[2]-1
  if(any(views=="x")) {
    int<-floor((marx[2]-marx[1]-1)/ns) #calculate slice interval
    mar<-floor(((marx[2]-marx[1])-(ns*int))/2) #calculate margin on each side
    slices<-seq(marx[1]+mar,marx[2]-mar-1,int) #generate slice sequence
    plotmap("x",slices,109,91)
  }
  if(any(views=="y")) {
    int<-floor((mary[2]-mary[1]-1)/ns)
    mar<-floor(((mary[2]-mary[1])-(ns*int))/2)
    slices<-seq(mary[1]+mar,mary[2]-mar,int)
    plotmap("y",slices,91,91)
  }
  if(any(views=="z")) {
    int<-floor((marz[2]-marz[1]-1)/ns)
    mar<-floor(((marz[2]-marz[1])-(ns*int))/2)
    slices<-seq(marz[1]+mar,marz[2]-mar,int)
    plotmap("z",slices,109,91)
  }
  
  dev.off()
  
  #print legend
  mat <- matrix(0,2,2)
  mat[1,1]<-min(map)
  mat[2,2]<-max(map)
  png(paste(output,"_legend.png",sep=""),height=400,width=500)
  par(mfrow=c(1,1),mar=c(2,2,2,2),oma=c(1,1,1,1))
  rmat <- raster(mat,xmn=1,xmx=200,ymn=1,ymx=200)
  plot(rmat,col=colmap,axes=F,legend=T,box=F)
  dev.off()
}