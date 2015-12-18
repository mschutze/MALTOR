#' MALTOR - Function to Import Data for Analysis
#' Manuel Schutze - August 2014 - manuels@ufmg.br
#' 
#' Required libraries: oro.nifti, raster, MASS

### Create new results folder (add _1, _2 etc if output folder exists) ###
newResultsFolder <- function(mainFolder,Script,Description) {
  origwd <- getwd()
  setwd(mainFolder)
  dir.create("results", showWarnings = FALSE)
  setwd("results")
  folderName <- paste(Script,Description,sep="_")
  fullOutputFolder <- file.path(getwd(), folderName)
  checkFolder <- dir(pattern = paste(folderName,"_[0-9]",sep=""))
  lastNum <- 0 
  if(length(checkFolder)>0)  {
    lastFolder <- checkFolder[length(checkFolder)]
    lastNum <- as.numeric(unlist(strsplit(lastFolder,"[_]"))[length(unlist(strsplit(lastFolder,"[_]")))])
  }
  fullOutputFolder <- paste(fullOutputFolder,"_", lastNum+1, sep="")
  dir.create(fullOutputFolder, showWarnings = FALSE, recursive = TRUE)
  setwd(origwd)
  cat(paste("\nCreated output folder:",fullOutputFolder,"\n"))
  return(fullOutputFolder)
}

### Print start/finish time ###
printTime <- function(type) {
  if(type=="start") {
    startTime <<- Sys.time()
    cat(paste("\n# Started at ",format(startTime, "%a %b %d %X %Y")," #\n",sep=""))
  }
  if(type=="end") {
    td <- as.numeric(floor(difftime(Sys.time(), startTime, units=c("secs"))))
    if(td/3600>1){th=floor(td/3600);td=td-(th*3600)}else{th=0};if(td/60>1){tm=floor(td/60);td=td-(tm*60)}else{tm=0}
    cat(paste("\n# Finished at ",format(Sys.time(), "%a %b %d %X %Y")," - Runtime: ",th,"h",tm,"m",td,"s #\n",sep=""))
  }
}

### Print histogram for imported images ###
printHist <- function(data,label) {
  temp.x <- NULL
  temp.y <- NULL
  for(i in 1:nrow(data)) {
    temp.x = c(temp.x, density(data[i,])$x)
    temp.y = c(temp.y, density(data[i,])$y)
  }
  xr <- range(temp.x)
  yr <- range(temp.y)
  plot(density(data[1,]), xlim = xr, ylim = yr, main = "Histogram for imported images")
  for(i in 1:nrow(data)) {
    lines(density(data[i,]), xlim = xr, ylim = yr, col = c(label)[i])
  }
  legend("topright", levels(label), fill=1+(0:nlevels(label)))
}

### Analytic method of thresholding weights (Gaonkar & Davatzikos, 2013) ###
thresholdWeights<- function(DATA,WEIGHTS,PTHRESH) {
  require("MASS")
  if((class(DATA)!="matrix")||(class(WEIGHTS)!="numeric")) stop("DATA has to be a matrix and w has to be a vector")
  threshold_z <- qnorm(1-PTHRESH/2) #p value >> z-score
  K <- DATA%*%t(DATA) 
  J <- matrix(1,nrow(DATA))
  KK <- ginv(K)
  Z <- KK+(KK%*%J%*%ginv(-t(J)%*%KK%*%J)%*%t(J)%*%KK)
  C <- t(DATA)%*%Z
  SD <- apply(C,1,function(x)sqrt(sum(x*x))) #get SD for each row
  SDt <- SD*threshold_z #threshold SD
  wt <- ifelse(abs(WEIGHTS)>SDt,WEIGHTS,0) #threshold weight vector
  return(wt)
}

### Script to write weights to a nifti file ###
writeMap <- function(w,output) {
  MAP <- maskNii
  MAP@.Data[maskIndex] <- w
  MAP@datatype <- 16
  MAP@bitpix <- 32
  MAP@cal_max <- max(w)
  MAP@cal_min <- min(w)
  MAP@scl_slope <- 1
  MAP@descrip <- "Created with MALTOR"
  writeNIfTI(MAP,output,gzipped=FALSE,verbose=FALSE) #write as nifti file
}

### Read Freesurfer binary MGH files (curvature, cortical thickness, etc) ###
readmgh <- function(path) {
  fid <- file(path, "rb")
  v <- readBin(fid, "integer", n=1, endian="big")
  dim <- readBin(fid, "integer", n=4, endian="big")
  dims <- dim[dim>1]
  skip <- readBin(fid, "integer", n=264, size=1, endian="big")
  data <- readBin(fid, "numeric", n=dims, size=4, endian="big")
  close(fid)
  return(data)
}

### Script to print weights to an image for visualization ###
# This version will accept images with sizes different from 91x109x91
# Works but needs improvement (align PET and MRI)
printMap <- function(views,sl,transp,w,mri,output) {
  require(raster)
  
  #clean up mri
  if(any(is.na(mri))) {  #remove NaNs
    mri <- apply(mri, MARGIN=c(1:3), FUN=function(x) ifelse(is.nan(x),0,x)) 
  }
  
  #create map
  map <- maskNii@.Data
  if(any(is.na(map))) {  #remove NaNs
    map <- apply(map, MARGIN=c(1:3), FUN=function(x) ifelse(is.nan(x),0,x)) 
  }
  if(max(map)==255) { #make voxels inside brain mask 0 or 1
    map <- map/255
  }
  map[maskIndex] <- w
  mval <- ifelse(max(w)>abs(min(w)),max(w),abs(min(w))) #gets absolute highest value (max or min)
  
  #define image colors
  palmap <- colorRampPalette(c("cyan", "blue", "black", "red", "yellow"), bias=1) #create ramp palette
  colmap <- palmap(21) #sample 21 colors from palette
  colmap <- paste(colmap,substring(rgb(0,0,transp,maxColorValue=100),6,7),sep="") #define alpha as HEX (ex: #FF00FF80)
  colmap[11] <- "#00000000" #make central color transparent
  palmri <- colorRampPalette(c("black", "white"), bias=1)
  colmri <- palmri(255)
  
  #define volume margins (based on map)
  dimx <- dim(map)[1]; dimy <- dim(map)[2]; dimz <- dim(map)[3]
  #marx <- c(11,82); mary <- c(11,100); marz <- c(4,77) #margins for 91x109x91
  #marx = 6,41 / mary = 6,50 / marz = 2,39
  marx <- c(round(0.121*dimx), round(0.901*dimx))
  mary <- c(round(0.101*dimy), round(0.917*dimy))
  marz <- c(round(0.217*dimz), round(0.901*dimz))
  
  #setup output device
  png(paste(output,".png",sep=""), height=dimx*4*length(views)*sl[1], width=dimy*4*sl[2])
  tiles <- c(length(views)*sl[1], sl[2]) #calculate ammount of tiles in image
  par(mfrow=tiles, mar=c(1,1,1,1), oma=c(1,1,1,1))
  
  #define ploting function
  plotmap <- function(view,slices,xmx,ymx) {
    for(s in slices) { #slices is a vector containing the coordinates to slice the image
      if(view=="x") { m<-mri[s,,]; z<-map[s,,] #get slice (m=mri / z=weight map)
                      m<-t(m[,ncol(m):1]); z<-t(z[,ncol(z):1]) } #flip image
      if(view=="y") { m<-mri[c(dimx:1),s,]; z<-map[c(dimx:1),s,] #get slice
                      m<-t(m[,ncol(m):1]); z<-t(z[,ncol(z):1]) } #flip image
      if(view=="z") { m<-mri[c(dimx:1),,s]; z<-map[c(dimx:1),,s] } #get slice
      z[1,1] <- -mval; z[1,2] <- mval #center color pallette
      z <- ifelse(z==0, NaN, z) #make voxels with value=0 -> NaN so they are transparent
      rz <- raster(z, xmn=1, xmx=xmx, ymn=1, ymx=ymx)
      rm <- raster(m, xmn=1, xmx=xmx, ymn=1, ymx=ymx)
      plot(rm, col=colmri, axes=FALSE, legend=FALSE, box=FALSE)
      plot(rz, col=colmap, axes=FALSE, add=TRUE, legend=FALSE, box=FALSE)
    }
  }
  
  #calculate slices
  ns <- sl[1]*sl[2]-1
  if(any(views=="x")) {
    int <- floor((marx[2]-marx[1]-1)/ns) #calculate slice interval
    mar <- floor(((marx[2]-marx[1])-(ns*int))/2) #calculate margin on each side
    slices <- seq(marx[2]-mar-1,marx[1]+mar,-int) #generate slice sequence (inverted for X: L>R)
    plotmap("x", slices, dimy, dimz)
  }
  if(any(views=="y")) {
    int <- floor((mary[2]-mary[1]-1)/ns)
    mar <- floor(((mary[2]-mary[1])-(ns*int))/2)
    slices <- seq(mary[1]+mar,mary[2]-mar,int)
    plotmap("y", slices, dimx, dimz)
  }
  if(any(views=="z")) {
    int <- floor((marz[2]-marz[1]-1)/ns)
    mar <- floor(((marz[2]-marz[1])-(ns*int))/2)
    slices <- seq(marz[1]+mar, marz[2]-mar, int)
    plotmap("z", slices, dimy, dimx)
  }
  
  dev.off()
  
  #print legend
  mat <- matrix(0,2,2)
  mat[1,1] <- -mval
  mat[2,2] <- mval
  png(paste(output,"_legend.png",sep=""), height=1000, width=1000)
  par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(1,1,1,1))
  rmat <- raster(mat, xmn=1, xmx=600, ymn=1, ymx=600)
  plot(rmat, col=colmap, axes=FALSE, legend=TRUE, box=FALSE)
  dev.off()
}

### Old version from printMap - works ONLY with image dimensions 91x109x91 ###
printMapOLD <- function(views,sl,transp,w,mri,output) {
  require(raster)
  
  #create map
  map <- maskNii@.Data
  map[maskIndex] <- w
  mval <- ifelse(max(w)>abs(min(w)),max(w),abs(min(w))) #gets absolute highest value (max or min)
  
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
  png(paste(output,".png",sep=""), height=91*4*length(views)*sl[1], width=109*4*sl[2])
  tiles <- c(length(views)*sl[1], sl[2]) #calculate ammount of tiles in image
  par(mfrow=tiles, mar=c(1,1,1,1), oma=c(1,1,1,1))
  
  #define ploting function
  plotmap <- function(view,slices,xmx,ymx) {
    for(s in slices) { #slices is a vector containing the coordinates to slice the image
      if(view=="x") { m<-mri[s,,]; z<-map[s,,] #get slice
                      m<-t(m[,ncol(m):1]); z<-t(z[,ncol(z):1]) } #flip image
      if(view=="y") { m<-mri[c(91:1),s,]; z<-map[c(91:1),s,] #get slice
                      m<-t(m[,ncol(m):1]); z<-t(z[,ncol(z):1]) } #flip image
      if(view=="z") { m<-mri[c(91:1),,s]; z<-map[c(91:1),,s] } #get slice
      z[1,1] <- -mval; z[1,2] <- mval #center color pallette
      z <- ifelse(z==0, NaN, z) #make voxels with value=0 -> NaN so they are transparent
      rz <- raster(z, xmn=1, xmx=xmx, ymn=1, ymx=ymx)
      rm <- raster(m, xmn=1, xmx=xmx, ymn=1, ymx=ymx)
      plot(rm, col=colmri, axes=FALSE, legend=FALSE, box=FALSE)
      plot(rz, col=colmap, axes=FALSE, add=TRUE, legend=FALSE, box=FALSE)
    }
  }
  
  #calculate slices
  ns <- sl[1]*sl[2]-1
  if(any(views=="x")) {
    int <- floor((marx[2]-marx[1]-1)/ns) #calculate slice interval
    mar <- floor(((marx[2]-marx[1])-(ns*int))/2) #calculate margin on each side
    slices <- seq(marx[2]-mar-1,marx[1]+mar,-int) #generate slice sequence (inverted for X: L>R)
    plotmap("x", slices, 109, 91)
  }
  if(any(views=="y")) {
    int <- floor((mary[2]-mary[1]-1)/ns)
    mar <- floor(((mary[2]-mary[1])-(ns*int))/2)
    slices <- seq(mary[1]+mar,mary[2]-mar,int)
    plotmap("y", slices, 91, 91)
  }
  if(any(views=="z")) {
    int <- floor((marz[2]-marz[1]-1)/ns)
    mar <- floor(((marz[2]-marz[1])-(ns*int))/2)
    slices <- seq(marz[1]+mar, marz[2]-mar, int)
    plotmap("z", slices, 109, 91)
  }
  
  dev.off()
  
  #print legend
  mat <- matrix(0,2,2)
  mat[1,1] <- -mval
  mat[2,2] <- mval
  png(paste(output,"_legend.png",sep=""), height=400, width=500)
  par(mfrow=c(1,1), mar=c(2,2,2,2), oma=c(1,1,1,1))
  rmat <- raster(mat, xmn=1, xmx=200, ymn=1, ymx=200)
  plot(rmat, col=colmap, axes=FALSE, legend=TRUE, box=FALSE)
  dev.off()
}