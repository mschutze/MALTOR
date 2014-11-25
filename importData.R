#### Function to Import Data for Analysis ####
# Manuel Schütze - August 2014

#Libraries required
require(AnalyzeFMRI)

### TEMP ###
mainFolder <- "_age"
input <- list(
  x=list(type="image", path="img", mask="gm_mask.nii"),
  y=list(type="text", path="age.tsv")
)
############

importData <- function(mainFolder,input) {

  #Validate input
  if(typeof(input)!="list") stop("Input is not defined or is not a list")
  if(length(input)<2) stop("Input should contain at least two lists/variables")
  for(i in 1:length(input)) {
    if(typeof(input[[i]])!="list") stop(paste("Input variable",i,"is not defined or is not a list"))
  }
  
  ## Function to read in text ##
  readText <- function(path) {
    cat("Importing text file... ")
    path <- paste(getwd(),mainFolder,path,sep="/")
    if(length(grep("csv",path))>0) { #read text as csv - comma as separator, dot as decimal
      textData <- read.csv(path, header=T, sep=",", dec=".")
    } else { #read text as text file (tab/space separated values)
      textData <- read.table(path, header=T, na.strings=c("\t","-","x","NA"))
    }
    #textData <- as.matrix(textData)
    cat("Read",nrow(textData),"rows and",ncol(textData),"columns.\n")
    returnList <- list(id=textData[,1],data=textData)
    return(returnList)
  }
  
  ## Function to read in images ##
  readImages <- function(path,maskfile) {
    maskfile <- paste(getwd(),mainFolder,maskfile,sep="/")
    
    #first read mask into array, apply formating if needed
    if(length(grep("nii",maskfile))>0) { 
      rawMask <- f.read.nifti.volume(maskfile)  
      } else { 
        rawMask <- f.read.analyze.volume(maskfile) 
    }
    if(any(is.na(rawMask))) {  #remove NaNs
      rawMask <- apply(rawMask, MARGIN=c(1:4), FUN=function(x) ifelse(is.nan(x),0,x)) 
    }
    if(max(rawMask)==255) { #make voxels inside brain mask 0 or 1
      rawMask <- rawMask/255
    }
    maskArray <<- rawMask #X*Y*Z*1 array - global variable
    maskIndex <<- which(maskArray!=0,arr.ind=TRUE) #gets index of non-zero voxels in mask - global variable
    features <- sum(maskArray)
    maskDim <- dim(maskArray)
    
    #find images in folder
    imgFolder <- paste(getwd(),mainFolder,path,sep="/") #full path to image folder
    files <- list.files(path=imgFolder,pattern="*.nii")
    N <- length(files) #number of images

    #read in images into matrix
    imgData <- matrix(0,N,features) #create DATA matrix and fill it with zeros
    cat(paste("Importing",N,"images... "))
    for (i in 1:N) {
      cat(i,ifelse(i==N,".",","),sep="") #prints current image being processed
      currentFile <- paste(imgFolder,files[i],sep="/")
      currentVolume <- f.read.nifti.volume(currentFile)
      if(all(maskDim!=dim(currentVolume))) { 
        stop(paste("Image dimensions for image",i,"don't match mask dimensions!"))
      }
      imgData[i,] <- currentVolume[maskIndex] #Get masked voxels from current volume
    }
    cat("\n")
    returnList <- list(id=files,data=imgData)
    return(returnList)
  }
  
  #Read in data according to data type
  output <- input
  nr <- vector()
  for(i in 1:length(input)) {
    if(input[[i]]$type=="text") output[[i]] <- readText(input[[i]]$path)
    if(input[[i]]$type=="image") output[[i]] <- readImages(input[[i]]$path,input[[i]]$mask)
    nr <- c(nr,length(output[[i]]$id))
  }
  
  #Validate some more
  if(max(nr)!=min(nr)) stop("The number of rows in the input variables does not match")
  ids <- matrix(0,nr[1],length(nr))
  for(i in 1:length(output)) ids[,i] <- as.character(output[[i]]$id)
  colnames(ids) <- paste("Var",names(input))
  rownames(ids) <- paste(1:nrow(ids),". ",sep="")
  cat("----------\nThe IDs for the imported variables are:\n")
  print.table(ids)
  cat("----------\n")
  
  #Clean up and return list
  returnOutput <- output
  for(i in 1:length(output)) returnOutput[[i]] <- output[[i]]$data
  return(returnOutput)
}

### TEMP ###
DATA <- importData(mainFolder,input)