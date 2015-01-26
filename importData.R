#### MALTOR - Function to Import Data for Analysis ####
## Manuel Schütze - August 2014 ##

#Libraries required
require(AnalyzeFMRI)

### TEMP ###
#mainFolder <- "nf1"
#input <- list(
#  x=list(type="image", path="gc_img", mask="gm_mask2.nii"),
#  y=list(type="text", path="dados.tsv")
#)
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
    fullpath <- paste(mainFolder,path,sep="/")
    cat("Importing text file",fullpath,"\n")
    if(length(grep("xlsx",fullpath))>0) { #read text as Excel file (xlsx) - only sheet 1 will be read!
      require("xlsx")
      textData <- read.xlsx(fullpath, sheetIndex=1, as.data.frame=TRUE, header=TRUE)
    } else {
      if(length(grep("csv",fullpath))>0) { #read text as csv - comma as separator, dot as decimal
        textData <- read.csv(fullpath, header=TRUE, sep=",", dec=".")
      } else { #read text as text file (tab/space separated values)
        textData <- read.table(fullpath, header=TRUE, na.strings=c("\t","x","NA"))
      }
    }
    cat("Read",nrow(textData),"rows and",ncol(textData),"columns.\n----------\n")
    returnList <- list(id=textData[,1],data=textData)
    return(returnList)
  }
  
  ## Function to read in images ##
  readImages <- function(path,maskfile) {
    maskfile <- paste(mainFolder,maskfile,sep="/")
    
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
    imgFolder <- paste(mainFolder,path,sep="/") #full path to image folder
    files <- list.files(path=imgFolder,pattern="*.nii")
    N <- length(files) #number of images

    #read in images into matrix
    imgData <- matrix(0,N,features) #create DATA matrix and fill it with zeros
    cat("Importing",N,"images from folder",imgFolder,"\n")
    for (i in 1:N) {
      cat(i,ifelse(i==N,". Done!",","),sep="") #prints current image being processed
      currentFile <- paste(imgFolder,files[i],sep="/")
      currentVolume <- f.read.nifti.volume(currentFile)
      if(all(maskDim!=dim(currentVolume))) { 
        stop(paste("Image dimensions for image",i,"don't match mask dimensions!"))
      }
      imgData[i,] <- currentVolume[maskIndex] #Get masked voxels from current volume
    }
    cat("\n----------\n")
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
  cat("The IDs for the imported variables are:\n")
  print.table(ids)
  cat("----------\n")
  
  #Clean up and return list
  returnOutput <- output
  for(i in 1:length(output)) returnOutput[[i]] <- output[[i]]$data
  return(returnOutput)
}

### TEMP ###
#DATA <- importData(mainFolder,input)