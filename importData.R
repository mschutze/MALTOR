#' MALTOR - Function to Import Data for Analysis
#' Manuel Schutze - August 2014 - manuels@ufmg.br
#' 
#' Required libraries: oro.nifti, xlsx
#'
#' How to use:
#'  mainFolder <- "data/sample"
#'  input1 <- list(type="image", path="img", mask="mask.nii")
#'  input2 <- list(type="table", path="labels.xlsx")
#'  images <- importData(mainFolder,input1)
#'  labels <- importData(mainFolder,input2)
#' 
#' Output: 
#'  list containing 'ids' and 'data' (eg: images$ids and images$data)

importData <- function(mainFolder,input) {

  #Validate input
  if(typeof(input)!="list") stop("input is not defined or is not a list")
  if(!"type"%in%names(input)) stop("input type is not defined")
  if(!"path"%in%names(input)) stop("input path is not defined")
  
  #Check if file/folder in path exists
  fullpath <- paste(getwd(),mainFolder,input$path,sep="/")
  if(!file.exists(fullpath)) stop(paste("input path does not exist:",fullpath))
  
  #Check if mask file exists
  if("mask"%in%names(input)) {
    maskfile <- paste(getwd(),mainFolder,input$mask,sep="/")
    if(!file.exists(maskfile)) stop(paste("mask file does not exist:",maskfile))
  }
  
  ## Function to read in a table of values ##
  readTable <- function(fullpath) {
    cat("\nImporting table file:",fullpath,"\n")
    if(length(grep("xlsx",fullpath))>0) { 
      #read table as Excel file (xlsx) - only sheet 1 will be read!
      require("xlsx")
      tableData <- read.xlsx(fullpath, sheetIndex=1, as.data.frame=TRUE, header=TRUE)
    } else {
      if(length(grep("csv",fullpath))>0) { 
        #read table as csv - comma as separator, dot as decimal
        tableData <- read.csv(fullpath, header=TRUE, sep=",", dec=".")
      } else { 
        #read table as text file (tab/space separated values)
        tableData <- read.table(fullpath, header=TRUE, na.strings=c("\t","x","NA"))
      }
    }
    cat("Read",nrow(tableData),"rows and",ncol(tableData),"columns. Done!\n")
    #Check if collumn 1 is the id and include in return list
    if(!colnames(tableData)[1]=="id")  {
      ids <- rep(NA,nrow(tableData))
      warning("collumn 1 is not named 'id'", immediate. = TRUE)
    } else {
      ids <-  as.character(unlist(tableData[1]))
    }
    returnList <- list(ids=ids,data=tableData,type=input$type)
    return(returnList)
  }
  
  ## Function to read in Nifti images ##
  readImages <- function(fullpath,maskfile) {
    require("oro.nifti")
    
    #first read mask into array, apply formating if needed
    cat("\nReading in mask file:",maskfile,"\n")
    if(length(grep("nii",maskfile))>0) { 
      maskNii <<- readNIfTI(maskfile) #object containing mask info (global)
      rawMask <- maskNii@.Data  #get volume
    } else { 
      stop("Mask should be in nifti format!")
    }
    if(any(is.na(rawMask))) {  #remove NaNs
      rawMask <- apply(rawMask, MARGIN=c(1:3), FUN=function(x) ifelse(is.nan(x),0,x)) 
    }
    if(max(rawMask)>1) { #make voxels inside brain mask 0 or 1
      rawMask <- apply(rawMask, MARGIN=c(1:3), FUN=function(x) ifelse(x>1,1,0))
    }
    maskArray <- rawMask #X*Y*Z array
    maskIndex <<- which(maskArray!=0,arr.ind=TRUE) #gets index of non-zero voxels in mask (global)
    features <- sum(maskArray)
    maskDim <- dim(maskArray)
    cat("Mask dimensions (xyz):",maskDim,"| Total features:",features,"\n")
    
    #find images in folder
    files <- list.files(path=fullpath,pattern="*.nii")
    N <- length(files) #number of images

    #read in images into matrix
    imgData <- matrix(0,N,features) #create DATA matrix and fill it with zeros
    ids <- 1:N
    cat("Importing",N,"images from folder:",fullpath,"\n")
    for (i in 1:N) {
      cat(i,ifelse(i==N,". Done!\n",","),sep="") #prints current image being processed
      currentFile <- paste(fullpath,files[i],sep="/")
      currentnii <- readNIfTI(currentFile)
      currentVolume <- currentnii@.Data
      if(all(maskDim!=dim(currentVolume))) { 
        stop(paste("Image dimensions for image",i,"don't match mask dimensions!"))
      }
      ids[i] <- unlist(strsplit(files[i],".",fixed=T))[1] #Get image name (id)
      imgData[i,] <- currentVolume[maskIndex] #Get masked voxels from current volume
      if(any(is.na(imgData[i,]))) {  #remove NaNs
        imgData[i,] <- apply(imgData[i,], MARGIN=c(1:3), FUN=function(x) ifelse(is.nan(x),0,x)) 
      }
    }
    ids <- as.character(ids)
    returnList <- list(ids=ids,data=imgData,type=input$type)
    return(returnList)
  }
  
  ## Function to read in FreeSurfer .mgh binary files ##
  # TO DO
  
  #Read in data according to data type
  if(input$type=="table") output <- readTable(fullpath)
  if(input$type=="image") output <- readImages(fullpath,maskfile)
  
  #Return output
  return(output)
}
