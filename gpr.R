#' MALTOR - Gaussian Process Regression (GPR)
#' Manuel Schutze - March 2014 - manuels@ufmg.br
#' 
#' Required libraries: kernlab, raster
#'
#' How to use:
#'  change variables in INPUT and click "Source file"
#' 
#' Output: 
#'  folder results/GPR_[description] containing:
#'  - log.txt (script log)
#'  - map.png / map.nii (weight map) or weights.pdf (for table)
#'  - map_th.png / map_th.nii (thresholded map) or weights.pdf (for table)
#'  - map_legend.png / map_th_legend.png (legend for map colors) for image
#'  - plots.pdf (generated plots)

#Setup
require(kernlab)
source("functions.R")
source("importData.R")

### INPUT (change) ####
mainFolder <- "data/sample"
input1 <- list(type="image", path="img", mask="mask.nii")
#input1 <- list(type="table", path="dados.xlsx")
images <- importData(mainFolder,input1)
input2 <- list(type="table", path="gpr.xlsx")
labels <- importData(mainFolder,input2)
description <- "controls_X_patients"
x <- images$data
xid <- images$ids
xtype <- images$type
y <- labels$data
yid <- labels$ids
#######################
  
#Create new results folder (`mainFolder`/results/GPR_`description`)
#if one already exists (append _2, _3 etc at the end)
resultsFolderPath <- newResultsFolder(mainFolder,"GPR",description)

# Sink output to file
sink(paste(resultsFolderPath,"log.txt",sep="/"),split=TRUE)


cat("\n### Running Gaussian Process Regression Script ###\n")

printTime("start")

#Check if input data matches
if(nrow(x)!=nrow(y)) stop("The number of rows in the input variables does not match")
cat("\nThe IDs for the imported variables are:\n")
cat("Order |   X   |   Y  \n")
for(i in 1:nrow(x)) cat(i,"-",xid[i],"|",yid[i],"\n")

#Remove id collumn from x if it contains it (eg x is a table)
if(class(x)=="data.frame") {
  if(colnames(x)[1]=="id") {
    x <- x[,-1]
  }
}

#Remove subjects from data if y collumn select contains blank group cells (eg only using subset of samples)
if(anyNA(y$select)) {
  cat("\nRemoving subjects:",as.character(y$id[is.na(y$select)]),"\n")
  x <- subset(x, !is.na(y$select))
  y <- subset(y, !is.na(y$select))
}

#Select only variables from y (remove collumns id, select, adj, etc)
Yvars <- y[, !grepl("id|select|adj_|nr_", colnames(y)), drop=FALSE]

#get IDS of patients (first collumn)
ids <- as.character(y[,1])

# If adj collumns are present, run linear regression and get residuals
if(any(grepl("adj_",colnames(y)))) {
  X <- as.matrix(x)
  Adj <- colnames(y)[grepl("adj_",colnames(y))]
  adjFormula <- paste("X ~",paste(rep("y$",length(Adj)),Adj,collapse=":",sep=""))
  lmod <- lm(as.formula(adjFormula))
  #the formula will look something like "X ~ y$adj_age" or "X ~ y$adj_age:y$adj_education", which means that
  #it will create a linear modell that tries to explain the values of X at each voxel based on the 
  #values of the variable(s) in `adj`.
  filteredX <- residuals(lmod)
  #by running ML on the residuals, we are basically removing the effect of our variable(s) defined
  #in `adj` prior to our analysis. this is particullarly important if we have a variable (eg age)
  #that has a big impact on the images being tested.
  cat("\nRemoving the effect of:",Adj,"\n")
} else {
  filteredX <- as.matrix(x) #don't do any correction on X
}

# Create some blank variables
ns <- nrow(filteredX)
nf <- ncol(filteredX)
nv <- ncol(Yvars)
plots <- vector(nv,mode='list')
results <- matrix(NA,nv,5)
colnames(results) <- c("Variable","Subj","Correlation","p","p.adj")

### Run for each variable in Y ###
for(v in 1:nv) {
  curVar <- Yvars[v][,1]
  curName <- names(Yvars[v])
  
  #Test if variable is numeric
  if(class(curVar)!="numeric") stop(paste("Variable",curName,"is not numeric!"))
  
  # Compensates for missing values
  if(NA %in% curVar) { #has missing values in NP
    missVals <- which(is.na(curVar),arr.ind=TRUE) #gets missing values
    X.new <- filteredX[-c(missVals),]
    Y.new <- curVar[-c(missVals)]
    ids.new <- ids[-c(missVals)]
    ns.new <- nrow(X.new)
  } else { #has no missing values in NP
    X.new <- filteredX
    Y.new <- curVar
    ids.new <- ids
    ns.new <- nrow(filteredX)
  }
  
  w <- matrix(0,ns.new,nf) #weight matrix
  pred <- 1:ns.new #predicted values
  
  ## Leave-one-out cross-validation for X.new ##
  cat("\nRunning variable",v,"from",nv,":",curName,"\n")
  
  for(i in 1:ns.new) { #For each subject
    cat("Subject",i,"of",ns.new)
    # Train model
    trainX <- X.new[-i,]
    trainY <- Y.new[-i]
    # Try running GP
    gp <- try(gausspr(trainX,trainY,kernel="vanilladot",variance.model=TRUE))
    # If there is an error, break out of loop
    if (class(gp) == "try-error") break
    # Predict class probability for current subject
    testX <- X.new[i,,drop=FALSE]
    pred[i] <- predict(gp,testX,type="probabilities")
    # Calculate the weights of each voxel using the model parameters (alpha)
    wtemp<-t(t(alpha(gp))%*%trainX)
    norm_a<-c(sqrt(t(wtemp)%*%wtemp))
    w[i,]<-wtemp/norm_a
  }
  
  # Test for NA in weights
  if(any(is.na(w))) {
    warning(paste("Weight matrix w contains NA. Please check values for variable",curName,"\n"), immediate.=TRUE)
  } else {
    # Calculate correlation
    test <- cor.test(Y.new,pred)
    corv <- as.numeric(round(test$estimate,4))
    corp <- as.numeric(round(test$p.value,4))
    
    cat("\nRESULTS:\n")
    cat("Correlation between observed and predicted values:",corv,"( p =",corp,")")
    
    # Generate MAP if p<0.05
    if(corp < 0.05) {
      cat("\nCreating weight map...\n")
      meanw <- apply(w,2,mean)
      if(xtype == "image") {
        #create nifti volume
        meanw_thresh <- thresholdWeights(as.matrix(X.new),meanw,0.05)
        writeMap(meanw,paste(resultsFolderPath,"/map_",curName,"_",sep=""))
        writeMap(meanw_thresh,paste(resultsFolderPath,"/map_th_0.05_",curName,sep=""))
        #create image
        mriNii <- readNIfTI("inc/ref_mri.nii") #reference mri
        mri <- mriNii@.Data
        views <- c("x","y","z") #views to print
        sl <- c(2,8) #rows,colums for each view
        transp <- 70 #0=transparent and 100=opaque
        printMap(views,sl,transp,meanw,mri,paste(resultsFolderPath,"/map_",curName,sep=""))
        printMap(views,sl,transp,meanw_thresh,mri,paste(resultsFolderPath,"/map_th_0.05_",curName,sep=""))
      }
      if(xtype == "table") {
        weights <- cbind(colnames(X),as.numeric(meanw),abs(as.numeric(meanw)))
        orderedWeights <- weights[order(weights[,3]),] #order weights according to absolute value
        graphics.off()
        pdf(paste(resultsFolderPath,"/weights_",curName,".pdf",sep=""), width=8, height=4)
        par(mar=c(10,5,3,2))
        plot(orderedWeights[,2], type="hist", panel.first=abline(v=c(1:nf),col="#DDDDDD"),
             main=paste("GPR for",curName), ylab="Weights", xaxt="n", xlab="",
             col=ifelse(orderedWeights[,2]>0,"blue","red"), lwd=2)
        abline(h=0)
        axis(1, at=c(1:nf), labels=orderedWeights[,1], las=2, cex.axis=0.6)
        dev.off()
        graphics.off()
      }
    }
    
    # Save results
    results[v,] <- c(curName,ns.new,corv,corp,"")
    
    # Plot results
    graphics.off()
    pdf(paste(resultsFolderPath,"/plot_",curName,".pdf",sep=""))
    par(mar=c(5.5,4.5,3,1))
    plot(Y.new, pred, xlab="Observed Values", ylab="Predicted Values", col="black", cex=0.8, pch=19,
         xlim=c(min(Y.new)-2,max(Y.new)+2), ylim=c(min(pred)-2,max(pred)+2),
         main=paste("GPR for ",curName," (n=",ns.new,")",sep=""),
         sub=paste("Correlation:",corv,"| p:",corp))
    text(Y.new, pred, labels=ids.new, pos=1, col="blue", cex=0.5) #plots patient IDs
    lines(c(min(Y.new),max(Y.new)),c(min(Y.new),max(Y.new)),col="red") #what predicted values should be
    arrows(Y.new,pred,Y.new,Y.new,code=0)
    dev.off()
    graphics.off()
    cat("\n")
  }
}

cat("\nSaving plots to PDF...\n")

#prints plots to pdf
if(TRUE) {
  graphics.off()
  pdf(paste(resultsFolderPath,"/plots_",description,".pdf",sep=""), onefile=TRUE)
  par(mar=c(5.5,4.5,3,1))
  for (p in plots) {
    if(class(p)=="recordedplot") replayPlot(p) 
  }
  dev.off()
  graphics.off()
}

cat("Saving results to txt...\n")

#Print results to txt output file
if(nv>1) {
  results.ordered <- results[order(results[,4]),]
  results.ordered[,5] <- p.adjust(results.ordered[,4], method="bonferroni", n=nv)
} else {
  results.ordered <- results
}
write.table(results.ordered,file=paste(resultsFolderPath,"/results_",description,".txt",sep=""),
            row.names=FALSE)

printTime("end")

#Close sink
sink()