#' MALTOR - Support Vector Machines (SVM)
#' Manuel Schutze - March 2014 - manuels@ufmg.br
#' 
#' Required libraries: e1071, raster, MASS
#'
#' How to use:
#'  change variables in INPUT and click "Source file"
#' 
#' Output: 
#'  folder results/SVM_[description] containing:
#'  - log.txt (script log)
#'  - map.png / map.nii (weight map) or weights.pdf (for table)
#'  - map_th.png / map_th.nii (thresholded map) or weights.pdf (for table)
#'  - map_legend.png / map_th_legend.png (legend for map colors) for image
#'  - plots.pdf (generated plots)

#Setup
require(e1071)
source("importData.R")
source("functions.R")

### INPUT (change) ####
mainFolder <- "data/sample"
input1 <- list(type="image", path="img", mask="mask.nii")
images <- importData(mainFolder,input1)
input2 <- list(type="table", path="labels.xlsx")
labels <- importData(mainFolder,input2)
description <- "controls_X_patients"
x <- images$data
xid <- images$ids
xtype <- images$type
y <- labels$data
yid <- labels$ids
#######################

#Create new results folder ('mainFolder'/results/SVM_'description')
#if one already exists (append _2, _3 etc at the end)
resultsFolderPath <- newResultsFolder(mainFolder,"SVM",description)

# Sink output to file
sink(paste(resultsFolderPath,"log.txt",sep="/"),split=TRUE)

cat("\n### Running Suport Vector Machine Script ###\n")

printTime("start")

#Check if input data matches
if(nrow(x)!=nrow(y)) stop("The number of rows in the input variables does not match")
cat("\nThe IDs for the imported variables are:\n")
cat("Order |   X   |   Y  \n")
for(i in 1:nrow(x)) cat(i,"-",xid[i],"|",yid[i],"\n")
cat("\n")

#Remove subjects from data if y contains blank group cells (eg only using subset of samples)
if(anyNA(y$group)) {
  x <- subset(x, !is.na(y$group))
  y <- subset(y, !is.na(y$group))
}

#Remove id collumn from x if it contains it (eg x is a table)
if(class(x)=="data.frame") {
  if(colnames(x)[1]=="id") {
    x <- x[,-1]
  }
}

#Make id collumn in y character if it is factor or number
if(typeof(y[,1])!="character") {
  y[,1] <- as.character(y[,1])
}

#Get info about Data
ns <- nrow(x) #number of subjects in study
nf <- ncol(x) #number of features in x

#Check if each group has the same number of subjects
if((ns %% 2)!=0) stop("The two groups have differing sizes")

#Order Data according to group (and pairs if present)
tempm <- cbind(y,x)
if("pair" %in% names(y)) {
  tempm <- tempm[order(tempm$group,tempm$pair),]
} else {
  tempm <- tempm[order(tempm$group),]
}
X <- tempm[,c((1+ncol(y)):(nf+ncol(y)))]
Y <- tempm[,c(1:(ncol(y)))]

#Create some blank variables
runs <- ns/2
xTest <- matrix(NA,2,nf) #Matrix to store test data
yTrain <- factor(Y$group)
w <- matrix(NA,(ns/2),nf)
class <- as.data.frame(matrix(NA,ns,3))
colnames(class) <- c("id","group","predclass")
class[,1] <- as.character(Y$id)
class[,2] <- Y$group
plots <- vector(2,mode='list')

#Print LOOCV order
cat('LOOCV order:\nOrder |  Group 1  |  Group 2\n')
for(i in 1:(runs)) { 
  cat(i, '-', as.character(Y$id[i]), as.character(Y$group[i]), " | ",
      as.character(Y$id[i+runs]), as.character(Y$group[i+runs]), "\n") 
}

### Run leave-one-out cross validation ###
cat("\nRunning leave-one-out cross-validation... this might take some time!\n")

for(i in 1:runs) {
  cat("Training set",i,"from",runs,"\n")
  
  #Define train data
  xTrain <- X[-c(i,i+runs),] #remove test subjects from data
  yTrain <- factor(Y$group[-c(i,i+runs)]) #remove test subjects from labels
  
  #Build model
  svmmodel <- svm(xTrain, yTrain, type="C-classification", kernel="linear")
 
  #Predict groups
  yPred1 <- predict(svmmodel, X[i,], type="response")
  yPred2 <- predict(svmmodel, X[(i+runs),], type="response")
  
  #Save prediction
  group1Name <- as.character(Y$group[1])
  group2Name <- as.character(Y$group[1+runs])
  class[i,3] <- yPred1
  class[(i+runs),3] <- yPred2
  
  #Calculate the weights of each voxel for this training set using the model parameters (alpha)
  wtemp<-t(t(svmmodel$coefs)%*%svmmodel$SV)
  norm_a<-c(sqrt(t(wtemp)%*%wtemp))
  w[i,]<-wtemp/norm_a
}

cat("\nRESULTS:\n")

#Build confusion matrix, calculate statistics
observedClass <- class[,2]
predictedClass <- factor(class[,3], labels=levels(class[,2]))
cmx <- table(observedClass,predictedClass)
names(dimnames(cmx)) <- c("Observed","Predicted")
sens <- cmx[1,1]/sum(cmx[1,]) #sensibilidy
spec <- cmx[2,2]/sum(cmx[2,]) #specificity
acc.test <- binom.test(c((cmx[1,1]+cmx[2,2]),(cmx[1,2]+cmx[2,1])), p=0.5, 
                       alternative="greater", conf.level=0.95)
acc <- acc.test$estimate #overall accuracy
acc.p <- acc.test$p.value #accuracy p
acc.ci <- acc.test$conf.int #accuracy confidence interval

#Plot confusion matrix
fourfoldplot(cmx, color = c("#CC6666", "#99CC99"),
             conf.level = 0, margin = 1, main = "Confusion Matrix")
plots[[1]] <- recordPlot() #saves plot

#Print info
cat("Subjects: ",sum(cmx),", True group ",group1Name,": ",cmx[1,1],", False group ", group1Name,": ",cmx[1,2],
    ", True group ",group2Name,": ",cmx[2,2],", False group ",group2Name,": ",cmx[2,1],"\n",sep="")
cat("Sensibility for group ", group2Name,": ",round(sens*100,2),"\nSpecificity for group ",group2Name,": ",
    round(spec*100,2),"\nOverall Acc: ",round(acc*100,2)," (CI=",acc.ci[1],"-",acc.ci[2],", p=",round(acc.p,3),")\n",sep="")

#Plot results
resplot<-class[,3]
plot(resplot,c(1:runs,1:runs),col=c(rep("blue",runs),rep("red",runs)),pch=19,yaxt="n",
     main="SVM classification",ylab="",xlab="class probabiliy",
     ylim=c(0.5,runs+0.5), xlim=c(0.8,2.2))
text(resplot,c(1:runs,1:runs),Y$id,cex=0.6,pos=c(rep(4,runs),rep(2,runs)))
#text(resplot,c(1:runs,1:runs),Y$age,cex=0.6,pos=4)
abline(v=1.5)
plots[[2]] <- recordPlot() #saves plot

cat("\nSaving plots to PDF...\n")

#prints plots to pdf
graphics.off()
pdf(paste(resultsFolderPath,"/plots_",description,".pdf",sep=""), onefile=TRUE)
par(mar=c(5.5,4.5,3,1))
for (p in plots) { 
  replayPlot(p) 
}
dev.off()
graphics.off()

if(xtype == "image") {
  #Create weight map for x=image
  cat("Creating weight map...\n")
  #create nifti volume
  meanw <- apply(w,2,mean)
  meanw_thresh <- thresholdWeights(as.matrix(X),meanw,0.05)
  writeMap(meanw,paste(resultsFolderPath,"/map_GPC_",description,sep=""))
  writeMap(meanw_thresh,paste(resultsFolderPath,"/map_th.05_GPC_",description,sep=""))
  #create image
  mrinii <- readNIfTI("inc/ref_mri.nii")
  mri <- mrinii@.Data #reference mri
  views <- c("x","y","z") #views to print
  sl <- c(2,8) #rows,colums for each view
  #views <- c("z") #views to print
  #sl <- c(4,4) #rows,colums for each view
  transp <- 70 #0=transparent and 100=opaque
  printMapOLD(views,sl,transp,meanw,mri,paste(resultsFolderPath,"/map_GPC_",description,sep=""))
  printMapOLD(views,sl,transp,meanw_thresh,mri,paste(resultsFolderPath,"/map_th.05_GPC_",description,sep=""))
}

if(xtype == "table") {
  #Create weight map for x=table
  cat("Creating weight map...\n")
  #Define variables
  meanw <- apply(w,2,mean)
  weightorder <- matrix(NA,length(meanw),3)
  weightorder[,1] <- colnames(x)
  weightorder[,2] <- meanw
  weightorder[,3] <- abs(meanw)
  weightorder <- weightorder[order(weightorder[,3]),] #order weights from small to big
  #Plot
  graphics.off()
  pdf(paste(resultsFolderPath,"/weights_",description,".pdf",sep=""))
  par(mar=c(8,4.5,3,1))
  plot(weightorder[,2], type="h", col=ifelse(weightorder[,2]<0,"red","blue"),
       panel.first=abline(v=c(1:nf),col="#DDDDDD"), lwd=2,
       main="Weights for SVM", xlab="", ylab="Weights", xaxt="n")
  axis(1, at=c(1:nf),labels=weightorder[,1], cex.axis=0.5, las=2)
  abline(h=0)
  dev.off()
  graphics.off()
}

#Finish
printTime("end")

#Close sink
sink()