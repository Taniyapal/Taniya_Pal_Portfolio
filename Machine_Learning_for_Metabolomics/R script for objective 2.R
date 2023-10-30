#Preparing the environment

rm(list=ls())
graphics.off()

#loading the data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
Sensory_data=read.table("Sensory_score.csv", sep=",", header=T)

#Objective 2

#using knn algorithm

#loading the libraries
library(class)
library(gmodels)
library(caret)
library(rpart)
library(mixOmics)
library(rpart.plot)


# Match rows from hplc data to rows from sensory
merged <- merge(HPLC_data, Sensory_data, by="row.names")

#remove row names column
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
#removing duplicate samples column

AllData=AllData[ ,-c(1,21)]

#converting classes of sensory data to factors
AllData$sensory=factor(AllData$sensory, levels = c("1", "2", "3"))

#Partitioning data 
set.seed(8)
trainIndex <- createDataPartition(AllData$sensory, p = .7, 
                                  list = FALSE, 
                                  times = 1)

#Generating Test and Training set
trainSet <- AllData[trainIndex,]
testSet <- AllData[-trainIndex,]

#Separating the classes of train and test set
trainCl <- trainSet[,ncol(trainSet)]
testCl <- testSet[,ncol(testSet)]

#removing class values from training and test set
trainSet <- trainSet[,1:(ncol(trainSet)-1)]
testSet <- testSet[,1:(ncol(testSet)-1)]

#auto scaling the train and test set
preProcValues <- preProcess(trainSet, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, trainSet)
testTransformed <- predict(preProcValues, testSet)

#Generating the model with k ranging from 1 to 20
k.results<-function(n, trS, tstS, trCl, tstCl){
  k.accuracy<-c()
  for (k in 1:n) {
    model.k<-knn(trS, tstS, trCl, k)
    confusion.matrix <- confusionMatrix(model.k, tstCl, positive="3")
    a<-confusion.matrix$overall[1]
    k.accuracy<-c(k.accuracy, a)
  }
  return(k.accuracy)
}

k.values=c(1:20)
k.accuracy<-k.results(20, trainTransformed, testTransformed, trainCl, testCl)
plot(k.values, k.accuracy, type='l', ylab="accuracy")
axis(1, at=1:20)
##From the plot we can see the accuracy is highest at value of k at 6.


#Creating a function to identify misclassified samples
calc_error_rate <- function(predicted.value, true.value){
  return(nrow(trainTransformed[true.value!=predicted.value, ])) 
}
#Testing over 100 data partitions
k.accuracy=c()
misclassified.samples=c()
trainIndex1 <- createDataPartition(AllData$sensory, p = .7, 
                                   list = FALSE, 
                                   times = 100)
for (i in 1:100){
  trainSet <- AllData[trainIndex1[,i],]
  testSet <- AllData[-trainIndex1[,i],]
  trainCl <- trainSet[,ncol(trainSet)]
  testCl <- testSet[,ncol(testSet)]
  trainSet <- trainSet[,1:(ncol(trainSet)-1)]
  testSet <- testSet[,1:(ncol(testSet)-1)]
  model.k<-knn(trainSet, testSet, trainCl, 6)
  cross.table <- CrossTable(testCl, model.k, prop.chisq=FALSE, prop.t=FALSE, prop.c=FALSE, prop.r=FALSE)
  confusion.matrix <- confusionMatrix(model.k, testCl, positive="3")
  a<-confusion.matrix$overall[1]
  k.accuracy<-c(k.accuracy, a)
  m <- calc_error_rate(predicted.value=model.k, true.value=testCl)
  misclassified.samples=c(misclassified.samples,m)
}

#plotting the accuracies

iterations=c(1:100)
plot(iterations, k.accuracy, type="l")
##We see that around the 50th iteration, the accuracy is the highest

#Calculating and plotting the cumulative mean accuracies
cummulative.mean.accuracy = c()

for(i in 1:length(k.accuracy)) {
  cummulative.mean.accuracy <- c(cummulative.mean.accuracy, mean(k.accuracy[1:i]))
}

plot(iterations, cummulative.mean.accuracy, type='l', ylab="cumulative mean accuracy")

#plotting the misclassified samples per class over the 100 iterations

barplot(misclassified.samples, col = "lightblue" , names.arg=c(1:100),xlab="Iterations", ylab="Number of missclassified samples per class")

#using svm algorithm

#preparing the environment
rm(list=ls())
graphics.off()

#loading the libraries
library(caret)
library(kernlab)
library(gmodels)

#loading the data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
Sensory_data=read.table("Sensory_score.csv", sep=",", header=T)

# Match rows from hplc data to rows from sensory
merged <- merge(HPLC_data, Sensory_data, by="row.names")

#remove row names column
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
#removing duplicate samples column

AllData=AllData[ ,-c(1,21)]

#converting classes of sensory data to factors
AllData$sensory=factor(AllData$sensory, levels = c("1", "2", "3"))

#Partitioning data 
set.seed(8)
trainIndex <- createDataPartition(AllData$sensory, p = .7, 
                                  list = FALSE, 
                                  times = 1)

#Generating Test and Training set
trainSet <- AllData[trainIndex,]
testSet <- AllData[-trainIndex,]

#Separating the classes of train and test set
trainCl <- trainSet[,ncol(trainSet)]
testCl <- testSet[,ncol(testSet)]

#removing class values from training and test set
trainSet <- trainSet[,1:(ncol(trainSet)-1)]
testSet <- testSet[,1:(ncol(testSet)-1)]

#auto scaling the train and test set
preProcValues <- preProcess(trainSet, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, trainSet)
testTransformed <- predict(preProcValues, testSet)

#generating the svm model
model.svm<-ksvm(trainCl ~ ., data=trainTransformed, kernel = "rbfdot", C=1)

predicted <- predict(model.svm, testTransformed, type="response") 
cross.table <- CrossTable(testCl, predicted, prop.chisq=FALSE, prop.t=FALSE, prop.c=FALSE, prop.r=FALSE) ##Note here we give the predicted values as argument unlike knn where we give the model
confusion.matrix <- confusionMatrix(predicted, testCl, positive="3")
confusion.matrix
cat('SVM accuracy: ', confusion.matrix$overall[1])

summary(model.svm)



#Creating a function to identify misclassified samples
calc_error_rate <- function(predicted.value, true.value){
  return(nrow(trainTransformed[true.value!=predicted.value, ])) 
}

#Testing over 100 data partitions
accuracies<-c()
misclassified.samples.svm=c()
trainIndex1 <- createDataPartition(AllData$sensory, p = .7, 
                                   list = FALSE, 
                                   times = 100)
for (i in 1:100){
  trainSet <- AllData[trainIndex1[,i],]
  testSet <- AllData[-trainIndex1[,i],]
  trainCl <- trainSet[,ncol(trainSet)]
  testCl <- testSet[,ncol(testSet)]
  trainSet <- trainSet[,1:(ncol(trainSet)-1)]
  testSet <- testSet[,1:(ncol(testSet)-1)]
  model_svm <- ksvm(trainCl ~ ., data=trainSet, C=1)
  predicted <- predict(model_svm, testSet, type="response")
  confusion.matrix <- confusionMatrix(predicted, testCl, positive="3")
  current.accuracy <-confusion.matrix$overall[1]
  accuracies<-c(accuracies, current.accuracy)
  m <- calc_error_rate(predicted.value=predicted, true.value=testCl)
  misclassified.samples.svm=c(misclassified.samples.svm,m)
}

iteration <- as.array(c(1:100))
plot(iteration, accuracies, type='l', ylim=c(0,1))


#Calculate the cummulative mean accuracy for each iteration
cummulative.mean.accuracy = c()
for(i in 1:length(accuracies)) {
  cummulative.mean.accuracy <- c(cummulative.mean.accuracy, mean(accuracies[1:i]))
}

##Plot the cummulative mean accuracy
plot(iteration, cummulative.mean.accuracy, type='l', ylim=c(0, 1))

#plotting the misclassified samples per class over the 100 iterations

barplot(misclassified.samples.svm, col = "lightblue" , names.arg=c(1:100),xlab="Iterations", ylab="Number of missclassified samples per class")




#using random forest model

#Preparing the environment
rm(list=ls())
graphics.off()


#loading libraries
library(mlr)
library(mlbench)
library(tidyverse)
library(parallelMap)

#loading the data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
Sensory_data=read.table("Sensory_score.csv", sep=",", header=T)

# Match rows from hplc data to rows from sensory
merged <- merge(HPLC_data, Sensory_data, by="row.names")

#remove row names column
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
#removing duplicate samples column

AllData=AllData[ ,-c(1,21)]

#converting classes of sensory data to factors
AllData$sensory=factor(AllData$sensory, levels = c("1", "2", "3"))

#Partitioning data 
set.seed(8)
trainIndex <- createDataPartition(AllData$sensory, p = .7, 
                                  list = FALSE, 
                                  times = 1)

#Generating Test and Training set
trainSet <- AllData[trainIndex,]
testSet <- AllData[-trainIndex,]

#Separating the classes of train and test set
trainCl <- trainSet[,ncol(trainSet)]
testCl <- testSet[,ncol(testSet)]

#removing class values from training and test set
trainSet <- trainSet[,1:(ncol(trainSet)-1)]
testSet <- testSet[,1:(ncol(testSet)-1)]

#auto scaling the train and test set
preProcValues <- preProcess(trainSet, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, trainSet)
testTransformed <- predict(preProcValues, testSet)

data<- mutate_all(AllData, as.factor)


#Create task
task<- makeClassifTask(data = data, target = "sensory")

#create RF learner
forest <- makeLearner("classif.randomForest")
#Tune the model
forestParamSpace <- makeParamSet(
  makeIntegerParam("ntree", lower = 300, upper = 300), 
  makeIntegerParam("mtry", lower = 6, upper = 12),
  makeIntegerParam("nodesize", lower = 1, upper = 5),
  makeIntegerParam("maxnodes", lower = 5, upper = 20))

randSearch <- makeTuneControlRandom(maxit = 100)       # Define a random search method with 100 iterations

cvForTuning <- makeResampleDesc("CV", iters = 5)       # Define a 5-fold cross-validation strategy

#Tune the hyperparameters
tunedForestPars <- tuneParams(forest, task = task ,
                              resampling = cvForTuning,
                              par.set = forestParamSpace,
                              control = randSearch)

# train a final model to make a learner with the tuned hyperparameters
tunedForest <- setHyperPars(forest, par.vals = tunedForestPars$x)
tunedForestModel <- mlr::train(tunedForest, task)           

#plot the mean out-of-bag error against tree number to see if we included enough trees
forestModelData <- getLearnerModel(tunedForestModel)
plot(forestModelData)
species <- colnames(forestModelData$err.rate)
legend("topleft", species,
       col = 1:length(species),
       lty = 1:length(species), cex=0.4)

#plotting the importance of variables

barplot(forestModelData$importance, beside = T,legend.text = T , names.arg="variables")
