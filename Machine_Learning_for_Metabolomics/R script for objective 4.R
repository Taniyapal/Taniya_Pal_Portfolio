#Objective 4

#Preparing the environment
rm(list=ls())
graphics.off()

#loading libraries
library(caret)
library(class)
library(gmodels)
library(rpart)
library(mixOmics)
library(rpart.plot)
library(mlr)

#loading and preparing data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
bact_count_data=read.table("Bacterial_Counts.csv", sep=",", header=T)
merged <- merge(HPLC_data,bact_count_data, by="row.names")
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
AllData=AllData[ ,-c(1,10)]
TVC_data=AllData[,-c(21,19)]
Pseu_data=AllData[,-c(20,19)]


#TVC
set.seed(90)
train.index <- createDataPartition(TVC_data$TVC, p = .7,
                                   list = FALSE, times = 1)
data.train <- TVC_data[train.index,]
data.test <- TVC_data[-train.index,]

#scaling the data

preProcValues <- preProcess(data.train, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, data.train)
testTransformed <- predict(preProcValues, data.test)

# Fit model using all variables
model.fit <- lm(TVC ~ ., data=trainTransformed)
summary(model.fit)

#Inspecting variable importance
varImp(model.fit)

# k-nearest neighbours

model.fit <- caret::train(TVC ~ ., method='knn', data=trainTransformed,tuneGrid=expand.grid(k=1:20))

#plotting the RMSE values of the knn model against k values
plot(model.fit$results$k, model.fit$results$RMSE, xaxt="n",
     ylab="RMSE", xlab="k", main="RMSE for k-Nearest Neighbours regression")
axis(1, at=1:20)
##We can see the RMSE is lowest at k=14. So, we will take k=14.



test.predictions <- predict(model.fit, testTransformed)
model.rmse <- RMSE(testTransformed$TVC, test.predictions)
plot(test.predictions, main=paste('RMSE:', model.rmse), col="red", type="l",ylim=c(-0.9,1.5))
lines(testTransformed$TVC, type="l", col="blue" )
legend("topright",legend=c("actual", "predicted"), col =c("blue","red"), lty=c(1,1),cex=0.5)


#For 100 iterations

trainIndex1 <- createDataPartition(TVC_data$TVC, p = .7, 
                                   list = FALSE, 
                                   times = 100)                              

for (i in 1:100){
  trainSet <- TVC_data[trainIndex1[,i],]
  testSet <- TVC_data[-trainIndex1[,i],]
  trainCl <- trainSet[,ncol(trainSet)]
  testCl <- testSet[,ncol(testSet)]
  preProcValues <- preProcess(trainSet, method = c("center", "scale"))
  trainTransformed <- predict(preProcValues, trainSet)
  testTransformed <- predict(preProcValues, testSet)
  model.fit <- caret::train(TVC ~ ., method='knn', data=trainTransformed,k=14)
  test.predictions <- predict(model.fit, testTransformed)
  model.rmse <- caret::RMSE(testTransformed$TVC, test.predictions)
}

#plotting and comparing actual and predicted sensory values
plot(test.predictions, main=paste('RMSE:', model.rmse), col="red", type="l",ylim=c(-1.5,1.8))
lines(testTransformed$TVC, type="l", col="blue" )
legend("topright",legend=c("actual", "predicted"), col =c("blue","red"), lty=c(1,1),cex=0.7)


#Random Forest

#preparing the environment
rm(list=ls())
graphics.off()

#loading libraries
library(mlr)
library(mlbench)
library(tidyverse)
library(parallelMap)
library(randomForest)
library(caret)
library(MASS)

#loading and preparing data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
bact_count_data=read.table("Bacterial_Counts.csv", sep=",", header=T)
merged <- merge(HPLC_data,bact_count_data, by="row.names")
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
AllData=AllData[ ,-c(1,10)]
TVC_data=AllData[,-c(21,19)]
Pseu_data=AllData[,-c(20,19)]


#Tuning the model
control <- trainControl(method="repeatedcv", number=100, p=0.7)
rf_tuning <- train(TVC~., data=TVC_data, method="rf", tuneGrid=expand.grid(.mtry=c(6:12)), trControl=control )
plot(x=rf_tuning$results$mtry, y=rf_tuning$results$RMSE, xlab="mtry", ylab="RMSE")
#We got mtry=8 as the optimum hyperparameter after tuning

model.fit <- randomForest(TVC ~ ., data = TVC_data, mtry = 8,
                          importance = TRUE, na.action = na.omit)
plot(model.fit, main="Selection of number of trees")
test.predictions <- predict(model.fit, testTransformed)
model.rmse <- RMSE(testTransformed$TVC, test.predictions)

#over 100 iterations

trainIndex1 <- createDataPartition(TVC_data$TVC, p = .7, 
                                   list = FALSE, 
                                   times = 100)                              

for (i in 1:100){
  trainSet <- TVC_data[trainIndex1[,i],]
  testSet <- TVC_data[-trainIndex1[,i],]
  trainCl <- trainSet[,ncol(trainSet)]
  testCl <- testSet[,ncol(testSet)]
  model.fit <- randomForest(TVC ~ ., data = TVC_data, mtry = 7,
                            importance = TRUE, na.action = na.omit)
  test.predictions <- predict(model.fit, testSet)
  model.rmse <- RMSE(testSet$TVC, test.predictions)
}

plot(test.predictions, main=paste('RMSE:', model.rmse), col="red", type="l")
lines(testSet$TVC, type="l", col="blue" )
legend("topright",legend=c("actual", "predicted"), col =c("blue","red"), lty=c(1,1),cex=0.7)

#using PLS-R algorithm

#Preparing the environment
rm(list=ls())
graphics.off()

#loading libraries
library(tidyverse)
library(caret)
library(pls)

#loading and preparing data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
bact_count_data=read.table("Bacterial_Counts.csv", sep=",", header=T)
merged <- merge(HPLC_data,bact_count_data, by="row.names")
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
AllData=AllData[ ,-c(1,10)]
TVC_data=AllData[,-c(21,19)]
Pseu_data=AllData[,-c(20,19)]

#TVC
set.seed(90)
train.index <- createDataPartition(TVC_data$TVC, p = .7,
                                   list = FALSE, times = 1)
data.train <- TVC_data[train.index,]
data.test <- TVC_data[-train.index,]

#scaling the data

preProcValues <- preProcess(data.train, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, data.train)
testTransformed <- predict(preProcValues, data.test)

#generating the model
model <- train(TVC~., data = trainTransformed, method = "pls", tuneLength=10)
# Plot model RMSE vs different values of components
plot(model, xlab="number of components", ylab="RMSE")
# Print the best tuning parameter ncomp that minimize the cross-validation error, RMSE
model$bestTune

##We see here ncomp=1 has lowest value of RMSE

# Make predictions 
trainIndex1 <- createDataPartition(TVC_data$TVC, p = .7, 
                                   list = FALSE, 
                                   times = 100)                              

for (i in 1:100){
  trainSet <- TVC_data[trainIndex1[,i],]
  testSet <- TVC_data[-trainIndex1[,i],]
  trainCl <- trainSet[,ncol(trainSet)]
  testCl <- testSet[,ncol(testSet)]
  
  model <- train(TVC~., data = trainTransformed, method = "pls",ncomp=1)
  test.predictions <- predict(model, testTransformed)
  model.rmse <- RMSE(testTransformed$TVC, test.predictions)
}

predictions <- model %>% predict(testTransformed)
# Model performance metrics

RMSE = caret::RMSE(predictions, testTransformed$TVC)


plot(predictions, type="l", col="red",main=paste('RMSE:', RMSE), ylim=c(-2.0, 2.0))
lines(testTransformed$TVC, type="l", col="blue")
legend("topright", legend=c("actual", "predicted"), col=c("blue","red"), lty=c(1,1), cex=0.4)

#Pseudomonads

#Preparing the environment
rm(list=ls())
graphics.off()

#loading libraries
library(caret)
library(class)
library(gmodels)
library(rpart)
library(mixOmics)
library(rpart.plot)
library(mlr)

#loading and preparing data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
bact_count_data=read.table("Bacterial_Counts.csv", sep=",", header=T)
merged <- merge(HPLC_data,bact_count_data, by="row.names")
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
AllData=AllData[ ,-c(1,10)]
TVC_data=AllData[,-c(21,19)]
Pseu_data=AllData[,-c(20,19)]

set.seed(90)
train.index <- createDataPartition(Pseu_data$Pseudomonads, p = .7,
                                   list = FALSE, times = 1)
data.train <- Pseu_data[train.index,]
data.test <- Pseu_data[-train.index,]

#scaling the data

preProcValues <- preProcess(data.train, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, data.train)
testTransformed <- predict(preProcValues, data.test)

# Fit model using all variables
model.fit <- lm(Pseudomonads ~ ., data=trainTransformed)
summary(model.fit)

#Inspecting variable importance
varImp(model.fit)

# k-nearest neighbours

model.fit <- caret::train(Pseudomonads ~ ., method='knn', data=trainTransformed,tuneGrid=expand.grid(k=1:20))

#plotting the RMSE values of the knn model against k values
plot(model.fit$results$k, model.fit$results$RMSE, xaxt="n",
     ylab="RMSE", xlab="k", main="RMSE for k-Nearest Neighbours regression")
axis(1, at=1:20)
##We can see the RMSE is lowest at k=16. So, we will take k=16.



test.predictions <- predict(model.fit, testTransformed)
model.rmse <- RMSE(testTransformed$Pseudomonads, test.predictions)
plot(test.predictions, main=paste('RMSE:', model.rmse), col="red", type="l",ylim=c(-0.9,1.5))
lines(testTransformed$Pseudomonads, type="l", col="blue" )
legend("topright",legend=c("actual", "predicted"), col =c("blue","red"), lty=c(1,1),cex=0.5)


#For 100 iterations

trainIndex1 <- createDataPartition(Pseu_data$Pseudomonads, p = .7, 
                                   list = FALSE, 
                                   times = 100)                              

for (i in 1:100){
  trainSet <- Pseu_data[trainIndex1[,i],]
  testSet <- Pseu_data[-trainIndex1[,i],]
  trainCl <- trainSet[,ncol(trainSet)]
  testCl <- testSet[,ncol(testSet)]
  preProcValues <- preProcess(trainSet, method = c("center", "scale"))
  trainTransformed <- predict(preProcValues, trainSet)
  testTransformed <- predict(preProcValues, testSet)
  model.fit <- caret::train(Pseudomonads~ ., method='knn', data=trainTransformed,k=2)
  test.predictions <- predict(model.fit, testTransformed)
  model.rmse <- caret::RMSE(testTransformed$Pseudomonads, test.predictions)
}

#plotting and comparing actual and predicted sensory values
plot(test.predictions, main=paste('RMSE:', model.rmse), col="red", type="l",ylim=c(-1.5,1.8))
lines(testTransformed$Pseudomonads, type="l", col="blue" )
legend("topright",legend=c("actual", "predicted"), col =c("blue","red"), lty=c(1,1),cex=0.7)

#Random Forest

#preparing the environment
rm(list=ls())
graphics.off()

#loading libraries
library(mlr)
library(mlbench)
library(tidyverse)
library(parallelMap)
library(randomForest)
library(caret)
library(MASS)

#loading and preparing data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
bact_count_data=read.table("Bacterial_Counts.csv", sep=",", header=T)
merged <- merge(HPLC_data,bact_count_data, by="row.names")
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
AllData=AllData[ ,-c(1,10)]
TVC_data=AllData[,-c(21,19)]
Pseu_data=AllData[,-c(20,19)]


set.seed(90)
train.index <- createDataPartition(Pseu_data$Pseudomonads, p = .7,
                                   list = FALSE, times = 1)
data.train <- Pseu_data[train.index,]
data.test <- Pseu_data[-train.index,]

#scaling the data

preProcValues <- preProcess(data.train, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, data.train)
testTransformed <- predict(preProcValues, data.test)



#Tuning the model
control <- trainControl(method="repeatedcv", number=100, p=0.7)
rf_tuning <- train(Pseudomonads~., data=Pseu_data, method="rf", tuneGrid=expand.grid(.mtry=c(6:12)), trControl=control )
plot(x=rf_tuning$results$mtry, y=rf_tuning$results$RMSE, xlab="mtry", ylab="RMSE")
#We got mtry=8 as the optimum hyperparameter after tuning

model.fit <- randomForest(Pseudomonads ~ ., data = Pseu_data, mtry = 9,
                          importance = TRUE, na.action = na.omit)
plot(model.fit, main="Selection of number of trees")
test.predictions <- predict(model.fit, testTransformed)
model.rmse <- RMSE(testTransformed$Pseudomonads, test.predictions)

#over 100 iterations

trainIndex1 <- createDataPartition(Pseu_data$Pseudomonads, p = .7, 
                                   list = FALSE, 
                                   times = 100)                              

for (i in 1:100){
  trainSet <- Pseu_data[trainIndex1[,i],]
  testSet <- Pseu_data[-trainIndex1[,i],]
  trainCl <- trainSet[,ncol(trainSet)]
  testCl <- testSet[,ncol(testSet)]
  model.fit <- randomForest(Pseudomonads ~ ., data = Pseu_data, mtry = 9,
                            importance = TRUE, na.action = na.omit)
  test.predictions <- predict(model.fit, testSet)
  model.rmse <- RMSE(testSet$Pseudomonads, test.predictions)
}

plot(test.predictions, main=paste('RMSE:', model.rmse), col="red", type="l")
lines(testSet$Pseudomonads, type="l", col="blue" )
legend("topright",legend=c("actual", "predicted"), col =c("blue","red"), lty=c(1,1),cex=0.7)

#using PLS-R algorithm

#Preparing the environment
rm(list=ls())
graphics.off()
HPLC_data=read.table("HPLC_data.csv", sep=",", header=T)
bact_count_data=read.table("Bacterial_Counts.csv", sep=",", header=T)
merged <- merge(HPLC_data,bact_count_data, by="row.names")
AllData<-as.data.frame(merged[,-1])
rownames(AllData) = AllData[,1]
AllData=AllData[ ,-c(1,10)]
TVC_data=AllData[,-c(21,19)]
Pseu_data=AllData[,-c(20,19)]


#loading libraries
library(tidyverse)
library(caret)
library(pls)

#loading and preparing data
setwd("/Users/taniyapal/Downloads/ML ASSIGNMENT DATA")


set.seed(90)
train.index <- createDataPartition(Pseu_data$Pseudomonads, p = .7,
                                   list = FALSE, times = 1)
data.train <- Pseu_data[train.index,]
data.test <- Pseu_data[-train.index,]

#scaling the data

preProcValues <- preProcess(data.train, method = c("center", "scale"))
trainTransformed <- predict(preProcValues, data.train)
testTransformed <- predict(preProcValues, data.test)


#generating the model
model <- train(Pseudomonads~., data = data.train, method = "pls", tuneLength=10)
# Plot model RMSE vs different values of components
plot(model, xlab="number of components", ylab="RMSE")
# Print the best tuning parameter ncomp that minimize the cross-validation error, RMSE
model$bestTune

##We see here ncomp=1 has lowest value of RMSE

# Make predictions 
trainIndex1 <- createDataPartition(Pseu_data$Pseudomonads, p = .7, 
                                   list = FALSE, 
                                   times = 100)                              

for (i in 1:100){
  trainSet <- Pseu_data[trainIndex1[,i],]
  testSet <- Pseu_data[-trainIndex1[,i],]
  trainCl <- trainSet[,ncol(trainSet)]
  testCl <- testSet[,ncol(testSet)]
  
  model <- train(Pseudomonads~., data = data.train, method = "pls",ncomp=1)
  test.predictions <- predict(model, testTransformed)
  model.rmse <- RMSE(testTransformed$Pseudomonads, test.predictions)
}

predictions <- model %>% predict(testTransformed)
# Model performance metrics

RMSE = caret::RMSE(predictions, data.test$Pseudomonads)


plot(predictions, type="l", col="red",main=paste('RMSE:', RMSE))
lines(data.test$Pseudomonads, type="l", col="blue")
legend("topright", legend=c("actual", "predicted"), col=c("blue","red"), lty=c(1,1), cex=0.4)


