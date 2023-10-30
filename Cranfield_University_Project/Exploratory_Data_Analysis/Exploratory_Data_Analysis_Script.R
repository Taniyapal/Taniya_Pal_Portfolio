
#clearing the environment
rm(list=ls())

#clearing all plots
graphics.off()

#setting working directory
setwd("/Users/taniyapal/Documents/Module 2 Repeat Assignment")

###Assignment Task 1

#loading the dataset

ASP=read.csv("A_CS.csv")

##Quality control and initial exploratory data analysis

#converting ASP data frame to numeric matrix

ASP.matrix=as.matrix(ASP[,4:ncol(ASP)])
class(ASP.matrix)="numeric"
ASP.matrix=na.omit(ASP.matrix)


#creating histogram for each variable

histplot=function(varname){
  n=which(colnames(ASP.matrix)==varname)
  hist(ASP.matrix[,n], main=paste("Histogram of", varname), xlab=varname)
  
}

histplot("Fructose")
histplot("Glucose")
histplot("Sucrose")
histplot("TSS")
histplot("ABA")
histplot("DPA")
histplot("PA")
histplot("X7OH_ABA")
histplot("Moisture.loss")

#removing NAs

ASP=na.omit(ASP)

#Descriptive Statistics: mean, sd and median
column.names=c()
mean.columns=c()
colnames=colnames(ASP)
for (i in 3:ncol(ASP)){
  mean.columns=c(mean.columns,mean(ASP[,i]))
}
mean.dataset=data.frame(colnames[3:length(colnames)], mean.columns)
colnames(mean.dataset)=c("Variable Name", "Mean")

##creating new dataframe ASPext

#calculating the sum of sugars from mean (23 is the number of observations in each column of sugar so
# multiplying the mean with number of observation to get the sum of each column of sugar)

row.values=mean.dataset$Mean[2:5]
sum_sugar=c()
for (i in row.values){
  sum_sugar=c(sum_sugar, (23*i))
}

#calculating the sum of ABA and metabolites following similar steps as sum_sugar

metab.values=mean.dataset$Mean[6:9]
ABA_metab=c()
for (i in metab.values){
  ABA_metab=c(ABA_metab, (23*i))
}

#merging the two columns into a data frame
ASP_ext=data.frame("sum_sugar"=sum_sugar,"ABA_metab"=ABA_metab)

##Creating box plots

#making the length of column of sum_sugar and ABA_metab same as TSS and Moisture.loss
sugar_rep=rep(mean(sum_sugar), times=19)
metab_rep=rep(mean(ABA_metab), times=19)
sum_sugar=c(sum_sugar, sugar_rep)
ABA_metab=c(ABA_metab, metab_rep)

#merging all the required columns for boxplot into a data frame
boxplot.dataframe=data.frame(sum_sugar, ABA_metab, ASP$TSS, ASP$Moisture.loss)

#plotting the box plot
boxplot=boxplot(boxplot.dataframe, axes=T, names=c("sugar", "ABA", "TSS", "Moisture"), cex.names=0.4)

#separating out the outliers
outliers.in.boxplot=boxplot$out

##Performing two way ANOVA

library(car)

anova_cal=function(var){
  anova2=aov(var~as.factor(ASP$Treatment)*as.factor(ASP$Time),data=ASP)
  #retrieve the residuals
  resid<-anova2$residuals
  #check the residuals distribution
  hist=hist(resid,main="Histogram of residuals",xlab="Residuals")
  x <- seq(-1, 1, by=0.001)
  rmean<-mean(resid);rsd<-sd(resid)
  y <- dnorm(x = x, mean = rmean, sd = rsd)
  lines(x = x, y = 13*y, col = "blue") ##fit a normal curve over the histogram 
 
  #Perform Tukey test
  tukey.test.result=TukeyHSD(anova2)
  
  
  return(tukey.test.result)
}

tukey.result.sugar=anova_cal(sum_sugar)
tukey.result.ABA=anova_cal(ABA_metab)
tukey.result.TSS=anova_cal(ASP$TSS)
tukey.result.moisture=anova_cal(ASP$Moisture.loss)

#segregating the statistical significant interactions from Tukey test result
library(dplyr)
which((tukey.result.sugar$`as.factor(ASP$Treatment):as.factor(ASP$Time)`[,4])>0.8)
which((tukey.result.ABA$`as.factor(ASP$Treatment):as.factor(ASP$Time)`[,4])>0.8)
which((tukey.result.TSS$`as.factor(ASP$Treatment):as.factor(ASP$Time)`[,4])>0.8)
which((tukey.result.moisture$`as.factor(ASP$Treatment):as.factor(ASP$Time)`[,4])>0.8)

#creating the correlation plots
sugar_rep=rep(mean(sum_sugar), times=19)
sum_sugar=c(sum_sugar, sugar_rep)
ASP_matrix=cbind(ASP.matrix, sum_sugar)

par(mfrow=c(2,2))

#regression plot for TSS and Sum of sugars  
mod_11=lm(ASP[,13]~ASP[,7], data=ASP) #sum_of_sugars(dependant)~TSS(independent)
modsum_11=summary(mod_11)
r2 = modsum_11$adj.r.squared
my.p = modsum_11$coefficients[2,4]
plot(ASP[,7], ASP[,13], xlab="sum of sugars", ylab="TSS", main=c(paste(c("R-squared_value=",r2)), paste(c("p-value",my.p))))
abline(mod_11, col="blue")


#regression plot for TSS and Fructose
mod_1=lm(ASP[,4]~ASP[,7], data=ASP)
modsum_1=summary(mod_1)
r2 = modsum_1$adj.r.squared
my.p = modsum_1$coefficients[2,4]
plot( ASP[,7],ASP[,4], xlab="Fructose", ylab="TSS",main=c(paste(c("R-square=",r2)), paste(c("p-value",my.p))))
abline(mod_1, col="blue")



#regression plot for TSS and Glucose 
mod_2=lm(ASP[,5]~ASP[,7], data=ASP)
modsum_2=summary(mod_2)
r2 = modsum_2$adj.r.squared
my.p = modsum_2$coefficients[2,4]
plot( ASP[,7],ASP[,5], xlab="Glucose", ylab="TSS", main=c(paste(c("R-square=",r2)), paste(c("p-value",my.p))))
abline(mod_2, col="blue")



#regression plot for TSS and Sucrose 
mod_3=lm(ASP[,6]~ASP[,7], data=ASP)
modsum_3=summary(mod_3)
r2 = modsum_3$adj.r.squared
my.p = modsum_3$coefficients[2,4]
plot( ASP[,7],ASP[,6], xlab="Sucrose", ylab="TSS",main=c(paste(c("R-square=",r2)), paste(c("p-value",my.p))))
abline(mod_3, col="blue")

##PCA

#performing PCA
library(FactoMineR)
ASP.pca=PCA(ASP[4:ncol(ASP)], scale=T)

#plotting individual and biplots
library(factoextra)
fviz_pca_ind(ASP.pca, geom="point", geom.ind="text", repel=T, habillage=as.factor(c(ASP$Treatment)))
fviz_pca_biplot(ASP.pca, geom="point", geom.ind="text", repel=T, habillage=as.factor(c(ASP$Treatment)), pallete=c("red", "green", "brown", "blue"), addEllipses = T, ggtheme=theme_minimal()) 

fviz_pca_ind(ASP.pca, geom="point", geom.ind="text", repel=T, habillage=as.factor(c(ASP$Time)))
fviz_pca_biplot(ASP.pca, geom="point", geom.ind="text", repel=T, habillage=as.factor(c(ASP$Time)), pallete=c("red", "green", "brown", "blue"), addEllipses = T, ggtheme=theme_minimal()) 

#plotting scree plot
fviz_screeplot(ASP.pca)

#plotting bar plots of contribution of variables in each PCs
barplot(ASP.pca$var$contrib[, "Dim.1"], names.arg=colnames(ASP[,4:ncol(ASP)]), axes=F, cex.names = 0.4,main="Variable contribution score for PCA Dimension 1" , col=c("blue", "green", "red", "yellow", "orange", "brown", "grey", "black", "light blue"))
barplot(ASP.pca$var$contrib[, "Dim.2"], names.arg=colnames(ASP[,4:ncol(ASP)]), axes=F, cex.names = 0.4,main="Variable contribution score for PCA Dimension 2" , col=c("blue", "green", "red", "yellow", "orange", "brown", "grey", "black", "light blue"))

##HCA

#loading the libraries
library(factoextra)
library(gplots)
library(RColorBrewer)
library(pvclust)
library(dendextend)

#creating function for HCA 
hcluster <- function(X, distance, linkage, title) {
  
  # Plot dendrogram
  result <- pvclust(X, method.dist=distance, 
                    method.hclust=linkage, nboot=100)
  # pvclust and dendextend
  # result %>% as.dendrogram %>% plot
  # #set("branches_k_color", k = k, value = c(1:k)) %>%
  # result %>% text
  # result %>% pvrect
  # cat(c("HCA for ", distance, " distance and ", linkage, " linkage for ", title, " data."))
  # 
  plot(result, main=title)
  pvrect(result, alpha=0.95)
}	
#generating a matrix containing scores for first 3 PCA dimensions
pca_matrix=ASP.pca$var$contrib
class(pca_matrix)="numeric"

#generating dendograms
hcluster(pca_matrix, "euclidean", "ward.D", "mean-centered")
hcluster(pca_matrix, "manhattan", "single", "mean-centered")
hcluster(pca_matrix, "euclidean", "single", "mean-centered")
hcluster(pca_matrix, "manhattan", "ward.D", "mean-centered")

#k-means on raw data

#kmeans cluster on raw data, k=3
x11()
km.res <- kmeans(t(ASP.matrix[,1:9]), 3, nstart = 10)
fviz_cluster(km.res, data=t(ASP.matrix[,1:9]), geom="text", repel=T)

#kmeans cluster on scaled data, k=3
Xm<-scale(t(ASP.matrix[,1:9]), center=T, scale=F)
km.res_scaled <- kmeans(Xm, 3, nstart = 25)
fviz_cluster(km.res_scaled, data=Xm, repel=T)

###Task 2

##loading the SNPs.csv file

snps=read.csv("SNPs_21.csv")

##converting the columns to factors

for (i in 1:ncol(snps)){
  snps[,i]=as.factor(snps[,i])
}

##applying mca on data

#loading the required libraries
library(FactoMineR)
library(factoextra)

#performing the mca
snps.mca=MCA(snps[,2:ncol(snps)])

#plotting the mca individual scores plot
fviz_mca_ind(snps.mca, geom.ind = "", habillage=snps$BMI, pallete=c("#00AFBB", "#E7B800"), repel=F)

#plotting the scores biplot 
fviz_mca_biplot(snps.mca, geom.var=c("point", "text"), habillage=snps$BMI,geom.ind="text",repel = TRUE, contrib=25,addEllipses = T,                
                palette  =  c("#00AFBB",  "#E7B800" ,"#00045"),  ggtheme  = 
                  theme_minimal()) 

###Task 3

#loading the data
df=read.csv("Food_Poisoning.csv")

#loading required packages
require(signal)
require(zoo)
require(ggplot2)

df.ma <- rollmean(df,k = 5) #apply Moving Average on raw data

df.ma_2<- rollmean(df, k=11)

#plotting raw data
plot(df, type="l", ylab="number of people", xlab="year", main="Numbers of people hospitalized from 1950 to 2008")

#filtering the data frame with a threshold of 600
filtered_df=df[which(df$X18>600),]

#plotting the filtered data frame
plot(filtered_df, type="l", ylab="number of people", xlab="year", main="Numbers of people hospitalized from 1950 to 2008")

#Draw raw, Moving Average at wzise 5 and 11 graphs superimposed
plot(df, type="l", xaxt="n", ylab="number of people", xlab="year", main="Numbers of people hospitalized from 1950 to 2008")
axis(1, df[,1], at = seq(1950,2008 ,15),labels=seq(1950, 2008, 15), lty = "solid")
lines(df.ma, lwd=2, lty =1, col="red")
lines(df.ma_2, lwd=2, lty = 2, col="green")
legend(1985, 1400, c("Raw signal", "Centered MA=5","Centered MA=11"),col=c("black", "red", "green"), lty=c(1,1,2), bty="n", cex=0.6)

#kernel smoothing at bandwidth 2
kernel.smoothing.1<-ksmooth(df$X1950, df$X18, 'normal', bandwidth=2)

#kernel smoothing at bandwidth 5
kernel.smoothing.2<-ksmooth(df$X1950, df$X18, 'normal', bandwidth=5)

#plot raw, kernel smoothing at bandwidth 2 and bandwidth 5 plots
plot(df, type="l", xaxt="n", ylab="number of people", xlab="year", main="Numbers of people hospitalized from 1950 to 2008")
axis(1, df[,1], at = seq(1950,2008 ,15),labels=seq(1950, 2008, 15), lty = "solid")
lines(kernel.smoothing.1, lwd=2, lty=1, col="green") 
lines(kernel.smoothing.2, lwd=2, lty=2, col="red")
legend(1985, 1400, c("Raw Signal", "Kernel Smoothing at bandwidth=2", "kernel Smoothing at bandwidth=5"), col=c("black", "green", "red"), lty=c(1,1,2), bty="n", cex=0.48) 

