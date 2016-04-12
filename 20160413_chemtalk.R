# example from my talk on IIGW 2016
# we classify hot water based on its hydrochemistry


# Paper title: CLASSIFYING HOT WATER CHEMISTRY: APPLICATION OF MULTIVARIATE STATISTICS 
# Authors: Prihadi Sumintadireja, Dasapta Erwin Irawan, Yuano Rezky, Prana Ugiana Gio, Anggita Agustin, Ali Lukman
# Team leader: Prihadi S.
# Data contributor: EBTKE hot water Gorontalo (Yuanno Rezky)
# Code and analysis: Dasapta Erwin Irawan and Prana Ugi
# Field team: Anggita Agustin, Ali Lukman


# Note: remove the '#' symbol only if you haven't done it previously

# Set up working folder

# install packages
## from CRAN server
#install.packages("cluster") 
#install.packages("readxl") 
#install.packages("rpart")
#install.packages("party")
#install.packages("vegan")
#install.packages("mgcv")
#install.packages("gam")
#install.packages("psych")


## from Bioconductor server
#pcaMethods package from Bioconductor server
#source("http://bioconductor.org/biocLite.R") 
#biocLite("pcaMethods")

# Loading library
library(pcaMethods) # for pcaMethods package
library(cluster) # for cluster analysis
#library(readxl) # for opening data directly from xls format
library(rpart)
library(party)
library(vegan)
library(gam)
library(mgcv)
library(psych)

# Load data
df <- read.csv("data_geothermal.csv", header=T) # load data 
row.names(df) <- df$ID # setting row names
df <- df[2:21]

#attach(df)
str(df) # observe data structure (str)
head(df) #observe data 
tail(df) # observe data
is.na(df) # checking NA if any

#df <- na.omit(df) # omitting NA if any
row.names(df)
str(df)# viewing data structure

# Exploratory using pairs() function
# Assesing data patterns
cor.tab <- cor(df)
write.csv(cor.tab, "cortab.csv")
univar.tab <- describe(df)
write.csv(univar.tab, "univar.csv")

pairs(df[2:6],
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)
pairs(df[7:21],
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)
pairs(df[11:21],
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)

# Run PCA 
# Run PCA (using pcamethods package)
## svdImpute = standard pca, with imputation, standardised, method univariate (uv)
pca <- pca(df, 
           method = "svdImpute", 
           scale = "uv",
           center = T,
           nPcs = 5,
           evalPcs = 1:5)
summary(pca)


## Evaluating results
plotPcs(pca, Pcs = 1:5)
plot(pca, type="lines")
slplot(pca) # default function in pcamethods but not big enough
loadings(pca) # loadings of each variables
write.csv(loadings(pca), "loadings.csv")
pca
scores(pca) # scores of each samples respectively to each variables
dev.off() # to reset plot display

## Plotting results (loadings and scores)
plot.new()
par(mfrow=c(1,2))
plot(loadings(pca), 
     pch = 20,
     main = "Variable loadings",
     sub = "Gorontalo area")
text(loadings(pca), 
     row.names(loadings(pca)),
     cex=0.6, 
     pos=1, 
     col="blue")
abline(v=0, h=0, 
       col = "gray70")

plot(scores(pca), 
     main = "Case scores",
     sub = "Gorontalo area")
abline(v=0, h=0, 
       col = "gray70")
text(scores(pca), row.names(df),
     cex=0.6, 
     pos=2, 
     col="blue", offset=0)
dev.off()

# PCA using'vegan' package
vegan.pca <- rda(df)
vegan.pca
plot(vegan.pca)
biplot(vegan.pca, scaling = -1)
ordilabel(plot(vegan.pca), display="species", font=1, col="gray70") # Add some frame on the label
orditkplot(plot(vegan.pca)) # you see that you can move the labels.

# Cluster analysis (using 'cluster' package)
d <- dist(df, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") # fitting distance to cluster, method Ward
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
summary(fit)

# Cluster anlaysis (using 'vegan' package)
d2 <- vegdist(df)
fit2 <- hclust(d2, "ward.D")
fit3 <- hclust(d2, "complete")
fit4 <- hclust(d2, "average")

par(mfrow=c(1,3))
plot(fit2)
plot(fit3)
plot(fit4)
dev.off()

# Regression tree
tree.fit <- rpart(ec ~ ., data = df, 
                  control=rpart.control(minsplit=2, minbucket=1, cp=0.001))
printcp(tree.fit) # display the results
plotcp(tree.fit) # visualize cross-validation results
summary(tree.fit) # detailed summary of splits


# plot tree
plot(tree.fit, uniform=TRUE,
     main="Tree classification of hot water samples: Gorontalo")
text(tree.fit, use.n=TRUE, all=TRUE, cex=.8)
mod_lm <- gam(ec ~ ., data = df)
summary(mod_lm)


# Ref:
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/CA_lab.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio3.pdf
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/PCA_lab1.pdf
## http://bioconductor.wustl.edu/bioc/vignettes/pcaMethods
## https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
## https://www3.nd.edu/~mclark19/learn/GAMS.pdf


# here we classify some groundwater springs based on their hydrochemistry

#####################################################
#####################################################
#####################################################
#####################################################

# PCA of Cisanti Area
# Data: PKM Project in Cisanti Area, Bandung
# Team leader: Arif Susanto
# Code and analysis: Dasapta Erwin Irawan

# install biocLite script for package installation from 
# biocondictor (not CRAN) repo 
# load library
source("http://bioconductor.org/biocLite.R") 
biocLite("pcaMethods")
library(pcaMethods) # for pcaMethods package
library(cluster) # for cluster analysis
install.packages("readxl")
library(readxl) # for opening data directly from xls format
install.packages('xtable')
library(xtable)

# Load data
# I'm reading from xls (usually from csv) due to 'unicode' failures 
# (don't know exactly what happened)
df <- read.csv("data_cisanti.csv")
#df <- as.data.frame(df<-read_excel("data2.xls")) 
dim(df) # viewing data frame dimension
row.names(df)<-df$ID # setting row names 
rownames(df)

df <- (df[4:33]) # cutting data frame to only numeric values
df <- na.omit(df) # omitting NA if any
is.na(df) # checking NA if any
str(df)# viewing data structure -> all numerics

# Exploratory using pairs() function
# Assesing data patterns
pairs(df[4:10],
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)
pairs(df[11:15],
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)
pairs(df[16:20],
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)
pairs(df[21:25],
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)
pairs(df[26:33],
      lower.panel=panel.smooth, 
      upper.panel=NULL, 
      pch=20)

# Run PCA (done manually)

## Calculating covariance matrix
cov(df)

## Calculating eigenvalues of the covar matrix
eigenvalues<-eigen(cov(df))$values
eigenvectors<-eigen(cov(df))$vectors
mean(eigenvalues)

## Estimating PC's via matrix multiplication = df * eigenvectors
PC <- as.matrix(df) %*% eigenvectors

## Computing the covar matrix of PC. 
## The variances of cov(PC) should be equal to the Eigenvalues and the 
## covariances should be 0 (aside from rounding errors) 
cov(PC)
eigenvalues[1:3]
cov(PC)[1:3, 1:3]

## Calculating the proportions of the variation explained by 
## the various components (see the numbers of the first 3 PC's)
print(round(eigenvalues/sum(eigenvalues) * 100, digits = 2))
round(cumsum(eigenvalues)/sum(eigenvalues) * 100, digits = 2)


# Run PCA (using prcomp)
PCA1 <- prcomp(df)
PCA1

## Extracting the variances of the components 
PCA1.var <- PCA1$sdev^2 
PCA1.var[1:3]

## which are identical to the eigenvalues (in the  first try out) 
eigenvalues[1:3]

## Plotting PC's and rotations
plot(PCA1)
par(mfrow = c(2, 2))
plot(PCA1$rotation[, 1], ylim = c(-0.7, 0.7))
plot(PCA1$rotation[, 2], ylim = c(-0.7, 0.7))
plot(PCA1$rotation[, 3], ylim = c(-0.7, 0.7))
plot(PCA1$rotation[, 4], ylim = c(-0.7, 0.7))

# Run PCA (using princomp)
PCA2 <- princomp(df, cor = F)
# Apparently the 'princomp()' function needs more samples
# We only have 7 LOL

# Run PCA (using pcamethods package)
## svdImpute = standard pca, with imputation, standardised, method univariate (uv)
pca <- pca(df, 
           method = "svdImpute", 
           scale = "uv",
           center = T,
           nPcs = 5,
           evalPcs = 1:5)
summary(pca)


## Evaluating results
plotPcs(pca, Pcs = 1:5)
plot(pca, type="lines")
slplot(pca) # default function in pcamethods but not big enough
loadings(pca) # loadings of each variables
scores(pca) # scores of each samples respectively to each variables
row.names(scores(pca)) <- rownames(df$ID)

## Plotting results (loadings and scores)
plot.new()
par(mfrow=c(1,2))
plot(loadings(pca), 
     pch = 20,
     main = "Variable loadings",
     sub = "Cisanti area")
text(loadings(pca), 
     row.names(loadings(pca)),
     cex=0.6, 
     pos=1, 
     col="red")
abline(v=0, h=0, 
       col = "gray70")

plot(scores(pca), 
     main = "Case scores",
     sub = "Cisanti area")
abline(v=0, h=0, 
       col = "gray70")
text(row.names(df),
     cex=0.6, 
     pos=1, 
     col="red", offset=0)
dev.off()

vegan.pca <- rda(df[4:33])
vegan.pca
plot(vegan.pca)
dev.off()
biplot(vegan.pca, scaling = -1)


# Cluster analysis (using 'cluster' package)
d <- dist(df, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") # fitting distance to cluster, method Ward
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
summary(fit)

# Cluster anlaysis (using 'vegan' package)
install.packages("vegan")
library(vegan)
d2 <- vegdist(df)
fit2 <- hclust(d2, "ward.D")
fit3 <- hclust(d2, "complete")
fit4 <- hclust(d2, "average")

par(mfrow=c(1,3))
plot(fit2)
plot(fit3)
plot(fit4)
dev.off()

# Ref:
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/CA_lab.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio3.pdf
## http://www2.stat.unibo.it/montanari/Didattica/Multivariate/PCA_lab1.pdf
## http://bioconductor.wustl.edu/bioc/vignettes/pcaMethods
## https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
## http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
