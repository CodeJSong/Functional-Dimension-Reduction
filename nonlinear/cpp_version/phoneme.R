library(MASS)
library(e1071)

phoneme <- read.csv("phoneme.data")
source("nonlinear.r")
#str(phoneme)

##############################################
# Choose N random sample (N/5 per class) N/5 should be multiple of 2----
##############################################
pclass <- levels(phoneme$g)
# sam.phone - index set of N random sample
nclass <- length(pclass)
N <- 2000   ## multiple of 5*2    
Ni <- c(1,(N/nclass),(N/nclass)+1,(N/nclass)*2,(N/nclass)*2+1,(N/nclass)*3,(N/nclass)*3+1,(N/nclass)*4,(N/nclass)*4+1,(N/nclass)*5)
sam.phone <- rep(0,N)
sam.phone[Ni[1]:Ni[2]] <- sample(which(phoneme$g=="sh"), N/nclass)
sam.phone[Ni[3]:Ni[4]] <- sample(which(phoneme$g=="iy"), N/nclass)
sam.phone[Ni[5]:Ni[6]] <- sample(which(phoneme$g=="dcl"), N/nclass)
sam.phone[Ni[7]:Ni[8]] <- sample(which(phoneme$g=="aa"), N/nclass)
sam.phone[Ni[9]:Ni[10]] <- sample(which(phoneme$g=="ao"), N/nclass)

#################################################
# Divide it training and test set by n.tuning / n.test 
#################################################
n.tuning <- 150  # 150 per class for training set
train.index <- sam.phone[c(Ni[1]:(n.tuning),Ni[3]:(Ni[2]+(n.tuning)),Ni[5]:(Ni[4]+(n.tuning)),Ni[7]:(Ni[6]+(n.tuning)),Ni[9]:(Ni[8]+(n.tuning)))]
test.index <- sam.phone[-c(Ni[1]:(n.tuning),Ni[3]:(Ni[2]+(n.tuning)),Ni[5]:(Ni[4]+(n.tuning)),Ni[7]:(Ni[6]+(n.tuning)),Ni[9]:(Ni[8]+(n.tuning)))]

train.data <- as.matrix(phoneme[train.index,2:151])
train.class <- phoneme$g[train.index]
test.data <- as.matrix(phoneme[test.index,2:151])
test.class <- phoneme$g[test.index]

#################################################
# Apply to training set ----
#################################################
# adjust the data that is fitted with fdr function
#################################################
ttt <- 1:150/150/20
tt <- NULL
x <- NULL
for(i in 1:nrow(train.data)){
  x[[i]] <- as.numeric(train.data[i,])
  tt[[i]] <- ttt
}
y <- rep(1:5,each=n.tuning)
# Plot for paper
Ni <- c(1,(n.tuning)+1,(n.tuning)*2+1,(n.tuning)*3+1,(n.tuning)*4+1)
for(i in 1:5){
  filename <- paste(train.class[Ni[(i)]],".pdf",sep="")
  pdf(filename)
  plot(x[[Ni[(i)]]],type="l", col=1, main=train.class[Ni[(i)]],xlim=c(0,150),
       #cex.main=1, cex.lab=0.9,cex.axis=0.5,
       ylim=c(0,30), xlab="Frequency", ylab="log-periodogram")
  for(j in (Ni[(i)]+1):(Ni[(i)]+9))
  {
    lines(x[[j]],col=(j-Ni[i]+1))
  }  
  dev.off()
}
#################################################
# Find FDR predictor with trainin set (without tuning)
#################################################
# Estimation with Gaussian kernel
phoneme_fir_gauss <- fdr(x=x, y=y, tx=tt, xftn=TRUE, xtequal=TRUE, yclass=TRUE,  kernel="gaussian")
# Estimation with Brownian kernel
phoneme_fir_brown <- fdr(x=x, y=y, tx=tt, xftn=TRUE, xtequal=TRUE, yclass=TRUE, kernel="brownian")


#################################################
# Get predicted value of Test set
#################################################
ttt <- 1:150/150/20
tt.test <- NULL
x.test <- NULL
for(i in 1:nrow(test.data)){
  x.test[[i]] <- as.numeric(test.data[i,])
  tt.test[[i]] <- ttt
}

aaa<-phoneme_fir_gauss
predicted_fir_gauss <- fdr.predict(x.train=x, x.test=x.test, tt.train=tt,tt.test=tt.test, Mat=aaa$Mat, gt=aaa$gt, gx=aaa$gx,
                                   xftn=TRUE, yftn=FALSE, xtequal=TRUE, ytequal=FALSE, 
                                   dimx=1,
                                   option="fir", kernel="gaussian")


aaa<-phoneme_fir_brown
predicted_fir_brown <- fdr.predict(x.train=x, x.test=x.test, tt.train=tt,tt.test=tt.test, Mat=aaa$Mat, gt=aaa$gt, gx=aaa$gx,
                                   xftn=TRUE, yftn=FALSE, xtequal=TRUE, ytequal=FALSE, 
                                   dimx=1, 
                                   option="fir", kernel="brownian")

#################################################
# Do Classification (simplified version)
#################################################
misrate <- function(pred, class){
  n=length(class)
  return(sum(pred!=class)/n)
}
d = 5 # number of nonlinear SDR predictors to be used for classification

train.data.fir.gauss <- data.frame(phoneme_fir_gauss$pred[,1:d],train.class)
train.data.fir.brown <- data.frame(phoneme_fir_brown$pred[,1:d],train.class)

fx <- "X1"
for(i in 1:(d-1)) fx <- paste(fx,"+X",i+1,sep="")    # for modeling
train.lda.fir.gauss <- lda(as.formula(paste("train.class~",fx)), data=train.data.fir.gauss)
train.lda.fir.brown <- lda(as.formula(paste("train.class~",fx)), data=train.data.fir.brown)


lda.pred.training.fir.gauss <- predict(train.lda.fir.gauss, train.data.fir.gauss)$class
lda.pred.training.fir.brown <- predict(train.lda.fir.brown, train.data.fir.brown)$class

c(misrate(lda.pred.training.fir.gauss, train.class),
  misrate(lda.pred.training.fir.brown , train.class))

vec=c(1,2)
plotpred(phoneme_fir_gauss$pred, n.tuning, vec)
plotpred(phoneme_fir_brown$pred, n.tuning, vec)
