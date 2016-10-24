#####################################################################
# Data cleaning for Handwritten Data
#####################################################################
# handwrttein character can be downloaded at
# http://archive.ics.uci.edu/ml/machine-learning-databases/uji-penchars/version2/
# Data input
hand_alpha <- scan("~/CloudStation/Codes/Data/handwritten digit/ujipenchars2.txt", character(0))

#################################################################################
# UJI Hand.Characters ------
#################################################################################
N <- sum((hand_alpha=="WORD")) 
char.class <- hand_alpha[which(hand_alpha=="WORD")+1]
char.start <- which(hand_alpha=="NUMSTROKES")+5
char.end <- (which(hand_alpha=="WORD")-1)[-1]
char.end <- c(char.end,length(hand_alpha))
data.char <- NULL
for(i in 1:N){
  hand.char <- hand_alpha[char.start[i]:char.end[i]]
  if(("POINTS" %in% hand.char)){
    rmpoint <- c(which(hand.char=="POINTS"),which(hand.char=="POINTS")+1,which(hand.char=="POINTS")+2)
    hand.char <- (hand.char[-rmpoint])
  }
  rmpoint2 <- which(is.na(as.numeric(hand.char)))
  if(!is.na(rmpoint2[1]))  hand.char <- hand.char[-rmpoint2]
  ni <- length(hand.char)/2
  data.char$id <- c(data.char$id,(rep(i,ni)))
  data.char$tt <- c(data.char$tt,(1:ni/ni/10))
  hand.mat <- (matrix(as.numeric(hand.char), ni, 2, byrow=T))
  data.char$xcord <- c(data.char$xcord, hand.mat[,1])
  data.char$ycord <- c(data.char$ycord, -hand.mat[,2])
  data.char$class <- c(data.char$class, rep(char.class[i],ni))
}
data.ch <- data.frame(data.char)
##################################################################
# Normalize the data ------
##################################################################
# Scale each inividual character by [0, 100] X [0, 180] 
##################################################################
unique.id <- unique(data.ch$id)
N <- length(unique.id)
char.norm <- data.ch
for(i in 1:N){
  target <- char.norm[char.norm$id==i,]
  # normalize x coordinate
  temp <- target$xcord
  scale <- max(temp)-min(temp)
  target$xcord <- (temp-min(temp))/scale * 100
  # normalize y coordinate
  temp <- target$ycord
  scale <- max(temp)-min(temp)
  target$ycord <- (temp-min(temp))/scale * 180
  char.norm[char.norm$id==i,] <- target
}
##################################################################
# Make data file, ydata <-  a, b, c,  xdata <- 2,3,6
##################################################################

# Training data
fst <- (1:60)*2 
xdata.index <- c(unique(char.norm$id[char.norm$class=="2"])[fst],unique(char.norm$id[char.norm$class=="3"])[fst],unique(char.norm$id[char.norm$class=="6"])[fst])
ydata.index <- c(unique(char.norm$id[char.norm$class=="a"])[fst],unique(char.norm$id[char.norm$class=="b"])[fst],unique(char.norm$id[char.norm$class=="c"])[fst])
n <- length(xdata.index)
xdata <- NULL
ydata <- NULL
for(i in 1:n){
  xdata <- rbind(xdata, char.norm[char.norm$id==xdata.index[i],])
  ydata <- rbind(ydata, char.norm[char.norm$id==ydata.index[i],])
}
xdata<- data.frame(xdata)
ydata<- data.frame(ydata)
# Test data
x.test.index <- c(unique(char.norm$id[char.norm$class=="2"])[-fst],unique(char.norm$id[char.norm$class=="3"])[-fst],unique(char.norm$id[char.norm$class=="6"])[-fst])
y.test.index <- c(unique(char.norm$id[char.norm$class=="a"])[-fst],unique(char.norm$id[char.norm$class=="b"])[-fst],unique(char.norm$id[char.norm$class=="c"])[-fst])
n.test <- length(x.test.index)
x.test.data <- NULL
y.test.data <- NULL
for(i in 1:n.test){
  x.test.data <- rbind(x.test.data, char.norm[char.norm$id==x.test.index[i],])
  y.test.data <- rbind(y.test.data, char.norm[char.norm$id==y.test.index[i],])
}
x.test.data<- data.frame(x.test.data)
y.test.data<- data.frame(y.test.data)


######################### Estimation with Three  ###################################
source('nonlinear/nonlinear.R')
n.tuning <- 60
vec <- c(1,2)
three_hand_fir_gauss <- fdr(x=xdata, y=ydata, dimx=2, dimy=2, xftn=TRUE,yftn=TRUE, kernel="gaussian", option="fir")
ppp<-three_hand_fir_gauss
plotpred(ppp$pred, n.tuning, vec, xlab="1st Predictor", ylab="2nd Predictor", main="Training set (f-GSIR)")
three_hand_fave_gauss <- fdr(x=xdata, y=ydata, dimx=2, dimy=2, xftn=TRUE,yftn=TRUE, kernel="gaussian", option="FAVE")
ppp<-three_hand_fave_gauss
plotpred(ppp$pred, n.tuning, vec, xlab="1st Predictor", ylab="2nd Predictor", main="Training set (f-GSAVE)")



