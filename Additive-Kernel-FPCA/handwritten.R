load('ujipenchars_norm.RData')

# Choose Characters

ids <- c(0:9)
nclass <- length(ids)

fst <- (1:60)*2
xdata.index <- NULL
for(i in 1:nclass){
  xdata.index <- c(xdata.index, unique(char.norm$id[char.norm$class==as.character(ids[i])])[fst])
}
n <- length(xdata.index)
xdata <- NULL
for(i in 1:n){
  xdata <- rbind(xdata, char.norm[char.norm$id==xdata.index[i],])
}
xdata<- data.frame(xdata)
xdata$sep="train"
# Test data
x.test.index <- NULL
for(i in 1:nclass){
  x.test.index <- c(x.test.index, unique(char.norm$id[char.norm$class==as.character(ids[i])])[-fst])
}
n.test <- length(x.test.index)
x.test.data <- NULL
for(i in 1:n.test){
  x.test.data <- rbind(x.test.data, char.norm[char.norm$id==x.test.index[i],])
}
x.test.data<- data.frame(x.test.data)
x.test.data$sep = "test"

xtotal = rbind.data.frame(xdata,x.test.data)


n.test
x1 = lapply(split(xtotal,xtotal$id), function(xtotal) as.numeric(as.matrix(xtotal[,3])))
x2 = lapply(split(xtotal,xtotal$id), function(xtotal) as.numeric(as.matrix(xtotal[,4])))
sep = lapply(split(xtotal,xtotal$id), function(xtotal) as.character(as.matrix(xtotal[,6])))
sep = as.character(unlist(lapply(sep,function(x) x[[1]])))
y = lapply(split(xtotal,xtotal$id), function(xtotal) as.character(xtotal[,5]))
y = as.character(unlist(lapply(y,function(x) x[[1]])))
tx <-  lapply(split(xtotal,xtotal$id), function(xtotal) as.numeric(xtotal[,2]))
x=NULL
x$x1 = x1
x$x2 = x2
tt=NULL
tt$x1 = tx
tt$x2 = tx

train.ind=which(sep=="train")
test.ind = which(sep=="test")

tmp=nafpca(x, tt, p=2, nbasis=35,unbalanced=TRUE)
#####################
# Nonlinear AFPCA
#####################
plotpred(tmp$pred,vec=c(1,2), ind.info=y)


train.raw = tmp$pred[train.ind,]
train.y = y[train.ind]

test.raw = tmp$pred[test.ind,]
test.y = y[test.ind]

#####################
d=35

train.data <- data.frame(train.raw[,1:d],y=factor(train.y))
my.formula = class.formula(d)

##############################
# Train
##############################


classifiers = c("lda",  "svm")
n.method = length(classifiers)
train.mis.rates = rep(0, n.method)
for(i in 1:n.method){
  train.out = classify(my.formula, train.data,  method=classifiers[i])
  train.mis.rates[i] = misrate(train.out, train.data$y)
  cat(classifiers[i],' : ',train.mis.rates[i],'\n')
}

class.table(train.out,train.data$y)


##############################
# Test
##############################

test.data <- data.frame(test.raw[,1:d],y=factor(test.y))
my.formula = class.formula(d)
classifiers = c("lda",  "svm")
n.method = length(classifiers)
test.mis.rates = rep(0, n.method)
for(i in 1:n.method){
  test.out = classify(my.formula, train.data, test.data, method=classifiers[i])
  test.mis.rates[i] = misrate(test.out, test.data$y)
  cat(classifiers[i],' : ',test.mis.rates[i],'\n')
}
class.table(test.out,test.data$y)



##############################
##############################
#### FPCA
##############################
##############################
tmp2 = fpca(tmp$ftn)
dim(tmp2$pred)
pred.out = tmp2$pred[,1:35]
plotpred(pred.out,vec=c(1,2),ind.info=y)

train.raw = pred.out[train.ind,]
train.y = y[train.ind]

test.raw = pred.out[test.ind,]
test.y = y[test.ind]


train.data <- data.frame(train.raw[,1:d],y=factor(train.y))
my.formula = class.formula(d)

##############################
# Train
##############################


classifiers = c("lda",  "svm")
n.method = length(classifiers)
train.mis.rates = rep(0, n.method)
for(i in 1:n.method){
  train.out = classify(my.formula, train.data,  method=classifiers[i])
  train.mis.rates[i] = misrate(train.out, train.data$y)
  cat(classifiers[i],' : ',train.mis.rates[i],'\n')
}

class.table(train.out,train.data$y)


##############################
# Test
##############################

test.data <- data.frame(test.raw[,1:d],y=factor(test.y))
my.formula = class.formula(d)
classifiers = c("lda",  "svm")
n.method = length(classifiers)
test.mis.rates = rep(0, n.method)
for(i in 1:n.method){
  test.out = classify(my.formula, train.data, test.data, method=classifiers[i])
  test.mis.rates[i] = misrate(test.out, test.data$y)
  cat(classifiers[i],' : ',test.mis.rates[i],'\n')
}
class.table(test.out,test.data$y)









