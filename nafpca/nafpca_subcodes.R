


plotpred <- function(x, n.tuning, vec, xlab="x", ylab="y", main="", ind.info = NULL){
  n <- nrow(x)
  if(!is.null(ind.info)){
    ind.info <- as.numeric(factor(ind.info))  # change indicator number
    nclass <- length(unique(ind.info))
    plot(x[which(ind.info==1),vec[1]], x[which(ind.info==1),vec[2]],
         xlim=c(min(x[, vec[1]]),max(x[, vec[1]])), 
         ylim=c(min(x[, vec[2]]),max(x[, vec[2]])),
         xlab=xlab, ylab=ylab, main=main, pch=16)
    for(i in 2:(nclass)){
      points(x[which(ind.info==i),vec[1]], x[which(ind.info==i),vec[2]], col=i,pch=16)
    }
  }else {
    if(n %% n.tuning != 0 ) stop('n is not multiple of n.tuning')
    nclass <- n/n.tuning
    plot(x[1:n.tuning, vec[1]],x[1:n.tuning, vec[2]],
         xlim=c(min(x[, vec[1]]),max(x[, vec[1]])), 
         ylim=c(min(x[, vec[2]]),max(x[, vec[2]])),
         xlab=xlab, ylab=ylab, main=main,pch=paste(1))
    for(i in 1:(nclass-1)){
      points(x[(i*n.tuning+1):((i+1)*n.tuning), vec[1]],x[(i*n.tuning+1):((i+1)*n.tuning), vec[2]], col=i+1, pch=paste(i+1))
    }
  }
  
  
}

######################################
# misrate : mis classification rate of total
######################################
misrate <- function(pred, class){
  n=length(class)
  return(sum(pred!=class)/n)
}
eval.ftn = function(the_ftn, tt){
  coef=the_ftn$coef
  basis=the_ftn$basis
  return(t(eval.basis(tt, basis)%*% coef))
}


class.table = function(pred, true){
  n = length(unique(true))
  pred = factor(pred); true=factor(true)
  clss = sort(levels(true))
  out = matrix(0,n+2,n+2)
  
  for(i in 1:n){
    ind = which(true==true[i])
    out[i,1:n] = as.vector(table(true[pred==clss[i]]))
  }
  out[1:n,n+1] = apply(out[1:n,1:n], 1,sum)
  out[n+1,1:n] = apply(out[1:n,1:n], 2,sum)
  for(i in 1:n){
    out[i,n+2] = out[i,i]/out[i,n+1]*100
    out[n+2,i] = out[i,i]/out[n+1,i]*100
  } 
  out[n+2,n+2] = sum(diag(out[1:n,1:n]))/sum(out[n+1,1:n])*100
  return(out) 
}
