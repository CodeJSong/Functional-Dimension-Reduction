source('sub_codes.R')
source('nonlinear_subcodes.R')
source('nonlinearfdr.R')


source('simulation model.R')


mi1_extra <- model.pcaxfyr.1.extra(n=100,nt=10,p=100,sd=1,ssq=0.1)
# Balanced data
x <- makelist(mi1_extra$x.b)
tt <- makelist(mi1_extra$t.b)
y<-mi1_extra$y.b
# Unbalanced data
xub <- makelist(mi1_extra$x.ub)
ttub <- makelist(mi1_extra$t.ub)
yub <- mi1_extra$y.ub

# fgsir with tuning
mi1_fgsir <- fdr(x=x, y=mi1_extra$y.b, tx=tt, xftn=TRUE, xtequal=TRUE, cvt=TRUE, cvx=TRUE, ncv=20)

mi1_tune = c(mi1_fgsir$shx,   mi1_fgsir$shy, mi1_fgsir$ex, mi1_fgsir$et, mi1_fgsir$ey)


# fgsir/fgsave using tuned value  (same result)
mi1firb <- fdr(x=x, y=mi1$y.b, tx=tt, xftn=TRUE, xtequal=TRUE, 
               shx=mi1_tune[1], shy=mi1_tune[2], ex=mi1_tune[3], et=mi1_tune[4], ey=mi1_tune[5])
mi1faveb <- fdr(x=x, y=mi1$y.b, tx=tt, xftn=TRUE, xtequal=TRUE, option="FAVE",
                shx=mi1_tune[1], shy=mi1_tune[2], ex=mi1_tune[3], et=mi1_tune[4], ey=mi1_tune[5])

# dim search
aa1 <- dim.search2(mi1firb,ncvl=10)
aa2 <- dim.search2(mi1faveb,ncvl=10)
aa1$dim;aa2$dim
# Balanced data for test

mi1.test<-model.pcaxfyr.1(n=100,nt=10,p=100,sd=1,ssq=0.1)
x.test <- makelist(mi1.test$x.b)
tt.test <- makelist(mi1.test$t.b)

# predict test
mi1firb.test <- fdr.predict(x.train=x, x.test=x.test, tt.train=tt, tt.test=tt.test, Mat=mi1firb$Mat, gt=mi1firb$gt, gx=mi1firb$gx, et=mi1_tune[4],
                            dimx=1, xftn=TRUE, yftn=FALSE, xtequal=TRUE, ytequal=FALSE, )
mi1faveb.test <- fdr.predict(x.train=x, x.test=x.test, tt.train=tt, tt.test=tt.test, Mat=mi1faveb$Mat, gt=mi1faveb$gt, gx=mi1faveb$gx, et=mi1_tune[4],
                             dimx=1, xftn=TRUE, yftn=FALSE, xtequal=TRUE, ytequal=FALSE, )

# mult.cor.mult.rank
c(mult.cor(mult.rank(mi1.test$true.pred[2,]), mult.rank(mi1firb.test[,1:aa1$dim])),
   mult.cor(mult.rank(mi1.test$true.pred[2,]), mult.rank(mi1faveb.test[,1:aa2$dim])))



###################
## miii-1
##################
n <- 100; nt <- 10 ; p <- 10; mu = 0; sd <- 2; N=50; gamma=7; sigma <- 0.1
miii1_brown <- model.xfyf1(n=100, ntx=10, nty=10, p=100, sd=1, ssq=.1)

# Training
x <- makelist(miii1_brown$x.b)
tx <- makelist(miii1_brown$t.b)
ty <- makelist(miii1_brown$ty.b)
y <- miii1_brown$y.b

miii1_brown_firb <- fdr(x=x, y=y, tx=tx, ty=ty,xftn=TRUE,yftn=TRUE, xtequal=TRUE, ytequal=TRUE, kernel="brownian",
                        cvt=TRUE, cvx=TRUE, ncv=20)
miii1_brown_faveb <- fdr(x=x, y=y, tx=tx, ty=ty,xftn=TRUE,yftn=TRUE, xtequal=TRUE, ytequal=TRUE, option="FAVE", kernel="brownian", 
                         cvt=TRUE, cvx=TRUE, ncv=20)

aa1 <- dim.search2(miii1_brown_firb,ncvl=10)
aa2 <- dim.search2(miii1_brown_faveb,ncvl=10)


aa1$dim;aa2$dim
# Test
miii1_brown.test <-  model.xfyf1(n=100, ntx=10, nty=10, p=100, sd=1, ssq=.1)

x.test <- makelist(miii1_brown.test$x.b)
tx.test <- makelist(miii1_brown.test$t.b)

miii1_brown_firb.test <- fdr.predict(x.train=x, x.test=x.test, tt.train=tx, tt.test=tx.test, Mat=miii1_brown_firb$Mat, gt=miii1_brown_firb$gt, gx=miii1_brown_firb$gx,
                                     kernel="brownian", dimx=1, xftn=TRUE, yftn=TRUE, xtequal=TRUE, ytequal=TRUE )
miii1_brown_faveb.test <- fdr.predict(x.train=x, x.test=x.test, tt.train=tx, tt.test=tx.test, Mat=miii1_brown_faveb$Mat, gt=miii1_brown_faveb$gt, gx=miii1_brown_faveb$gx, 
                                      kernel="brownian", dimx=1, xftn=TRUE, yftn=TRUE, xtequal=TRUE, ytequal=TRUE)


# mult.cor.mult.rank
c(mult.cor(mult.rank(miii1_brown.test$true.pred), mult.rank(miii1_brown_firb.test[,1:aa1$dim])),
  mult.cor(mult.rank(miii1_brown.test$true.pred), mult.rank(miii1_brown_faveb.test[,1:aa2$dim])))
