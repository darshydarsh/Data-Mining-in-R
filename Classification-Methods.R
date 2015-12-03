
# install and use the MASS library

library("MASS")


# generating isotropic gaussians
# n data for each class 1,2,3 in p-dimension
isogs <- function(n,p) {
  X=matrix(rnorm(3*n*p),nrow=3*n,ncol=p)
  for (i in 1:n) {
    X[i,1]=X[i,1]+3;
    X[n+i,1]=X[n+i,1]-3;
    X[2*n+i,1:2]=X[2*n+i,1:2]+c(1,-1);
  }
  y=c(rep(1,n),rep(2,n),rep(3,n))
  data = list(X=X,y=y)
  return(data)
}

# generate anisotropic gaussians
# n data for each class 1,2,3 in p-dimension
anisogs <- function(n,p) {
  X=matrix(rnorm(3*n*p),nrow=3*n,ncol=p)
  X[,1]=X[,1]*0.5;
  X[,2]=X[,2]*1.5;
  X[,3:p]=X[,3:p]*matrix(runif(3*n*(p-2),min=0,max=3),nrow=3*n,ncol=p-2);
  for (i in 1:n) {
    X[i,1]=X[i,1]+3;
    X[n+i,1]=X[n+i,1]-3;
    X[2*n+i,1:2]=X[2*n+i,1:2]+c(1,-1);
  }
  y=c(rep(1,n),rep(2,n),rep(3,n))
  data = list(X=X,y=y)
  return(data)
}


# scatter plot data to 2-d
# X: data matrix
# y: label from 1-3
# v: v$v1 and v$v2 are two projection directions
# outfn: output file name
scatterplot <- function(X,y,v,outfn) {
  cat(paste("output plot to ",outfn,"\n"))
  Xp=pj(X,v)
  pdf(file=outfn)
  i1=which(y==1,TRUE);
  plot(Xp[i1,],xlim=c(min(Xp[,1]),max(Xp[,1])),ylim=c(min(Xp[,2]),max(Xp[,2])),type="p",pch='1',col='red')
  i2=which(y==2,TRUE);
  points(Xp[i2,],type="p",pch='2',col='green')
  i3=which(y==3,TRUE);
  points(Xp[i3,],type="p",pch='3',col='blue')
  dev.off()
}

# projection of data matrix X to 2-dimension [v$v1,v$v2]
pj <- function(X, v) {
  Xp=matrix(nrow=dim(X)[1],ncol=2)
  Xp[,1]=X%*%matrix(v$v1,ncol=1)
  Xp[,2]=X%*%matrix(v$v2,ncol=1)
  return(Xp)
}

# projection of data$X to 2-dimension [v$v1,v$v2] and align with data$y
pj2 <- function(data,v) {
  Xp=matrix(nrow=length(data$y),ncol=2)
  Xp[,1]=data$X%*%matrix(v$v1,ncol=1)
  Xp[,2]=data$X%*%matrix(v$v2,ncol=1)
  return (list(X=Xp,y=data$y))
}

# generating two orthogonal random-directions
randpj <- function(X) {
  p=dim(X)[2]
  v1=rnorm(p)
  v1=v1/sqrt(sum(v1*v1))
  v2=rnorm(p)
  v2=v2-sum(v1*v2)*v1
  v2=v2/sqrt(sum(v2*v2))
  v=list(v1=v1,v2=v2)
  return(v)
}

# top2-pca-directions, assuming centered X
pcapj <- function(X) {
  #Setting up the singular value decomposition
  sv <- svd(X)
  #Setting up the values of v
  svp <- sv$v
  #Looking at two largest directions
  v1 <- svp[,1]
  v2 <- svp[,2]
  #Returning as a list
  v = list(v1=v1, v2=v2)
  return(v)
}


# top2-lda-directions, assuming centered X
ldapj <- function(data) {
  
  # use lda() function, and return the first dimension as v$v1
  # and second as v$v2

  lv <- lda(data$y~data$X, data) 
  ldp <- predict(lv, data)
  v1 <- ldp$x[,1]
  v2 <- ldp$x[,2]
  v = list(v1=v1[1:200], v2=v2[1:200])
  
  return(v)
}

# one versus all training with generalized linear model
oneversusall.fit <- function(data,fam) {
  # w1: y=1 versus others
  # w2: y=2 versus others
  # w3: y=3 versus others
  
  w1= glm.fit(data$X, data$y[1:200] , family=fam, intercept = TRUE)
  w2= glm.fit(data$X, data$y[1:200] , family=fam, intercept = TRUE)
  w3= glm.fit(data$X, data$y[1:200] , family=fam, intercept = TRUE)
  w=list(w1=w1,w2=w2,w3=w3)
  return(w)
}


# one versus all classification 
predict.oneversusall <- function(w,X) {
  n=dim(X)[1]
  p=dim(X)[2]
  y=rep(1,n)
  for (i in 1:n) {
    xx=X[i,]
    # s1: score for class 1
    # s2: score for class 2
    # s3: score for class 3
    s1= predict.glm(w1, xx)
    s2= predict.glm(w2, xx)
    s3= predict.glm(w3, xx)

    # y is the label with maximum score
    y[i]= max(s1,s2,s3)
  }
  return(y)
}

#classification error between predicted label vector py and true label vector ty
# 
cerr <- function(py,ty) {
  err= sum(py != ty)/length(py)
  return(err)
}

#classification error between predicted label vector py and true label vector ty:
# try to assign the optimal label correspondence btween label values for py (clustering output) and ty
#
cerr.bestmatch <- function(py,ty) {
  K=3
  # find best correspondence between label value in py and ty
  # assuming K=3 classes
  for (i in (1:(K-1))) {
    ii=which(ty==i)
    for (j in ((i+1):K)) {
      if (sum(py[ii]==i) <sum(py[ii]==j)) {
        py[which(py==i)]=(K+1)
        py[which(py==j)]=i
        py[which(py==(K+1))]=j
      }
    }
  }
  return (cerr(py,ty))
}

implementation1 <- function(trn,tst,label) {

  # First Implementation
  #random projection
  #vrand <<- randpj(trn$X)
  #scatterplot(trn$X,trn$y,vrand, paste("prob1-",label,"-rand.pdf",sep=""))

   #pca projection
  vpca <<- pcapj(trn$X)
  cat("first 5 components of 1st PC\n")
  print(vpca)
  print(vpca$v1[1:5])
  cat("first 5 components of 2nd PC\n")
  print(vpca$v2[1:5])
  scatterplot(trn$X,trn$y,vpca,paste("prob1-",label,"-pca-trn.pdf",sep=""))
  scatterplot(tst$X,trn$y,vpca,paste("prob1-",label,"-pca-tst.pdf",sep=""))

   #lda projection
  vlda <<- ldapj(trn)
  cat("first 5 components of 1st LDC\n")
  print(vlda$v1[1:5])
  cat("first 5 components of 2nd LDC\n")
  print(vlda$v2[1:5])
  scatterplot(trn$X,trn$y,vlda,paste("prob1-",label,"-lda-trn.pdf",sep=""))
  scatterplot(tst$X,trn$y,vlda,paste("prob1-",label,"-lda-tst.pdf",sep=""))
}



Implentation2 <- function(trn,tst,label) {

  lvtrn <- lda(trn$y~trn$X, trn)
  lvtrn.pred <- predict(lvtrn, trn)
  fin.predtrn <- lvtrn.pred$class
  pytrn <- vector(,300) 
  for (i in 1:length(trn$y)) {
    pytrn[i]=  fin.predtrn[1]
  }
  err=cerr(pytrn,trn$y)
  cat("lda ",label," error =",err,"\n")
  
  lvtst <- lda(tst$y~tst$X, tst)
  lvtst.pred <- predict(lvtst, tst)
  fin.predtst <- lvtst.pred$class
  pytst <- vector(,300)
  for (i in 1:length(tst$y)) {
    pytst[i]= fin.predtst[1]
  }
  err=cerr(pytst,tst$y)
  cat("lda ",label," error =",err,"\n")
}

Implmentation3 <- function(trn,tst) {
  # least squares 
  w=oneversusall.fit(trn,gaussian())
  err=cerr(predict.oneversusall(w,trn$X),trn$y)
  cat(paste("  least squares training error =",err,"\n"))
  err=cerr(predict.oneversusall(w,tst$X),tst$y)
  cat(paste("  least squares test error =",err,"\n"))
  
  # logistic regression 
  w=oneversusall.fit(trn,binomial())
  err=cerr(predict.oneversusall(w,trn$X),trn$y)
  cat(paste("  logistic training error =",err,"\n"))
  err=cerr(predict.oneversusall(w,tst$X),tst$y)
  cat(paste("  logistic test error =",err,"\n"))
}

 Implementation4 <- function(trn, tst)
{
# use kmeans() to cluster it into three categories, and return py as cluster label
py= kmeans(tst$X,3, iter.max = 10, nstart=1)
# find the best label correspondence to true label tst$y
err=cerr.bestmatch(py,tst$y)
cat(paste("  kmeans clustering error =",err,"\n"))
}

# do the experiments with train data Xtrn and test data Xtst
# label is the dataset: isotropic or anisotropic
#
doexp <- function(trn,tst,label) {
  mu=colMeans(trn$X)
  trn$X=t(t(trn$X)-mu)
  tst$X=t(t(tst$X)-mu)
  
  cat(paste("---",label,"---\n"))

  implementation1(trn,tst,label)

  #implementation2(trn,trn,"training")
  #implementation2(trn,tst,"test")

  # without dimension reduction
  #cat("without dimension reduction\n")
  #implementation3(trn,tst)

  # random
  #cat("random dimension reduction\n")
  #implementation4(pj2(trn,vrand),pj2(tst,vrand))

  # pca
  #cat("pca dimension reduction\n")
  #implementation5(pj2(trn,vpca),pj2(tst,vpca))

  # lda
  #cat("lda dimension reduction\n")
  #implementation5(pj2(trn,vlda),pj2(tst,vlda))
  #cat("\n")
  
  #K-means
  #cat("k-means error")
  #implementation6(trn,tst)
}

#Setting p and n values
p=200
n=100

set.seed(76368)

#Command to create the isotropic gaussians
trn1=isogs(n,p)
tst1=isogs(n,p)
tst1
#Command to create the anisotropic gaussians
trn2=anisogs(n,p)
tst2=anisogs(n,p)

#Problem 1.1, 1.2, 1.3 and 1.4
#Command to have the first two orthogonal directions in R^p using previous functions.
doexp(trn1,tst1,"isotropic")
doexp(trn2,tst2,"anisotropic")

