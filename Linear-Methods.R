
# load lasso 
#

install.packages("lars")
library("lars")


# Compute normalization factor. The normalization factor was set so that the average sum
#squares was equal to 1.
normfact <- function (x) {
  nfact = sqrt(sum(x^2))/(253^(.5))
  return (nfact)
}


# Normalize x with the normalization factors obtained above.
normalize <- function (x, nf) {
  x = t(apply(x,1, function(x) x/nf))
  return(x)
}

# centering with respect to mean mu. Using the sweep command to set up the center
centering <- function (x, mu) {
    x = sweep(x, 2, mu)
  return(x)
}

# MSE (mean squared error) between predicted-y py and true-y ty
mse <- function(py,ty) {
  return(mean((py-ty)^2))
}

#Code to compute ridge regression weight
myridge.fit <- function(X,y,lambda) {
  w= rep(0,dim(X)[2])
  w= lm.ridge(y~X, lambda)
  return(w)
}

#Forward Stepwise Regression

model.type <- list(y=~X, lower=~.)
model.forward <- step(lm(y ~ 1, data=train), model.type, direction="forward")
residual.forward <- sum(model.forward$resid^2)/(model.forward$df.residual)

  
# compute predicted Y with coefficients w and data matrix X
predict.linear <- function (w,X) {
  n<-dim(X)[1]
  y<-c(array(0,dim=c(n,1)))
  for (i in 1:n) {
    y[i]=sum(w*X[i,])
  }
  return(y)
}

# compute mean squared error for linear weight w on data (X,y)
mse.linear <- function(w,X,y) {
  py=predict.linear(w,X)
  return (mse(py,y))
}

# select the best s features from path
features.from.path <- function(path, s) {
  k=0;
  bestj=1
  besti=path[1]
  mys=1;
  kk=rep(0,s);
  for (j in 1:length(path)) {

    if (path[j]>0) {
      k=k+1;
      kk[k]=path[j];
    }
    else {
      ik=which(kk[1:k]==-path[j]);
      kk[ik[1]]=kk[k];
      k=k-1;
    }
    if ((mys<k) & (k<=s)) {
      mys=k
    }
    if (k==mys) {
      besti=kk[1:k];
      bestj=j;
    }
  }
  return (besti);
}


# compute variable importance using F-score
feature.Fscore <- function (X,y) {
  p=dim(X)[2]
  score=rep(0,p)
  fit0 <- lm.ridge(y~1)
  fit1 <- lm.ridge(y~X)
  RSS0 <- sum(residuals(fit0)^2)
  RSS1 <- sum(residuals(fit1)^2)
  score <- (RSS0-RSS1)/(13-1)/RSS1*(253-2)
  return(score)
}

# compute variable importance using correlation

feature.Cor <- function (X,y) {
  p=dim(X)[2]
  score=rep(0,p)
  fit = lm.ridge(y~X)
  score = varImp(fit)
  return(score)
}

# plot training/test error versus lambda and save to outfile
#
myplot <- function(lambda,trnerr,tsterr,outfile) {
  ymin=min(min(trnerr),min(tsterr))
  ymax=max(max(trnerr),max(tsterr))
  
  pdf(file=outfile)
  plot(lambda,trnerr,ylim=c(ymin,ymax),ylab="mean squared error", type="l",log="x",lty=3)
  lines(lambda,tsterr,lty=1)
  legend(x="topleft",legend=c("training error","test error"), lty=c(3,1))
  dev.off()
}

#
# plot training and test error with respect to lambda
# for ridge regression and lasso, and output to pdf file
#
implement1 <- function(trnx, trny, tstx, tsty, lambda) {

  # ridge regression
  
  trnerr=lambda
  tsterr=lambda

  # loop through regualrization parameters
  for (i in 1:length(lambda)) {
    # compute coefficient using ridge regression
    ww<-myridge.fit(trnx,trny,lambda[i])
    # compute training error
    trnerr[i]<-mse.linear(ww,trnx,trny)
    # compute test error
    tsterr[i]<-mse.linear(ww,tstx,tsty)
  }

  # plot training/test error for ridge regression
  outfn=paste("prob1.1","-ridge.pdf",sep="")
  cat(paste("output plot to ",outfn,"\n"))
  myplot(lambda,trnerr,tsterr,outfn)
  
  # lasso regression
  
  trnerr2=lambda
  tsterr2=lambda

  # form lasso model
  lasso=0
  # loop through regualrization parameters
  for (i in 1:length(lambda)) {
    fitrr <- lars(trany, tranx, "lasso")
    fitts <- lars(testy, testx, "lasso")
    coetr <- coef.lars(fit)
    coets <- coef.lars(fit)
    ww=rep(0,dim(trnx)[2])
    trnerr2[i]<- mse.linear(fitrr, tranrx, coetr)
    tsterr2[i]<- mse.linear(fitts, tranrx, coets)
  }

  # plot training/test error for lasso
  outfn=paste("prob1.1","-lasso.pdf",sep="")
  cat(paste("output plot to ",outfn,"\n"))
  myplot(lambda,trnerr2,tsterr2,outfn)
  cat("\n\n")
}


implement2 <- function(trnx, trny) {
  cat("features ranked through F-score:\n")
  varimpF= feature.Fscore(trnx,trny)
  path.F <<- sort( varimpF ,decreasing=TRUE,index.return=TRUE)$ix
  print(path.F)
  cat("---\n\n")

  cat("features ranked through correlation:\n")
  varimpCor= feature.Cor(trnx,trny)
  path.Cor <<- sort( varimpCor ,decreasing=TRUE,index.return=TRUE)$ix
  print(path.Cor)
  cat("---\n\n")

  cat("features ranked through Least Squares coefficients:\n")
  varimpLS= abs(myridge.fit(trnx,trny,1e-10))
  path.LS <<- sort( varimpLS ,decreasing=TRUE,index.return=TRUE)$ix
  print(path.LS)
  cat("---\n\n")

  cat("feature ranked through Ridge coefficients (lambda=1):\n")
  varimpRidge= abs(myridge.fit(trnx,trny,1))
  path.Ridge <<- sort( varimpRidge,decreasing=TRUE,index.return=TRUE)$ix
  print(path.Ridge)
  cat("---\n\n")
}

implement3 <- function(trnx, trny) {
  lars("forward stagewise", trace = FALSE, normalise = TRUE, intercept = TRUE)
  path.forward = c(1:13)
  cat("features ranked through forward feature selection:\n")
  print(path.forward)
  cat("---\n\n")

}

implement4 <- function(trnx,trny) {
  lars("lasso", trace = FALSE, normalise = TRUE, intercept = TRUE)
  coef.lars(trainr, trainx)
  cat("lasso path\n")
  print(path.lasso)
  # find top three features
  f3=features.from.path(path.lasso,3)
  cat("top three features\n")
  print(f3)
  cat("---\n\n")

}
  
# evaluate best k features of path on data (X,y)
eval.path <- function(pa,k,trnx,trny,tstx,tsty) {
  # find the best k features
  best.features=features.from.path(pa,k)
  cat("features=[",best.features,"] ",sep=" ")

  
  # least squares fit on training data
  xp=trnx[,best.features];
  if (is.vector(xp)) {
    dim(xp) <- c(length(xp),1)
  }
  ww=myridge.fit(xp,trny,1e-10)

  # mean-squared error on training data
   cat("train-error=",mse.linear(ww,xp,trny)," ",sep="")

 # mean-squared error on test data
 # fill-in code: replace with mean squared error on test data
   cat("test-error=",mse.linear(ww,xp,trny)," ",sep="")
}

implement5 <- function(trnx,trny,tstx,tsty) {
  for (k in c(1,2,3,4,5)) {
    cat("k=",k,"\n",sep="")

    cat(" F-score: ",sep="")
    eval.path(path.F,k,trnx,trny,tstx,tsty)

    cat(" LS-weight: ",sep="")
    eval.path(path.LS,k,trnx,trny,tstx,tsty)

    cat(" ridge-weight: ",sep="")
    eval.path(path.Ridge,k,trnx,trny,tstx,tsty)

    cat(" foward: ",sep="")
    eval.path(path.forward,k,trnx,trny,tstx,tsty)

    cat(" lasso: ",sep="")
    eval.path(path.lasso,k,trnx,trny,tstx,tsty)

  }
}

#Downloading the Data Set

#Header to label the column names
header<-c("CRIM","ZN","INDUS","CHAS","NOX","RM","AGE","DIS","RAD","TAX","PTRATIO","B","LSTAT","MEDV")
#Downloading the data set
housing<-data.matrix(read.table("housing.data.txt",col.names=header))

#Splitting the training and test data. Also having the last
#column represent the response variable in both cases

#Training data split, first 253 rows.
trnx<-housing[1:253,1:13]
trny<-housing[1:253,14]

#Testing data split, second 253 rows.
tstx<-housing[254:506,1:13]
tsty<-housing[254:506,14]


#Normalization of training and test data.
#Command to set normalization coefficient for each feature
nfact=apply(trnx,2,normfact)
#Verifying that each the mean of the squares is equal to 1.
mean(trnx[1:253, 1]^2)
#Normalizing the Training Data
trnx=normalize(trnx,nfact)
#Normalizing the Test Data with same normalization coefficient as that of
#the training set
tstx=normalize(tstx,nfact)

#Centering the data with respect to the mean. For the y-variable, 
#we just subtracting the mean directly.

mux=apply(trnx,2,mean)
muy=mean(trny)
trnx=centering(trnx,mux)
trny= trny - muy
tstx=centering(tstx,mux)
tsty = tsty - muy

#regularization parameter lambda
lambda=c(1e-4,1e-3,1e-2,0.1,1,10,1e2, 1e3,1e4,1e5,1e6,1e7,1e8)
implement1(trnx,trny,tstx,tsty, lambda)
implement2(trnx,trny)
implement3(trnx,trny)
implement4(trnx,trny)
implement5(trnx,trny,tstx,tsty)


