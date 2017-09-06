library( boot )
library( hadron )
j <- 1
Lt=40

rawifit.tsboot <- array(0, dim=c(2000, 3))
etafit.tsboot <- array(0, dim=c(2000, 3))
etacritical.tsboot <- array(0, dim=c(2000))

chisqfn <- function(par,x,y,fitfn,Cinv,Lt){
  # matrix notation,
  ( y - fitfn(x,par,Lt) ) %*% Cinv %*% ( y - fitfn(x,par,Lt) )
}

linearfitfn <- function(x,par,Lt){
  par[1]+par[2]*x
}

constantfitfn <- function(x,par,Lt){
  par[1]
}


etaresults <- data.frame(eta=double(),
                        rawi=double(),
                        error=double(),
                        stringsAsFactors=FALSE
                       )


mu03results <- data.frame(mu03=double(),
                        rawi=double(),  
                        error=double(),   
                        stringsAsFactors=FALSE
                       )

ratio<-matrix( ncol=3, nrow=20)

#Getting x coordinate for the eta_cr fit

etaresults[[1,1]]<-c(-0.52)
etaresults[[2,1]]<-c(-0.56)
etaresults[[3,1]]<-c(-0.58)


# Getting x coordinate for the mu_03 limit

mu03results[[1,1]]<-0.0224
mu03results[[2,1]]<-0.0316
mu03results[[3,1]]<-0.0387


mj<-0
for (m in c( 520, 560, 580)){
   ij<-0
   mj<-mj+1
   for(i in c(224,316,387)) {

      ij<-ij+1

#reading the current density correlation function and symmetrize it

      filename <- paste("results_wigner_eta_m0.",m,"_M02_0000_mu03_0",i,"_rho1.000000_phi_double_smeared_L16T40_JV3DS3_avescalar", sep="")
      correlationfunctionjv3ds3=readtextcf(filename, T=40, check.t=1,skip=0,sym=FALSE, symmetrise=TRUE)

#taking the forward derivative

      tmp<-sapply(seq(from = 1, by = 1, length = (ncol(correlationfunctionjv3ds3$cf) - 1 ) ),function(a){correlationfunctionjv3ds3$cf[,a+1]-correlationfunctionjv3ds3$cf[,a]})
      correlationfunctionjv3ds3$cf<-tmp

#readint the density density correlation function and symmetrize it

      filename <- paste("results_wigner_eta_m0.",m,"_M02_0000_mu03_0",i,"_rho1.000000_phi_double_smeared_L16T40_DS3DS3_avescalar", sep="")
      correlationfunctionds3ds3=readtextcf(filename, T=40, check.t=1,skip=0,sym=TRUE, symmetrise=TRUE)

#removing its last column, to be able to form the ratio

      tmp<-correlationfunctionds3ds3$cf[,-21]
      correlationfunctionds3ds3$cf<-tmp

#Bootstrapping
      bootstrapjv3ds3<-bootstrap.cf(cf=correlationfunctionjv3ds3, boot.R=2000, boot.l=30, seed=1234)
      bootstrapds3ds3<-bootstrap.cf(cf=correlationfunctionds3ds3, boot.R=2000, boot.l=30, seed=1234)

#Taking the ratio of the 2 correlation function, and bootstrapping

      ratio[,1]<-c(1:20)
      ratio[,2]<-bootstrapjv3ds3$cf0/bootstrapds3ds3$cf0
      for ( k in c(1:20)){
        ratio[[k,3]]<-sqrt(sum(bootstrapjv3ds3$cf.tsboot$t[,k]^2/bootstrapds3ds3$cf.tsboot$t[,k]^2)/2000.-(sum(bootstrapjv3ds3$cf.tsboot$t[,k]/bootstrapds3ds3$cf.tsboot$t[,k])/2000.)^2)
      }

      idx<-c(13:19)
      Cinv<-diag(1.0/(ratio[idx,3])^2)
      pars.init<-c(0.1)
      fitfn <- constantfitfn

# fit a constant to the data and store the result of the fit in "fitresult"
   
      fitresult <- optim(par=pars.init,
                   # lambda function  
                   fn=chisqfn,
                   # here one could pass the first derivative of the chisqfn
                   # to make the fit behave better
                   gr=NULL,
                   # get some debug output from the algorithm
                   # finally, the remaining arguments for the chisqfn
                   x=ratio[idx,1],
                   y=ratio[idx,2],
                   fitfn=fitfn,
                   Cinv=Cinv,
                   Lt=Lt)

   

# Getting y coordinate for the mu_03 limit
      mu03results[[ij,2]]<-fitresult$par

#Do the bootstrapping for the error

      for(k in c(1:2000)) {

        fitresultsbootstrap <- optim(par=pars.init, fn = chisqfn,
                                 x=ratio[idx,1],
                                 y=bootstrapjv3ds3$cf.tsboot$t[k,idx]/bootstrapds3ds3$cf.tsboot$t[k,idx],
                                 fitfn=fitfn,
                                 Cinv=Cinv,
                                 Lt=Lt)
        rawifit.tsboot[k, ij  ] <- fitresultsbootstrap$par

      }

      mu03results[[ij,3]]<-(sqrt(sum(rawifit.tsboot[,ij]^2)/2000.-sum(rawifit.tsboot[,ij])/2000.*sum(rawifit.tsboot[,ij])/2000.))

#End for the for loop mu03
   }

#Determining the mu_03->0 limit

# fit a linear to the data and store the result of the fit in "fitresult"
   fitfn<-linearfitfn
   idx <- which( mu03results$mu03 >= 0.0 & mu03results$mu03 <= 0.05 )
   Cinv <- diag(1.0/(mu03results$error[idx])^2)
   pars.init<-c(0.1, 0.0)
   fitresult <- optim(par=pars.init,
                   # lambda function  
                   fn=chisqfn,
                   # here one could pass the first derivative of the chisqfn
                   # to make the fit behave better
                   gr=NULL,
                   # get some debug output from the algorithm
                   # finally, the remaining arguments for the chisqfn
                   x=mu03results$mu03[idx],
                   y=mu03results$rawi[idx],
                   fitfn=fitfn,
                   Cinv=Cinv,
                   Lt=Lt)

#Printing rawi at finite mu_03
   cat("rawi at finite_mu03\n")
   cat("eta=\n")
   print(m)
   print(mu03results)
   etaresults[[mj,2]]<-fitresult$par[1]

#Computing the error using bootstrap samples

   for(k in c(1:2000)) {
        
      fitresultsbootstrap <- optim(par=pars.init, fn = chisqfn,
                                   x=mu03results$mu03[idx],
                                   y=rawifit.tsboot[k,idx],
                                   fitfn=fitfn,
                                   Cinv=Cinv,
                                   Lt=Lt)
      etafit.tsboot[k, mj  ] <- fitresultsbootstrap$par[1]

   }

   etaresults[[mj,3]]<-sqrt(sum(etafit.tsboot[,mj]^2)/2000.-(sum(etafit.tsboot[,mj])/2000.)^2)

}

#Printing results corresponding to mu03=0
print(etaresults)

#Determining eta_cr corresponding to rawi(eta_cr)=0

fitfn<-linearfitfn
idx <- which( etaresults$eta <= -0.5 & etaresults$eta >= -0.6 )
Cinv <- diag(1.0/(etaresults$error[idx])^2)
pars.init<-c(0.1, 0.0)
fitresult <- optim(par=pars.init,
               # lambda function  
                  fn=chisqfn,
               # here one could pass the first derivative of the chisqfn
               # to make the fit behave better
                  gr=NULL,
               # get some debug output from the algorithm
               # finally, the remaining arguments for the chisqfn
                  x=etaresults$eta[idx],
                  y=etaresults$rawi[idx],
                  fitfn=fitfn,
                  Cinv=Cinv,
                  Lt=Lt)
cat("Value for eta_cr from the ensembles\n")
print(-1.*fitresult$par[1]/fitresult$par[2])

#Estimating the error from the bootstap samples

for(k in c(1:2000)) {

  fitresultsbootstrap <- optim(par=pars.init, fn = chisqfn,
                               x=etaresults$eta[idx],
                               y=etafit.tsboot[k,idx],
                               fitfn=fitfn,
                               Cinv=Cinv,
                               Lt=Lt)
  etacritical.tsboot[k] <- -1.*fitresultsbootstrap$par[1]/fitresultsbootstrap$par[2]

}

etacrerror<-sqrt(sum(etacritical.tsboot^2)/2000.-(sum(etacritical.tsboot)/2000.)^2)
print(etacrerror)



