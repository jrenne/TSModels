
Phi_x <- .99 * diag(2)
Sigma_x <- diag(2)
Mu_x <- matrix(0,2,1)

nb.periods <- 100

X <- simul.X(Phi_x,Mu_x,Sigma_x,nb.periods)

uncond.Variance.X <- matrix(solve(diag(2*2) - Phi_x %x% Phi_x) %*%
                              c(Sigma_x %*% t(Sigma_x)),2,2)

stdv.i      <- .01
stdv.lambda <- .02

a.dot <- matrix(c(stdv.i/sqrt(uncond.Variance.X[1,1]),0))
a     <- matrix(c(0,stdv.lambda/sqrt(uncond.Variance.X[2,2])))

b.dot <- .02
b     <- .0

H <- 5

res <- compute.condit.Exp(a,b,
                          a.dot,b.dot,
                          Phi_x,Mu_x,Sigma_x,
                          X=cbind(X[1,1,],X[2,1,]),
                          max.h = H)
res0 <- compute.condit.Exp(a*0.00001,b*0.00001,
                          a.dot,b.dot,
                          Phi_x,Mu_x,Sigma_x,
                          X=cbind(X[1,1,],X[2,1,]),
                          max.h = H)

y_RF <- - log(res0$E.n) / t(matrix(1:H,H,nb.periods))
spreads <- (log(res0$E.n) - log(res$E.n)) / t(matrix(1:H,H,nb.periods))

par(mfrow=c(2,2))
plot(X[1,1,],type="l")
plot(X[2,1,],type="l")
abline(h=0,col="grey")
lines(pmax(0,X[2,1,]),col="blue",lwd=2,lty=3)
plot(y_RF[,1],type="l")
lines(y_RF[,H],col="red")
plot(spreads[,1],type="l")
lines(spreads[,H],col="red")




# Check pricing formula based on simulations:

nb.replics <- 10000
t <- 12

Xt <- matrix(X[,1,t],ncol=1)

X <- simul.X(Phi_x,Mu_x,Sigma_x,nb.periods=H,X0 = Xt,nb.replics = nb.replics)
plot(X[1,1,1:H],type="l",ylim=c(min(X[1,,1:H]),max(X[1,,1:H])))
for(i in 1:100){
  lines(X[1,i,1:H])
}

aux1 <- apply(X,c(2,3),function(x){exp(- b.dot - t(x)%*%a.dot -
                                        pmax(0,b + t(x)%*%a))})

aux2         <- apply(aux1,1,cumprod)
simul.prices <- apply(aux2,1,mean)

rbind(simul.prices,res$E.n[t,])

