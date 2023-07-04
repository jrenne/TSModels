
Phi_x <- .9 * diag(2)
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
b     <- .02

H <- 10

res <- compute.condit.Exp(a,b,
                          a.dot,b.dot,
                          Phi_x,Mu_x,Sigma_x,X,max.h = 10)
res0 <- compute.condit.Exp(a*0.00001,b*0.00001,
                          a.dot,b.dot,
                          Phi_x,Mu_x,Sigma_x,X,max.h = 10)

y_RF <- - log(res0$E.n) / t(matrix(1:H,H,nb.periods))
spreads <- (log(res0$E.n) - log(res$E.n)) / t(matrix(1:H,H,nb.periods))

par(mfrow=c(2,2))
plot(X[,1],type="l")
plot(X[,2],type="l")
abline(h=0,col="grey")
lines(pmax(0,X[,2]),col="blue",lwd=2,lty=3)
plot(y_RF[,1],type="l")
lines(y_RF[,10],col="red")
plot(spreads[,1],type="l")
lines(spreads[,10],col="red")

