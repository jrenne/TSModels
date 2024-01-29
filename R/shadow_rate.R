simul.var <- function(Model,nb.sim,x0=NaN){
  n <- dim(Model$rho)[1]
  if(is.na(x0[1])){
    x0 <- solve(diag(n) - Model$rho) %*% Model$mu
  }

  X <- c(x0)
  x <- x0
  for(t in 2:nb.sim){
    x <- Model$mu + Model$rho %*% x + Sigma %*% rnorm(n)
    X <- rbind(X,c(x))
  }
  return(X)
}


compute.price.WX <- function(Model,X,max.H){
  # max.H is the maximum considered maturity, expressed in number of periods.

  SS <- Model$Sigma %*% t(Model$Sigma)

  n <- dim(Model$rho)[1]

  # For h=1:
  sum.rho.j_1 <- diag(n)
  rho.j <- diag(n)

  sigma.h.2 <- 0

  all.a.bar <- NULL
  all.b <- NULL
  all.a <- NULL
  all.sigma <- NULL

  for(h in 1:max.H){
    sigma.h.2 <- sigma.h.2 + matrix(Model$delta.1,nrow=1) %*% rho.j %*% SS %*%
      t(rho.j) %*% matrix(Model$delta.1,ncol=1)

    a.h.bar <- Model$delta.0 + matrix(Model$delta.1,nrow=1) %*% sum.rho.j_1 %*% Model$mu
    a.h <- a.h.bar - .5*matrix(Model$delta.1,nrow=1) %*% sum.rho.j_1 %*%
      SS %*% t(sum.rho.j_1) %*% matrix(Model$delta.1,ncol=1)

    rho.j <- rho.j %*% Model$rho

    b.h <- matrix(Model$delta.1,nrow=1) %*% rho.j

    all.a.bar <- c(all.a.bar,a.h.bar)
    all.b <- rbind(all.b, b.h)
    all.a <- c(all.a, a.h)
    all.sigma <- c(all.sigma,sqrt(sigma.h.2))

    sum.rho.j_1 <- sum.rho.j_1 + rho.j # for next iteration
  }

  # Computation of forward rates:

  vec.1 <- matrix(1,dim(X)[1],1)

  aux <- vec.1 %*% matrix(all.a,nrow=1) +  X %*% t(all.b) - Model$r.bar
  aux <- aux / (vec.1 %*% matrix(all.sigma,nrow=1))

  ggg <- function(z){
    return(z * pnorm(z) + dnorm(z))
  }

  vec.f <- Model$r.bar + (vec.1 %*% matrix(all.sigma,nrow=1)) * ggg(aux)

  vec.y <- t(apply(vec.f,1,cumsum))
  vec.y <- vec.y / (vec.1 %*% matrix(1:max.H,nrow=1))

  return(list(
    all.a.bar = all.a.bar,
    all.a = all.a,
    all.b = all.b,
    all.sigma = all.sigma,
    aux = aux,
    vec.f = vec.f,
    vec.y = vec.y
  ))
}



