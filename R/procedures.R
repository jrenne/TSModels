
simul.X <- function(Phi_x,Mu_x,Sigma_x,nb.periods){
  n <- dim(Phi_x)[1]
  X0 <- solve(diag(n) - Phi_x) %*% Mu_x
  X <- X0
  all.X <- c(X0)
  for(t in 2:nb.periods){
    X <- Mu_x + Phi_x %*% X + Sigma_x %*% rnorm(n)
    all.X <- rbind(all.X,t(X))
  }
  return(all.X)
}


make.Gamma <- function(Psi){
  # Psi is of dimension (H+1)x(H+1) (see paper)
  # The sum of the entries (j+1,1),(j+2,2),...,(n,n-j) of Psi is Gamma_{n,j}
  # This function computes the matrix Mat.Gamma whose component (i,j) is Gamma_{j,i-1}
  # That is, Gamma_{i+j-1,i-1} is the component (i,j) of Mat.Gamma,
  #   or Gamma_{n,k} is the component  (k+1,n-k) of Mat.Gamma

  n <- dim(Psi)[1]
  indic.entries.Psi.to.keep <- which(!upper.tri(Psi))

  AUX <- (!lower.tri(Psi))[,n:1]
  indic.entries.AUX.to.fill <- which(AUX)

  AUX <- matrix(0,n,n) #NaN
  AUX[indic.entries.AUX.to.fill] <- Psi[indic.entries.Psi.to.keep]

  Mat.Gamma <- t(apply(AUX,1,cumsum))
  return(Mat.Gamma)
}

make.indices <- function(M){
  n.m <- dim(M)[1]
  mat.indices <- matrix(NaN, (n.m-1), (n.m-1))
  for (i in 1:(n.m-1)){
    mat.indices[i:(n.m-1),i] <- 2:(n.m-i+1) + (i-1)*n.m
  }
  return(mat.indices)
}


compute.condit.Exp <- function(a,b,
                               a.dot,b.dot,
                               Phi_x,Mu_x,Sigma_x,
                               X,max.h,
                               u = 0){
  # This function computes two types of conditional expectations:
  #
  # E.n :
  # E_t(exp(b.dot + a.dot*X_{t+1} + ... + b.dot + a.dot*X_{t+h} +
  #         max(0,b + a*X_{t+1})+max(0,b + a*X_{t+1} + ... + b + a*X_{t+h})))
  # E.n.star :
  # E_t(exp(a.dot + b.dot*X_{t+1} + ... + a.dot + b.dot*X_{t+h} +
  #         max(0,b + a*X_{t+1})+max(0,b + a*X_{t+1} + ... + b + a*X_{t+h-1})))

  n_x <- dim(Phi_x)[1]
  T <- dim(X)[1]

  phi1 <- - matrix(a.dot, ncol=1)

  Phi.x.h   <- diag(n_x)
  MU        <- array(NaN, c(n_x, T, max.h))
  vec.1.T   <- matrix(1, T, 1)

  Omega <- Sigma_x%*%t(Sigma_x)

  GAMMA <- array(NaN, c(n_x,n_x, max.h)) # The i^th layer is Gamma_{i,0}. In particular, the first layer if Gamma_{1,0} = Omega

  gamma.h.0 <- matrix(0, n_x, n_x)

  KSI.a        <- matrix(NaN, n_x, max.h + 1) # the i^th column of KSI.a will be xi_{i-1}^a = t(Phi_x^i) x a
  KSI.a[,1]    <- a
  KSI.a.dot    <- matrix(NaN, n_x, max.h + 1)  # the i^th column of KSI.a.dot will be xi_{i-1}^a.dot = t(Phi_x^i) x a.dot
  KSI.a.dot[,1]<- a.dot

  SIGMA.lambda <- matrix(NaN,max.h,1)

  MU.lambda <- matrix(NaN, T, max.h)
  MU.i        <-  matrix(NaN, T, max.h)

  all.b.h <- NULL
  all.a.h <- NULL
  a.h <- matrix(0,n_x,1)
  b.h <- 0
  sum.Phi.x.h_1 <- 0

  for (h in 1:max.h){
    # This loop also computes the (a,b) loadings to compute E(exp( - Lambda (w_{t+1}+...+w_{t+h})))
    # It also computes matrices MU.lambda and MU.i

    b.h <- b.h + 1/2 * c(t( a.h + phi1 ) %*% Omega %*% ( a.h + phi1 )) +
      c(t( a.h + phi1 ) %*% Mu_x)
    a.h <- t(Phi_x) %*% ( a.h + phi1 )
    all.a.h <- cbind(all.a.h,a.h)
    all.b.h <- cbind(all.b.h,b.h)

    sum.Phi.x.h_1 <- sum.Phi.x.h_1 + Phi.x.h
    Phi.x.h <- Phi.x.h%*%Phi_x # Phi.x.h is Phi_x^h

    MU[,,h] <- (sum.Phi.x.h_1%*%Mu_x)%*%t(vec.1.T) +
      Phi.x.h%*%t(X)

    gamma.h.0 <- Omega + Phi_x%*%gamma.h.0%*%t(Phi_x)

    GAMMA[,,h] <- gamma.h.0

    SIGMA.lambda[h] <- sqrt(t(a)%*%gamma.h.0%*%a)

    MU.lambda[, h] <- b +  t(MU[,,h])%*%a
    MU.i[,h]        <- b.dot + t(MU[,,h])%*%a.dot

    KSI.a[,h+1]    <- Phi.x.h%*%a
    KSI.a.dot[,h+1]<- Phi.x.h%*%a.dot
  }

  SIGMA.lambda <- vec.1.T%*%t(SIGMA.lambda) # This is matrix of dimension T x (max.h - 1)

  AUX.lambda <-  MU.lambda/SIGMA.lambda

  P <- pnorm(AUX.lambda)
  PSI.a         <- t(KSI.a)%*%Omega%*%KSI.a
  PSI.a.dot     <- t(KSI.a.dot)%*%Omega%*%KSI.a.dot
  PSI.a.dot.a     <- t(KSI.a.dot)%*%Omega%*%KSI.a
  PSI.a.a.dot     <- t(KSI.a)%*%Omega%*%KSI.a.dot
  PSI.dot.dot     <- t(KSI.a + KSI.a.dot)%*%Omega%*%(KSI.a + KSI.a.dot)
  mat.GAMMA           <- make.Gamma(PSI.a)
  mat.GAMMA.dot       <- make.Gamma(PSI.a.dot)
  mat.GAMMA.left.dot  <- make.Gamma(PSI.a.dot.a)
  mat.GAMMA.right.dot <- make.Gamma(PSI.a.a.dot)
  mat.GAMMA.dot.dot   <- make.Gamma(PSI.dot.dot)

  MU.0 <- t(X)
  MU.i.lagged <- cbind(b.dot + t(MU.0)%*%a.dot,MU.i[,1:(max.h-1)])

  MU.i.lagged        <- MU.i

  F.n_1.n      <- matrix(NaN,T,max.h) # The first column is for (n-1,n) = (0,1)
  F.n_1.n.star <- matrix(NaN,T,max.h) # The first column is for (n-1,n) = (0,1)


  # Treat F.n_1.n: ------------------
  F.n_1.n        <- MU.i.lagged        + P*MU.lambda + dnorm(-AUX.lambda)*SIGMA.lambda

  F.n_1.n[,2:max.h] <- F.n_1.n[,2:max.h] -
    1/2*(P[,2:(dim(P)[2])]*(vec.1.T%*%matrix(mat.GAMMA.dot.dot[1,2:(dim(P)[2])],nrow=1))+
           (1-P[,2:(dim(P)[2])])*(vec.1.T%*%matrix(mat.GAMMA.dot[1,2:(dim(P)[2])],nrow=1)))


  # Treat F.n_1.n.star: ------------------
  F.n_1.n.star <- MU.i.lagged

  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] + P[,1:(dim(P)[2]-1)]*MU.lambda[,1:(dim(P)[2]-1)] +
    dnorm(-AUX.lambda[,1:(dim(P)[2]-1)])*SIGMA.lambda[,1:(dim(P)[2]-1)]

  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] -
    1/2*(P[,2:dim(P)[2]]*(vec.1.T%*%matrix(mat.GAMMA.dot[1,2:(dim(P)[2])],nrow=1))+
           (1-P[,2:dim(P)[2]])*(vec.1.T%*%matrix(mat.GAMMA.dot[1,2:(dim(P)[2])],nrow=1)))

  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] -
    P[,2:dim(P)[2]]*(vec.1.T%*%matrix(mat.GAMMA.right.dot[2,2:(dim(P)[2])],nrow=1))

  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] -
    1/2*(P[,2:dim(P)[2]]*(vec.1.T%*%matrix(mat.GAMMA[1,2:dim(P)[2]],nrow=1)))

  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] -
    vec.1.T%*%matrix(mat.GAMMA.dot[(dim(mat.GAMMA.dot)[1]),2:(dim(P)[2])],nrow=1) +
    P[,1]*( vec.1.T%*%matrix(mat.GAMMA.right.dot[(dim(mat.GAMMA.right.dot)[1]-1),1:(dim(P)[2]-1)],nrow=1)
    )


  SUM.GAMMA.P.n_1.n      <- matrix(0, T, max.h)
  SUM.GAMMA.P.n_1.n.star <- matrix(0, T, max.h)

  mat.indices <- make.indices(PSI.a.dot)

  for (n in 3:max.h){

    AUX.SUM.GAMMA.P.n_1.n <- P[,1:(n-2)]*vec.1.T%*%matrix(mat.GAMMA.dot.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1) +
      (1-P[,1:(n-2)])*(vec.1.T%*%matrix(mat.GAMMA.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1))

    AUX.SUM.GAMMA.P.n_1.n.star <- P[,1:(n-2)]*
      (
        vec.1.T%*%matrix(  mat.GAMMA.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1)+
          vec.1.T%*%matrix(mat.GAMMA.right.dot[mat.indices[(n-2),1:(n-2)]],nrow=1) +
          vec.1.T%*%matrix(mat.GAMMA.left.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1)+
          vec.1.T%*%matrix(mat.GAMMA[mat.indices[(n-2),1:(n-2)]],nrow=1)
      ) +
      (1-P[,1:(n-2)])*(vec.1.T%*%matrix(mat.GAMMA.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1))

    SUM.GAMMA.P.n_1.n[,n]      <- AUX.SUM.GAMMA.P.n_1.n%*%matrix(1,n-2,1)
    SUM.GAMMA.P.n_1.n.star[,n] <- AUX.SUM.GAMMA.P.n_1.n.star%*%matrix(1,n-2,1)
  }

  F.n_1.n        <- F.n_1.n        - SUM.GAMMA.P.n_1.n
  F.n_1.n.star   <- F.n_1.n.star   - SUM.GAMMA.P.n_1.n.star

  cumsum.f.n_1.n        <- t(apply(F.n_1.n, 1, cumsum))
  cumsum.f.n_1.n.star <- t(apply(F.n_1.n.star, 1, cumsum))
  E.n        <- exp(-cumsum.f.n_1.n)
  E.n.star   <- exp(-cumsum.f.n_1.n.star)

  return(list(E.n = E.n,
              E.n.star = E.n.star,
              all.a.h = all.a.h,
              all.b.h = all.b.h))
}




