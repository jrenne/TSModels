compute_lambdas <- function(AB,model){
  # A is of dimension N * 1
  # B is of dimension N * m

  m <- dim(model$Phi)[1]   # dimension of state vector
  N <- length(model$Alpha) # maximum maturity

  A      <- AB$A
  B      <- AB$B

  # Load relevant objects from model:
  gamma = model$gamma
  Sigma = model$Sigma
  Alpha = model$Alpha
  Beta  = model$Beta
  Zeta  = model$Zeta
  a1    = model$a1
  b1    = model$b1

  Theta      <- rbind(0,B[1:(N-1),])

  caligA <- Alpha + Zeta/(1:N) * A
  lambda <- gamma * t(Sigma) %*% t(Theta) %*% caligA

  caligB <- Beta + ((Zeta/(1:N)) %*% matrix(1,1,m)) * B
  Lambda <- gamma * t(Sigma) %*% t(Theta) %*% caligB

  # Compute alpha1 and beta1 needed to have sum(z)=1:
  alpha1 <- 1 - sum(caligA[2:N]) - Zeta[1]*a1
  beta1  <- - apply(caligB[2:N,],2,sum) - Zeta[1]*b1

  # Update caligA and caligB:
  Alpha[1] <- alpha1
  Beta[1,] <- beta1
  caligA <- Alpha + Zeta/(1:N) * A
  caligB <- Beta + ((Zeta/(1:N)) %*% matrix(1,1,m)) * B

  return(list(lambda = lambda,
              Lambda = Lambda,
              alpha1 = alpha1,
              beta1  = beta1,
              Alpha  = Alpha,
              Beta   = Beta,
              caligA = caligA,
              caligB = caligB))
}

compute_AB_VV <- function(model,lambdas){
  lambda <- lambdas$lambda
  Lambda <- lambdas$Lambda

  # Load relevant objects from model:
  mu_f  = model$mu_f
  Phi   = model$Phi
  Sigma = model$Sigma
  a1    = model$a1
  b1    = model$b1

  m <- dim(Phi)[2]

  Phi_Q  <- Phi  - Sigma %*% Lambda
  mu_f_Q <- mu_f - Sigma %*% lambda

  a_n_1 <- 0
  b_n_1 <- matrix(0,m,1)

  A      <- matrix(NaN,N,1)
  B      <- matrix(NaN,N,m)

  for(n in 1:N){
    a_n <- a1 + a_n_1 + t(b_n_1) %*% mu_f_Q +
      .5 * t(b_n_1) %*% Sigma %*% t(Sigma) %*% b_n_1
    b_n <- b1 + t(Phi_Q) %*% b_n_1

    a_n_1 <- a_n
    b_n_1 <- b_n

    A[n]  <- a_n
    B[n,] <- b_n
  }

  a <- - A/1:N
  b <- - B/matrix(1:N,N,m)

  return(list(A=A,B=B,a=a,b=b))
}


solve_model <- function(model,
                        ini.lambdas = list(lambda=NaN,Lambda=NaN),
                        nb.iter = 20,
                        indic.print = FALSE){

  # Solve for lambdas:
  if(is.na(ini.lambdas$lambda[1])){
    N <- dim(model$Beta)[1]
    m <- dim(model$Beta)[2]
    lambdas <- list(lambda = matrix(0,m,1),
                    Lambda = matrix(0,m,m))
  }
  for(i in 1:nb.iter){
    AB <- compute_AB_VV(model,lambdas)
    lambdas <- compute_lambdas(AB,model)
    if(indic.print){
      print(lambdas$lambda)
    }
  }

  return(list(
    AB = AB,
    lambdas = lambdas
  ))
}

