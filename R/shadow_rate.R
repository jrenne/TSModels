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

ggg <- function(z){
  return(z * pnorm(z) + dnorm(z))
}





Kalman_filter <- function(Y_t,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos=0,
                          Rfunction=Rf, Qfunction=Qf){
  # Measurement equations:
  #    y_t   = mu_t + G * rho_t + M * eps_t
  # Transition equations:
  #    rho_t = nu_t + H * rho_t-1 + N * xi_t
  # Y_t is a T*ny-dimensional vector
  # nu_t and mu_t depend on past only (predetermined at date t)
  #
  # The function returns:
  # .$r   filtered variables (rho_t|t) size T*nr
  # .$S  SIGMA t|t T*(nr*nr)

  # Number of observed variables:
  ny = NCOL(Y_t)
  # Number of unobs. variables:
  nr = NCOL(G)
  # Number of time periods:
  T = NROW(Y_t)

  # loglik.vector will contain the vector of date-specific log-likelihood
  loglik.vector <- NULL

  #Initilize output matrices:
  rho_tt    = matrix(0,T,nr)
  rho_tp1_t = matrix(0,T,nr)
  y_tp1_t   = matrix(0,T,ny)

  Sigma_tt    = matrix(0,T,nr*nr)
  Sigma_tp1_t = matrix(0,T,nr*nr)
  Omega_tt    = matrix(0,T,ny*ny)
  Omega_tp1_t = matrix(0,T,ny*ny)

  #Initilize log-Likelihood:
  logl = ny*T/2*log(2*pi)


  #Kalman algorithm:

  for (t in 1:T){

    # print(t)

    # ==========================================================================
    # Forecasting step (between t-1 and t):
    # ==========================================================================
    if(t==1){
      rho_tp1_t[1,] = nu_t[1,] + t(H %*% rho_0)
      R = Rfunction(M,rho_0) #Rfunction defined below (measurement eq.)
      Q = Qfunction(N,rho_0) #Qfunction defined below (transition eq.)
    } else {
      rho_tp1_t[t,] = nu_t[t,] + t(H %*% rho_tt[t-1,])
      R = Rfunction(M,rho_tt[t-1,],t)
      Q = Qfunction(N,rho_tt[t-1,],t)
    }
    if(sum(indic_pos==1)>0){
      rho_tp1_t[t,indic_pos==1] = pmax(rho_tp1_t[t,indic_pos==1],0)
    }
    y_tp1_t[t,] = mu_t[t,] + t(G %*% rho_tp1_t[t,])

    if(t==1){
      aux_Sigma_tp1_t = Q + H %*% Sigma_0 %*% t(H)
    } else {
      aux_Sigma_tp1_t = Q + H %*% matrix(Sigma_tt[t-1,],nr,nr) %*% t(H)
    }

    Sigma_tp1_t[t,] = matrix( aux_Sigma_tp1_t ,1,nr*nr)
    omega           = R + G %*% aux_Sigma_tp1_t %*% t(G)
    Omega_tp1_t[t,] = matrix(omega,1,ny*ny)


    # ==========================================================================
    # Updating step
    # ==========================================================================

    # Detect observed variables:
    vec.obs.indices <- which(!is.na(Y_t[t,])) # indices of observed variables
    ny.aux <- length(c(vec.obs.indices)) # number of observed variables
    # Resize matrices accordingly:
    G.aux <- matrix(G[vec.obs.indices,],nrow=ny.aux)
    R.aux <- R[vec.obs.indices,]
    omega <- omega[vec.obs.indices,]
    if(class(R.aux)[1]=="numeric"){
      R.aux <- R.aux[vec.obs.indices]
      omega <- omega[vec.obs.indices]
    }else{
      R.aux <- R.aux[,vec.obs.indices]
      omega <- omega[,vec.obs.indices]
    }
    R.aux <- matrix(R.aux,ny.aux,ny.aux)

    #Compute gain K:
    #print(G.aux)
    if(dim(G.aux)[1]>0){
      if(sum(eigen(R.aux + G.aux %*% aux_Sigma_tp1_t %*% t(G.aux))$values<=0)>0){
        print("ici")
      }
      K = aux_Sigma_tp1_t %*% t(G.aux) %*% ginv(R.aux + G.aux %*% aux_Sigma_tp1_t %*% t(G.aux))


      lambda_t   = Y_t[t,] - y_tp1_t[t,]
      lambda_t <- matrix(lambda_t[vec.obs.indices],ncol=1)
      rho_tt[t,] = t( rho_tp1_t[t,] + K %*% lambda_t )
      if(sum(indic_pos==1)>0){
        rho_tt[t,indic_pos==1] = pmax(rho_tt[t,indic_pos==1],0)
      }
      Id         = diag(1,nrow=nr,ncol=nr)
      Sigma_tt[t,] = matrix( (Id - K %*% G.aux) %*% aux_Sigma_tp1_t,1,nr*nr)

      #calcul de la log-vraisemblance
      #print(det(omega))
      if(length(c(omega))==1){
        det.omega <- omega
      }else{
        det.omega <- det(omega)
      }
      loglik.vector <- rbind(loglik.vector,
                             ny.aux/2*log(2*pi) - 1/2*(log(det.omega) +
                                                         t(lambda_t) %*% ginv(omega) %*% lambda_t)
      )
      logl <- logl + loglik.vector[t]

    }else{
      rho_tt[t,] = t( rho_tp1_t[t,])
      Sigma_tt[t,] = matrix( aux_Sigma_tp1_t,1,nr*nr)
    }

  }

  fitted.obs <- mu_t + rho_tt %*% t(G) # fitted observables

  output = list(r=rho_tt,Sigma_tt=Sigma_tt,loglik=logl,y_tp1_t=y_tp1_t,
                S_tp1_t=Sigma_tp1_t,r_tp1_t=rho_tp1_t,
                loglik.vector = loglik.vector,Omega_tp1_t=Omega_tp1_t,M=M,
                fitted.obs=fitted.obs)
  return(output)
}

Rf <- function(M,RHO,t=0){
  return(M %*% t(M))
}
Qf <- function(N,RHO,t=0){
  return(N %*% t(N))
}

Kalman_smoother <- function(Y_t,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos=0,
                            Rfunction=Rf, Qfunction=Qf){

  output_filter <- Kalman_filter(Y_t,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos,
                                 Rfunction,Qfunction)
  rho_tt        <- output_filter$r
  Sigma_tt      <- output_filter$Sigma_tt
  Sigma_tp1_t   <- output_filter$S_tp1_t
  rho_tp1_t     <- output_filter$r_tp1_t
  Omega_tp1_t   <- output_filter$Omega_tp1_t

  # Number of observed variables:
  ny = NCOL(Y_t)
  # Number of unobs. variables:
  nr = NCOL(G)
  # Number of time periods:
  TT = NROW(Y_t)

  if(class(M)[1]=="list"){
    M.aux <- M$sigmas
  }else{
    M.aux <- M
  }

  #R = M%*%t(M)
  #Q = N%*%t(N)

  #Initilize output matrices:
  rho_tT       <- matrix(0,TT,nr)
  rho_tT[TT,]   <- rho_tt[TT,]
  Sigma_tT     <- matrix(0,TT,nr*nr)
  Omega_tT     <- matrix(0,TT,ny*ny) # This is the covariance matrix of Y_t if Y_t is missing
  Sigma_tT[TT,] <- Sigma_tt[TT,]
  Sigma_tT.aux <- matrix(Sigma_tt[TT,],nr,nr)
  R = Rfunction(M,rho_0) #Rfunction defined below (measurement eq.)
  Omega_tT[TT,] <- c(
    R + G %*% Sigma_tT.aux %*% t(G)
  )
  Var.model.implied.obs <- NaN * Y_t # this matrix will contain the variances of G*rho_t
  Var.model.implied.obs[TT,] <- diag(
    G %*% matrix(Sigma_tT[TT,],nr,nr) %*% t(G))

  for(t in seq(TT-1,1,by=-1)){
    F <- matrix(Sigma_tt[t,],nr,nr) %*% t(H) %*% solve(matrix(Sigma_tp1_t[t+1,],nr,nr))
    rho_tT[t,] <- rho_tt[t,] + t( F %*% (rho_tT[t+1,]-rho_tp1_t[t+1,]) )
    if(sum(indic_pos==1)>0){
      rho_tT[t,indic_pos==1] = pmax(rho_tT[t,indic_pos==1],0)
    }
    Sigma_tT.aux <- Sigma_tt[t,] + F %*%
      matrix(Sigma_tT[t+1,]-Sigma_tp1_t[t+1,],nr,nr) %*% t(F)
    Sigma_tT[t,] <- c(Sigma_tT.aux)
    R = Rfunction(M,rho_0) #Rfunction defined below (measurement eq.)
    Omega_tT[t,] <- c(
      R + G %*% Sigma_tT.aux %*% t(G)
    )
  }
  fitted.obs <- mu_t + rho_tT %*% t(G) # fitted observables
  output = list(r_smooth=rho_tT,S_smooth=Sigma_tT,
                loglik=output_filter$loglik,
                loglik.vector=output_filter$loglik.vector,
                Omega_tp1_t=Omega_tp1_t,
                #Var.model.implied.obs=Var.model.implied.obs,
                M=M,
                Sigma_tt=Sigma_tt,
                fitted.obs=fitted.obs,G=G,mu_t=mu_t,
                Omega_tT = Omega_tT)
  return(output)
}


