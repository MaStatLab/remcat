remcat = function( xg.mat, xe.mat, y.vec, rho0 = 0.5, mode=1){

  is.inbound = function(x) {
    return (x == 0 || x == 1);
  }

  if (!(is.inbound(xe.mat))) {
    stop("Covarites must be binary (0 or 1)")
  }

  if (!(is.inbound(y.vec))) {
    stop("Response must be binary (0 or 1)")
  }



  ## mode specifies the
  ## mode = 1: rho = rho0
  ## mode = 2: rho = 1 - c*(1-rho0)^{level + 1}
  ## mode = 3: rho = 1 - (1-rho0) e^{-c*level}
  ## mode = 4: rho0 = rho0 e^{-c*level}

  ## current the constant c is set to 1.0

  ans = .Call('remcat_C', xg.mat, xe.mat, y.vec, ncol(xg.mat), ncol(xe.mat),rho0,mode)
  names(ans)=c("logrho","loggamma","logphi","logpsi")

  p = c(NA,NA,NA,NA)
  p[4] = exp(ans$logrho + ans$loggamma)
  p[3] = exp(ans$logpsi - ans$logphi) - p[4]
  p[1] = 1 - exp(ans$logrho) - p[3]
  p[2] = 1 - p[1] - p[3] - p[4]

  ans[[5]] = p
  names(ans)[5]=c("p")

  return(ans)
}


remcat.2toK.ind.test.bf = function(y1.mat,y2.mat,x.mat) {
  ## returns the BF for conditional independence

  ans=remcat(y1.mat,x.mat)$logphi  + remcat(y2.mat,x.mat)$logphi -
    remcat(cbind(y1.mat,y2.mat),x.mat)$logphi

  return(ans)
}
