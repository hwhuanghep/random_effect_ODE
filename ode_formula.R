#obtain the derivative vector based on given the spline coefficient x, ODE parameter x, 
# and spline matrix bmatrix  
odeformulas = function(x,theta,bmatrix){
  beta = theta[1];
  delta = theta[2];
  cw = theta[3];
  y=bmatrix%*%x
  du = -beta*y[,1]*y[,3];
  di = beta*y[,1]*y[,3]-delta*y[,2];
  dv = y[,2]-cw*y[,3];
  dudidv = cbind(du,di,dv)
  return(dudidv)
}

#matrix used in sampling the coefficient vector based on ODE observation
odematrix = function(y){
  m = matrix(0,3,3)
  m[1,] = c(-y[1]*y[3],0,0)
  m[2,] = c(y[1]*y[3],-y[2],0)
  m[3,] = c(0,0,-y[3]) 
  return(m)
}

#ODE contribution to the sampling of the first ODE variable. x-spline coefficient, 
#theta-ODE parameter, gamma-variance for ODE discrepancy, bmatrix-spline matrix
#for ODE variable, bmatrix1-spline matrix for derivative of ODE variable 
odevar1 <- function(x,theta,gamma,bmatrix,bmatrix1){
  beta = theta[1];
  delta = theta[2];
  c = theta[3];
  s=bmatrix%*%x
  ds=bmatrix1%*%x
  hmatrix <- t(bmatrix1+beta*bmatrix*s[,3])%*%(bmatrix1+beta*bmatrix*s[,3])/gamma[1]
  hmatrix=hmatrix+t(beta*bmatrix*s[,3])%*%(beta*bmatrix*s[,3])/gamma[2]
  mu <- t(beta*bmatrix*s[,3])%*%(ds[,2]+delta*s[,2])/gamma[2]
  return(list(matrix=hmatrix,mu=mu))
}

#ODE contribution to the sampling of the second ODE variable 
odevar2 = function(x,theta,gamma,bmatrix,bmatrix1){
  beta = theta[1];
  delta = theta[2];
  c = theta[3];
  s=bmatrix%*%x
  ds=bmatrix1%*%x
  hmatrix <- t(bmatrix1+delta*bmatrix)%*%(bmatrix1+delta*bmatrix)/gamma[2]
  hmatrix=hmatrix+t(bmatrix)%*%bmatrix/gamma[3]
  mu <- t(bmatrix1+delta*bmatrix)%*%(beta*s[,1]*s[,3])/gamma[2]
  mu=mu+t(bmatrix)%*%(ds[,3]+c*s[,3])/gamma[3]
  return(list(matrix=hmatrix,mu=mu))
}

#ODE contribution to the sampling of the third ODE variable 
odevar3 <- function(x,theta,gamma,bmatrix,bmatrix1){
  q=ncol(bmatrix)
  beta = theta[1];
  delta = theta[2];
  c = theta[3];
  s=bmatrix%*%x
  ds=bmatrix1%*%x
  hmatrix <- beta^2*t(s[,1]*bmatrix)%*%(s[,1]*bmatrix)*(1/gamma[1]+1/gamma[2])
  hmatrix=hmatrix+t(bmatrix1+c*bmatrix)%*%(bmatrix1+c*bmatrix)/gamma[3]
  mu <- -beta*t(s[,1]*bmatrix)%*%ds[,1]/gamma[1]
  mu=mu+beta*t(s[,1]*bmatrix)%*%(ds[,2]+delta*s[,2])/gamma[2]
  mu=mu+t(bmatrix1+c*bmatrix)%*%(s[,2])/gamma[3]
  return(list(matrix=hmatrix,mu=mu))
}
