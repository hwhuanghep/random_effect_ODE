source("ode_formula.R")

#updating ODE variable and parameter for one individual. var.i.c and alpha.c are 
#population level values
gpode = function(tobs,yobs,x.ini,theta.ini,var.i.c,alpha.c,knots){
  n=nrow(yobs)
  k=ncol(yobs)
  q=df
  #initial values for x 
  x=x.ini
  theta=theta.ini
  bmatrix=bSpline(tobs,knots=knots,degree=3,intercept =T) 
  bmatrix1=deriv(bmatrix)
  ##initial value for all parameters
  sigma = rep(1,k)   #control goodness of fit
  gamma = rep(1,k)   #control ODE
  tau = 1
  ##obtain ds given x
  ds=bmatrix1%*%x 
  #sample sigma
  squaresum = apply((yobs-bmatrix%*%x)^2,2,sum)
  sigma[1] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum[1])
  sigma[2] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum[2])
  sigma[3] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum[3])
  sigma=rep(1e-6,3)
  #sample tau
  D=1e-3*diag(q)
  squareterm = sum(diag(t(x)%*%D%*%x))
  tau = 1/rgamma(1,shape=q*k/2+1,scale=2/squareterm)
  tau = 1
  ##sample gamma
  muode = odeformulas(x,theta,bmatrix)
  squareterm = diag((t(ds-muode)%*%(ds-muode)))
  gamma[1] = 1/rgamma(1,shape=n/2+1,scale=2/squareterm[1]) + 1e-5    
  gamma[2] = 1/rgamma(1,shape=n/2+1,scale=2/squareterm[2]) + 1e-5    
  gamma[3] = 1/rgamma(1,shape=n/2+1,scale=2/squareterm[3]) + 1e-5        
  gamma=rep(1e-6,3)
  ##sample x[,1]
  ##get results from ODE
  odeout <- odevar1(x,theta,gamma,bmatrix,bmatrix1)
  hmatrix1 = odeout$matrix
  mu1 = odeout$mu
  ##get results from penalty on dx
  hmatrix2 = D/tau
  ##get results from observations
  hmatrix3 = t(bmatrix)%*%bmatrix/sigma[1]
  mu3 = t(bmatrix)%*%yobs[,1]/sigma[1]
  sigmax = solve(hmatrix1+hmatrix2+hmatrix3)
  mux = sigmax%*%(mu1+mu3)
  x[,1] = mux+chol(sigmax)%*%rnorm(q)
  ##sample x[,2]
  ##get results from ODE
  odeout <- odevar2(x,theta,gamma,bmatrix,bmatrix1)
  hmatrix1 = odeout$matrix
  mu1 = odeout$mu
  ##get results from penalty on dx
  hmatrix2 = D/tau
  ##get results from observations      
  hmatrix3 = t(bmatrix)%*%bmatrix/sigma[2]
  mu3 = t(bmatrix)%*%yobs[,2]/sigma[2]
  sigmax = solve(hmatrix1+hmatrix2+hmatrix3)
  mux = sigmax%*%(mu1+mu3)
  x[,2] = mux+chol(sigmax)%*%rnorm(q)
  ##sample x[,3]
  ##get results from ODE
  odeout <- odevar3(x,theta,gamma,bmatrix,bmatrix1)
  hmatrix1 = odeout$matrix
  mu1 = odeout$mu
  ##get results from from penalty on dx
  hmatrix2 = D/tau
  ##get results from observations      
  hmatrix3 = t(bmatrix)%*%bmatrix/sigma[3]
  mu3 = t(bmatrix)%*%yobs[,3]/sigma[3]
  sigmax = solve(hmatrix1+hmatrix2+hmatrix3)
  mux = sigmax%*%(mu1+mu3)
  x[,3] = mux+chol(sigmax)%*%rnorm(q)
  
  ##sample theta
  thetavar = solve(var.i.c)
  thetamu = thetavar%*%(alpha.c)
  s=bmatrix%*%x
  ds=bmatrix1%*%x
  for(i in 1:n){ 
    mm = odematrix(s[i,])
    thetavar = thetavar + t(mm)%*%diag(1/gamma)%*%mm
    thetamu = thetamu + t(mm)%*%diag(1/gamma)%*%(ds[i,]-c(0,0,s[i,2]))
  }
  thetavar = solve(thetavar+diag(1e-3,3))
  thetamu = thetavar%*%thetamu
  theta = thetamu + chol(thetavar)%*%matrix(rnorm(3),3,1)
  for(i in 1:3){
    if(theta[i]<0){
      mutemp <- thetamu[i,1] + thetavar[i,-i]%*%solve(thetavar[-i,-i])%*%(theta[-i]-thetamu[-i,])
      vartemp <- thetavar[i,i] - thetavar[i,-i]%*%solve(thetavar[-i,-i])%*%thetavar[-i,i]
      theta[i] = rtnorm(1,mean=mutemp,sd=sqrt(vartemp),lower=0)
    }
  }
  return(list(x=x,theta=theta))
}
