##this code give direct ode parameter estimation method for flu virus model
##based on spline Bayesian MCMC method including
##derivative term and random effect variable selection
rm(list=ls())
graphics.off(); #close all graphics windows
require(deSolve)  #loads ODE solver package
require(mvtnorm)
require(nloptr) #for fitting
require(GenSA)
require(msm)
require(splines2)
require(LaplacesDemon)

rk<-function(x,z){ 
  x=x/scale
  z=z/scale 
  ((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4-((abs(x-z)-0.5)^4-(abs(x-z)-0.5)^2/2+7/240)/24
}

rk1<-function(x,z){ 
  x=x/scale
  z=z/scale 
  out=((z-0.5)^2-1/12)*(x-0.5)/2-(4*(abs(x-z)-0.5)^3-(abs(x-z)-0.5))*sign(x-z)/24
  out/scale
}

# set up the penalized regression spline penalty matrix, given knot sequence xk
spl.S<-function(xk){ 
  q<-length(xk)+2;
  S<-matrix(0,q,q) # initialize matrix to 0
  S[3:q,3:q]<-outer(xk,xk,FUN=rk) # fill in non-zero part
  S
}

# set up the design matrix, given knot sequence xk
spl.X<-function(x,xk){
  n=length(x) 
  q<-length(xk)+2;
  S=matrix(1,n,q)
  S[,2]=x/scale
  S[,3:q]<-outer(x,xk,FUN=rk) 
  S
}

spl.X1<-function(x,xk){
  n=length(x) 
  q<-length(xk)+2;
  S=matrix(0,n,q)
  S[,2]=1/scale
  S[,3:q]<-outer(x,xk,FUN=rk1) 
  S
}

odeequations = function(t,y,x){
  kapa = x[1];
  delta = x[2];
  cw = x[3];
  du = -kapa*y[1]*y[3];
  di = kapa*y[1]*y[3]-delta*y[2];
  dv = y[2]-cw*y[3];
  dudidv = c(du,di,dv)
  return(list(dudidv))
}

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

odematrix = function(y){
  m = matrix(0,3,3)
  m[1,] = c(-y[1]*y[3],0,0)
  m[2,] = c(y[1]*y[3],-y[2],0)
  m[3,] = c(0,0,-y[3]) 
  return(m)
}

odevar1 <- function(x,theta,gamma,bmatrix,bmatrix1){
  q=ncol(bmatrix)
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

odevar2 = function(x,theta,gamma,bmatrix,bmatrix1){
  q=ncol(bmatrix)
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

fitfunction=function(pars.all,virusdata){
  atolv=1e-7; 
  rtolv=1e-7;
  Y0=c(pars.all[1],0,max(0.1,virusdata[1,2]))
  pars.ode=pars.all[-1]
  tobs=virusdata[,1]
  odestack = try(lsoda(Y0,tobs,odeequations,pars.ode,atol=atolv,rtol=rtolv));
  if(length(odestack)==1) {
    cat('!!!!!!!!!!!!!unresolvable integrator error early return\n')
    return(1e10);
  } 
  vir.SSR=sum((virusdata[,2]-odestack[,4])^2)
  Fobject = vir.SSR; #the full objective function to be returned        
  
  if(is.na(Fobject)){
    Fobject=1e10
  }
  return(Fobject)
}

gpode = function(tobs,yobs,x.ini,theta.ini,var.i.c,alpha.c,
                 missing,conditioning,knots){
  n=nrow(yobs)
  k=ncol(yobs)
  q=df
  #initial values for x 
  x=x.ini
  theta=theta.ini
  bmatrix=bSpline(tobs,knots=knots,degree=3,intercept =T) 
  bmatrix1=deriv(bmatrix)
  #make orthogonal matrix from bmatrix
  mothg=diag(df)
  a0=sqrt(sum(bmatrix[1,]^2))
  mothg[,1]=bmatrix[1,]/a0
  tmatrix=matrix(rnorm(q*q),q,q)
  for(i in 2:q){
    mothg[,i]=tmatrix[,i]
    for(j in 1:(i-1)){
      mothg[,i]=mothg[,i]-sum(mothg[,i]*mothg[,j])*mothg[,j] 
    }
    mothg[,i]=mothg[,i]/sqrt(sum(mothg[,i]^2))
  }
  ##initial value for all parameters
  sigma = rep(1,k)   #control goodness of fit
  gamma = rep(1,k)   #control ODE
  tau = 1
  ##obtain ds given x
  ds=bmatrix1%*%x 
  #sample sigma
  squaresum = apply((yobs-bmatrix%*%x)^2,2,sum)
  if(!missing[1]){
    sigma[1] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum[1])
  }
  if(!missing[2]){
    sigma[2] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum[2])
  }
  if(!missing[3]){
    sigma[3] = 1/rgamma(1,shape=n/2+1,rate=2/squaresum[3])
  }
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
  if(missing[1]|missing[2]){
    gamma=rep(1e-6,3)
  }
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
  if(missing[1]){
    sigmax = solve(hmatrix1+hmatrix2)
    mux = sigmax%*%mu1
  }else{
    sigmax = solve(hmatrix1+hmatrix2+hmatrix3)
    mux = sigmax%*%(mu1+mu3)
  }
  x[,1] = mux+chol(sigmax)%*%rnorm(q)
  if(missing[1]&conditioning[1]){
    mux=t(mothg)%*%mux
    sigmax=t(mothg)%*%sigmax%*%mothg
    mux1=mux[2:q]+sigmax[2:q,1]/sigmax[1,1]*(u0/a0-mux[1])
    sigma12=matrix(sigmax[2:q,1],ncol=1)
    sigmax1=sigmax[2:q,2:q]-sigma12%*%t(sigma12)/sigmax[1,1]
    xt=rep(0,q)
    xt[1]=u0/a0
    xt[2:q] = mux1+chol(sigmax1)%*%rnorm(q-1)
    x[,1]=mothg%*%xt      
  }    
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
  if(missing[2]){
    sigmax = solve(hmatrix1+hmatrix2)
    mux = sigmax%*%mu1
  }else{
    sigmax = solve(hmatrix1+hmatrix2+hmatrix3)
    mux = sigmax%*%(mu1+mu3)
  }
  x[,2] = mux+chol(sigmax)%*%rnorm(q)
  if(missing[2]&conditioning[2]){
    mux=t(mothg)%*%mux
    sigmax=t(mothg)%*%sigmax%*%mothg
    mux1=mux[2:q]+sigmax[2:q,1]/sigmax[1,1]*(i0/a0-mux[1])
    sigma12=matrix(sigmax[2:q,1],ncol=1)
    sigmax1=sigmax[2:q,2:q]-sigma12%*%t(sigma12)/sigmax[1,1]
    xt=rep(0,q)
    xt[1]=i0/a0
    xt[2:q] = mux1+chol(sigmax1)%*%rnorm(q-1)      
    x[,2]=mothg%*%xt      
  }    
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

setwd("E:/paper/duck_paper_2")
atolv=1e-7; 
rtolv=1e-7;
model.d=3

#model setting
missing=c(0,0,0)
conditioning=c(0,0,0)

##hyper parameters
hyp.a = 1e-5
hyp.b = 1e-5
hyp.nu = 1e-5
hyp.eta = rep(0,model.d)
hyp.lambda = diag(1e5,model.d);
hyp.omega = diag(1e5,model.d)

#mcmc procedure setting
icount = 0
nrep=25000
nburn=5000
nskip=10

#read data and extract information for strains
rawdata=read.csv("camille-duck-data.csv",header=TRUE)
#number of subjects in each cluster
cts=as.character(unique(rawdata$Strain))
#total number of clusters
nct=length(cts)
ncts=rep(0,nct)
for(ctl in 1:nct){
  ncts[ctl]=length(unique(rawdata[rawdata$Strain==cts[ctl],"DuckID"]))
}
bn=sum(ncts)
#impute the observed data to get the completed observation

data.c = list()
y0.c = list()
x.c = list()
theta.c = list()

tvec=seq(1,14,by=0.3)
knots=c(seq(1.1,2.9,by=0.2),seq(3.1,12,by=2))
bmatrix=bSpline(tvec,knots=knots,degree=3,intercept =T) 
df=ncol(bmatrix)

theta.s = list()
duck.id=0
for(ctl in 1:nct){
  y0.c[[ctl]]=rep(0,3)
  data.c[[ctl]] = list()
  x.c[[ctl]] = array(0,c(ncts[ctl],df,3))
  theta.c[[ctl]] = matrix(0,ncts[ctl],model.d)
  theta.s[[ctl]] = matrix(0,ncts[ctl],model.d)
  for(i in 1:ncts[ctl]){
    duck.id=duck.id+1
    tobs0 = rawdata[rawdata$Strain==cts[ctl]&rawdata$DuckID==duck.id,"DayPI"]
    yobs0 = as.vector(rawdata[rawdata$Strain==cts[ctl]&rawdata$DuckID==duck.id,"CloacaVirus"])
    tobs=tobs0[!is.na(yobs0)]
    yobs=yobs0[!is.na(yobs0)]
    virusdata=cbind(tobs,yobs)
    p.0=c(6,rep(0.5,3))
    plevel = 0
    max.stps = 1; 
    max.time = 5*24*60*60 
    xerr.global = 1e-7;
    xerr.local = 1e-6;   #used for multi-start solver MLS
    xerr.general = 1e-7; #used for all other solvers
    fres <- nloptr(x0=p.0,eval_f=fitfunction,
                   opts=list("algorithm"="NLOPT_LN_BOBYQA",
                             xtol_rel=xerr.global,maxtime=max.time,maxeval=-1,
                             print_level=plevel,
                             "local_opts"=list("algorithm"="NLOPT_LN_COBYLA",
                                               xtol_rel=xerr.local,maxeval=max.stps)),
                   virusdata=virusdata)  
    Y0=c(fres$solution[1],0,max(0.1,yobs[1]))
    pars.ode=fres$solution[-1]
    odestack = try(lsoda(Y0,tvec,odeequations,pars.ode,atol=atolv,rtol=rtolv))
    y0.c[[ctl]]=rbind(y0.c[[ctl]],Y0)
    data.c[[ctl]][[i]]=odestack
    x.c[[ctl]][i,,]=solve(t(bmatrix)%*%bmatrix+1e-7*diag(df))%*%t(bmatrix)%*%odestack[,c(2:4)]
    
    theta.c[[ctl]][i,]=pars.ode
  }
  y0.c[[ctl]]=y0.c[[ctl]][-1,]
}
#save x.c and theta.c into x.c.0 and theta.c.0
x.c.0=x.c
theta.c.0=theta.c
x.c=x.c.0
theta.c=theta.c.0
##initial values for x and theta
alpha.c <- matrix(0,nct,model.d)
for(ctl in 1:nct){
  alpha.c[ctl,]=apply(theta.c[[ctl]],2,sum)
}
mu.c=apply(alpha.c,2,sum)/bn
for(ctl in 1:nct){
  alpha.c[ctl,]=alpha.c[ctl,]-mu.c
}
p0=0.5
mu.0=0
sigma.0=100
lambda.c=rep(0,model.d)
for(d in 1:model.d){
  if(runif(1)<p0){
    lambda.c[d]=0
  }else{
    lambda.c[d]=rtnorm(1,mean=0,sd=sqrt(sigma.0),lower=0)
  }
}
var.i.c = diag(rep(1e-3,model.d))
gamma.c=rnorm(model.d,0,sqrt(0.5))
var.s.c.r=rbind(c(1,0,0),c(gamma.c[1],1,0),c(gamma.c[2],gamma.c[3],1))
var.s.c.l=diag(lambda.c)%*%var.s.c.r
var.s.c=var.s.c.l%*%t(var.s.c.l)

x.s = list()
for(ctl in 1:nct){
  x.s[[ctl]]=matrix(0,df,3)
}
mu.s=NULL
alpha.s <- matrix(0,nct,model.d)
alpha.s.s <- matrix(0,nct,model.d)
index.s=rep(0,8)
for(ii in 1:nrep){  
  #update x and theta
  for(ctl in 1:nct){
    for(i in 1:ncts[ctl]){
      odestack=data.c[[ctl]][[i]]
      tobs=odestack[,1]
      yobs=odestack[,-1]
      x.ini = x.c[[ctl]][i,,]
      theta.ini = theta.c[[ctl]][i,]
      odeout = gpode(tobs,yobs,x.ini,theta.ini,var.i.c,alpha.c[ctl,]+mu.c,
                     missing,conditioning,knots)
      x.c[[ctl]][i,,] = odeout$x
      theta.c[[ctl]][i,] = odeout$theta
    }
  }
  ##update alpha.c and b.c
  b.c=NULL
  for(ctl in 1:nct){
    sum.theta <- apply(theta.c[[ctl]],2,sum)-ncts[ctl]*mu.c
    norm.var <- solve(ncts[ctl]*t(var.s.c.l)%*%solve(var.i.c)%*%var.s.c.l+diag(model.d))
    norm.mean <- norm.var%*%t(var.s.c.l)%*%solve(var.i.c)%*%sum.theta
    temp=as.vector(norm.mean+chol(norm.var)%*%matrix(rnorm(model.d),model.d,1))
    b.c=rbind(b.c,temp)
    alpha.c[ctl,]=var.s.c.l%*%temp
  }  
  
  ##update mu.c
  sum.theta=rep(0,model.d)
  for(ctl in 1:nct){
    sum.theta <- sum.theta+apply(theta.c[[ctl]],2,sum)-ncts[ctl]*alpha.c[ctl,]
  }  
  norm.var = solve(bn*solve(var.i.c)+solve(hyp.lambda))
  norm.mean = norm.var%*%(solve(hyp.lambda)%*%hyp.eta+solve(var.i.c)%*%sum.theta)
  mu.c = as.vector(norm.mean+chol(norm.var)%*%matrix(rnorm(model.d),model.d,1))
  for(d in 1:model.d){
    if(mu.c[d]<0){
      mu.temp <- norm.mean[d]+norm.var[d,-d]%*%solve(norm.var[-d,-d])%*%(mu.c[-d]-norm.mean[-d])
      var.temp <- norm.var[d,d]-norm.var[d,-d]%*%solve(norm.var[-d,-d])%*%norm.var[-d,d]
      mu.c[d]=rtnorm(1,mean=mu.temp,sd=sqrt(var.temp),lower=0)
    }
  }
  
  ##update var.i 
  sum.theta = matrix(0,model.d,model.d)
  for(ctl in 1:nct){
    ttvec <- t(theta.c[[ctl]]) - alpha.c[ctl,] - mu.c
    sum.theta <- sum.theta + ttvec%*%t(ttvec)
  }
  wishart.nu = bn + hyp.nu
  wishart.sigma = sum.theta + solve(hyp.omega)
  var.i.c =rinvwishart(wishart.nu, wishart.sigma)
  
  ##update var.s
  
  #first update lambda
  var.lambda.in=diag(0,model.d)
  sum.lambda=rep(0,model.d)
  for(ctl in 1:nct){
    b.m=diag(as.vector(var.s.c.r%*%b.c[ctl,]))
    var.lambda.in=var.lambda.in+ncts[ctl]*b.m%*%solve(var.i.c)%*%b.m
    sum.lambda <- sum.lambda + b.m%*%solve(var.i.c)%*%(apply(theta.c[[ctl]],2,sum)-ncts[ctl]*mu.c)
  }
  var.lambda=solve(var.lambda.in)
  mean.lambda=var.lambda%*%sum.lambda
  #sample from truncated normal
  for(d in 1:model.d){
    mu.temp <- mean.lambda[d]+var.lambda[d,-d]%*%solve(var.lambda[-d,-d])%*%(lambda.c[-d]-mean.lambda[-d])
    var.temp <- var.lambda[d,d]-var.lambda[d,-d]%*%solve(var.lambda[-d,-d])%*%var.lambda[-d,d]
    sigma.bar=var.temp*sigma.0/(var.temp+sigma.0)
    mu.bar=sigma.bar*(mu.temp/var.temp+mu.0/sigma.0)
    fact=sqrt(sigma.bar/sigma.0)*pnorm(mu.bar/sqrt(sigma.bar))/pnorm(mu.0/sqrt(sigma.0))
    p0.bar=p0/(p0+(1-p0)*fact*exp(mu.bar^2/sigma.bar/2-mu.0^2/sigma.0/2))
    if(runif(1)<p0.bar){
      lambda.c[d]=0
    }else{
      lambda.c[d]=rtnorm(1,mean=mu.bar,sd=sqrt(sigma.bar),lower=0)
    }
  }
  #second update gamma
  var.gamma.in=diag(0,model.d)
  sum.gamma=rep(0,model.d)
  for(ctl in 1:nct){
    a.m=cbind(c(0,lambda.c[2]*b.c[ctl,1],0),c(0,0,lambda.c[3]*b.c[ctl,1]),c(0,0,lambda.c[3]*b.c[ctl,2]))
    var.gamma.in=var.gamma.in+ncts[ctl]*t(a.m)%*%solve(var.i.c)%*%a.m
    sum.gamma <- sum.gamma + t(a.m)%*%solve(var.i.c)%*%(apply(theta.c[[ctl]],2,sum)
                                                        -ncts[ctl]*(lambda.c*b.c[ctl,]+mu.c))
  }
  if(lambda.c[2]==0&lambda.c[3]==0){
    gamma.c=rep(0,model.d)
  }else if(lambda.c[2]==0){
    var.gamma.in=var.gamma.in[-1,-1]
    var.gamma=solve(var.gamma.in)
    mean.gamma=var.gamma%*%sum.gamma[-1]
    gamma.c[1]=0
    gamma.c[2:3]=as.vector(mean.gamma+chol(var.gamma)%*%matrix(rnorm(model.d-1),model.d-1,1))
  }else if(lambda.c[3]==0){
    var.gamma.in=var.gamma.in[1,1]
    var.gamma=1/var.gamma.in
    mean.gamma=var.gamma*sum.gamma[1]
    gamma.c[1]=mean.gamma+sqrt(var.gamma)*rnorm(1)
    gamma.c[2:3]=0
  }else{
    var.gamma=solve(var.gamma.in)
    mean.gamma=var.gamma%*%sum.gamma
    gamma.c=as.vector(mean.gamma+chol(var.gamma)%*%matrix(rnorm(model.d),model.d,1))
  }
  var.s.c.r=rbind(c(1,0,0),c(gamma.c[1],1,0),c(gamma.c[2],gamma.c[3],1))
  var.s.c.l=diag(lambda.c)%*%var.s.c.r
  var.s.c=var.s.c.l%*%t(var.s.c.l)
  
  #save the results
  if(ii>nburn&((ii-nburn)/nskip==as.integer((ii-nburn)/nskip))){
    print(mu.c)
    mu.s = rbind(mu.s,mu.c)
    for(ctl in 1:nct){
      theta.s[[ctl]]=theta.s[[ctl]]+theta.c[[ctl]]
      alpha.s[ctl,]=alpha.s[ctl,]+mu.c+alpha.c[ctl,]
      alpha.s.s[ctl,]=alpha.s.s[ctl,]+(mu.c+alpha.c[ctl,])^2
      for(i in 1:ncts[ctl]){
        x.s[[ctl]]=x.s[[ctl]]+x.c[[ctl]][i,,]
      }
    }
    digit=4*(lambda.c[1]!=0)+2*(lambda.c[2]!=0)+(lambda.c[3]!=0)+1
    index.s[digit]=index.s[digit]+1
    do.digit=0
    if(do.digit){
      if(lambda.c[1]==0&lambda.c[2]==0&lambda.c[3]==0){
        index.s[1]=index.s[1]+1
      }else if(lambda.c[1]==0&lambda.c[2]==0){
        index.s[2]=index.s[2]+1
      }else if(lambda.c[1]==0&lambda.c[3]==0){
        index.s[3]=index.s[3]+1
      }else if(lambda.c[2]==0&lambda.c[3]==0){
        index.s[4]=index.s[4]+1
      }else if(lambda.c[1]==0){
        index.s[5]=index.s[5]+1
      }else if(lambda.c[2]==0){
        index.s[6]=index.s[6]+1
      }else if(lambda.c[3]==0){
        index.s[7]=index.s[7]+1
      }else{
        index.s[8]=index.s[8]+1
      }
    }
    icount = icount + 1
  }
}
index.s=index.s/icount
alpha.s=alpha.s/icount
alpha.s.s=alpha.s.s/icount
sd.s=sqrt(alpha.s.s-alpha.s^2)
for(ctl in 1:nct){
  x.s[[ctl]]=x.s[[ctl]]/icount
  theta.s[[ctl]]=theta.s[[ctl]]/icount
}

library(HDInterval)
hdi(mu.s[,3])
apply(mu.s,2,mean)
apply(mu.s,2,sd)
apply(mu.s,2,hdi)

require(ggplot2)
require(reshape2)
require(xtable)
library(vcd) 
virusload=NULL
for(ctl in 1:nct){
  for(i in 1:ncts[ctl]){
    if(ctl==1){
      id=i
    }else{
      id=i+sum(ncts[1:(ctl-1)])
    }
    tobs0 = rawdata[rawdata$Strain==cts[ctl]&rawdata$DuckID==id,"DayPI"]
    yobs0 = as.vector(rawdata[rawdata$Strain==cts[ctl]&rawdata$DuckID==id,"CloacaVirus"])
    #impute the observed data
    tobs0=tobs0[!is.na(yobs0)]
    yobs0=yobs0[!is.na(yobs0)]
    odesol=try(lsoda(y0.c[[ctl]][i,],tobs0,odeequations,theta.s[[ctl]][i,],atol=atolv,rtol=rtolv))
    virusload=c(virusload,odesol[,4])
  }
}

shapevals=c(1:7)
duckdata=rawdata[!is.na(rawdata$CloacaVirus),]
duckdata$virusload=virusload
ww=19.5/2.54; wh=ww/1; 
#size is 17.5cm, needs to be in inches 17.5/2.54 approx 7
windows(width=ww,height=wh)
strainnames=c("H3N8","H3N8*","H4N6","H4N8","H6N1","H6N2","H6N8")
pl1 <- ggplot(duckdata,aes(DayPI,CloacaVirus,colour=Strain,shape=Strain))
pl1=pl1+geom_point(size=2.5)+scale_shape_manual(values=shapevals)
pl1 <- pl1+facet_wrap(~DuckID,nrow=7,as.table=TRUE)+scale_x_continuous("Days",limits=c(0,15))
pl1 <- pl1+scale_y_continuous("virus titer (log10 EID/mL)",limits=c(0,8.5)) 
pl1 <-pl1+geom_line(aes(DayPI,virusload),size=1.5)
pl1 <- pl1 + theme_bw() + theme(strip.background = element_blank(),strip.text.x = element_blank())
print(pl1)
dev.copy2pdf(file="duck_flu.pdf")
