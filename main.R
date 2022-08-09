##this code give direct ode parameter estimation method for flu virus model
##based on spline Bayesian MCMC method including
##derivative term and random effect variable selection
rm(list=ls())
graphics.off(); #close all graphics windows
require(deSolve)  #loads ODE solver package
require(mvtnorm) #multi-variate norm package 
require(nloptr) #Solve optimization problems for fitting
require(GenSA) #Performs search for global minimum of a very complex non-linear objective function 
require(splines2) #Constructs basis matrix of B-splines
require(LaplacesDemon) #Provides a complete environment for Bayesian inference
require(msm) #Functions for fitting continuous-time Markov Models

#directory of the data file
setwd("E:/paper/jcgs_paper/code")

#input ODE fit 
source("ode_fit.R")

#input ODE MCMC 
source("ode_mcmc.R")


#tolerance for ODE solution
atolv=1e-7; 
rtolv=1e-7;

#number of ODE variables
model.d=3

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

#knots can be chosen in different way according to the distribution of the observation
knots=c(seq(1.1,2.9,by=0.3),seq(3.1,12,by=3))
bmatrix=bSpline(tvec,knots=knots,degree=3,intercept =T) 
df=ncol(bmatrix)

#initialization
theta.s = list()
duck.id=0
for(ctl in 1:nct){
  y0.c[[ctl]]=rep(0,3)
  data.c[[ctl]] = list()
  x.c[[ctl]] = array(0,c(ncts[ctl],df,3))
  theta.c[[ctl]] = matrix(0,ncts[ctl],model.d)
  theta.s[[ctl]] = matrix(0,ncts[ctl],model.d)
  #ODE variable and parameter are initialized using the optimal solution based on the fit
  # to the observed data
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
##initial values for x and theta at population level
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

#MCMC procedure for total number of nrep iterations 
for(ii in 1:nrep){  
  #update x and theta for each individual
  for(ctl in 1:nct){
    for(i in 1:ncts[ctl]){
      odestack=data.c[[ctl]][[i]]
      tobs=odestack[,1]
      yobs=odestack[,-1]
      x.ini = x.c[[ctl]][i,,]
      theta.ini = theta.c[[ctl]][i,]
      odeout = gpode(tobs,yobs,x.ini,theta.ini,var.i.c,alpha.c[ctl,]+mu.c,knots)
      x.c[[ctl]][i,,] = odeout$x
      theta.c[[ctl]][i,] = odeout$theta
    }
  }
  ##update population level alpha.c and b.c
  b.c=NULL
  for(ctl in 1:nct){
    sum.theta <- apply(theta.c[[ctl]],2,sum)-ncts[ctl]*mu.c
    norm.var <- solve(ncts[ctl]*t(var.s.c.l)%*%solve(var.i.c)%*%var.s.c.l+diag(model.d))
    norm.mean <- norm.var%*%t(var.s.c.l)%*%solve(var.i.c)%*%sum.theta
    temp=as.vector(norm.mean+chol(norm.var)%*%matrix(rnorm(model.d),model.d,1))
    b.c=rbind(b.c,temp)
    alpha.c[ctl,]=var.s.c.l%*%temp
  }  
  
  ##update population level mu.c
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
  
  ##update population level var.i 
  sum.theta = matrix(0,model.d,model.d)
  for(ctl in 1:nct){
    ttvec <- t(theta.c[[ctl]]) - alpha.c[ctl,] - mu.c
    sum.theta <- sum.theta + ttvec%*%t(ttvec)
  }
  wishart.nu = bn + hyp.nu
  wishart.sigma = sum.theta + solve(hyp.omega)
  var.i.c =rinvwishart(wishart.nu, wishart.sigma)
  
  ##update population level var.s
  
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
  
  #sample from truncated normal for strain level random effect variable selection
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
    #index.s represent the model for different combination of random effect
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

#summary based on MCMC sampling results
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

#plot the fitted curve based on the Bayesian MCMC over the observed data for each strain
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
 
#save the plot into a file: duck_flu.pdf 
dev.copy2pdf(file="duck_flu.pdf")
