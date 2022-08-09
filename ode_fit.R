#ODE equation which can be solved numerically. t-time, y-observation, x-parameter
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

#ODE fit for the current parameter estimation pars.all
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

