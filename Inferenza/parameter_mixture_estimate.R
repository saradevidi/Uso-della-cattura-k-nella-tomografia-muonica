EM.mixture = function(start, data, eps = 1e-5, max.iter = 100)
{
  abs.e = 1
  estimate = start # (k1, k2, k3, k4, delta, beta, mu1, sigma1, mu2, sigma2, mu3, sigma3) oppure
                   # (k1, k2, k3, delta, beta, mu1, sigma1, mu2, sigma2) 
  new.estimate = numeric(length(start))
  iter = 0
  
  if(length(estimate)==12)
  {# se la distribuzione ha 3 picchi
    
    
    # densita' dati osservati
    p.theta = function(theta, data)
    {
      out=  as.numeric( theta[1]*(dgamma(data,shape=theta[5],rate = theta[6])*I(data>0)) + 
                            theta[2]*((1/pnorm(q=(-theta[7]/theta[8]),lower.tail = FALSE))*dnorm(data, mean = theta[7], sd=theta[8])*I(data>0)) +
                            theta[3]*((1/pnorm(q=(-theta[9]/theta[10]),lower.tail = FALSE))*dnorm(data, mean = theta[9], sd=theta[10])*I(data>0)) +
                            theta[4]*((1/pnorm(q=(-theta[11]/theta[12]),lower.tail = FALSE))*dnorm(data, mean = theta[11], sd=theta[12])*I(data>0)) ) 
      out
    }
    
    # funzione che calcola costanti da utilizzare nell'M-step
    const = function(thetag,data)
    {
        mat = matrix(NA,nrow=length(data),ncol=4)
        mat[,1]= ((dgamma(data,shape=thetag[5],rate = thetag[6])*thetag[1])/p.theta(thetag,data))
        mat[,2]= (((1/pnorm(q=(-thetag[7]/thetag[8]),lower.tail = FALSE))*dnorm(data, mean = thetag[7], sd=thetag[8])*thetag[2])/p.theta(thetag,data))
        mat[,3]= (((1/pnorm(q=(-thetag[9]/thetag[10]),lower.tail = FALSE))*dnorm(data, mean = thetag[9], sd=thetag[10])*thetag[3])/p.theta(thetag,data))
        mat[,4]= (((1/pnorm(q=(-thetag[11]/thetag[12]),lower.tail = FALSE))*dnorm(data, mean = thetag[11], sd=thetag[12])*thetag[4])/p.theta(thetag,data))
        mat
    }
    
    # funzione da massimizzare per via numerica nell'M-step
    function.maximize = function(par, data, const.calc) # par  = (delta, beta, mu1, sigma1, mu2, sigma2, mu3, sigma3) 
    {

        -sum( (par[1]*log(par[2]) - log(gamma(par[1])) + par[1]*log(data) - par[2]*data)*const.calc[,1] +
                (-0.5*log(par[4]^2)-0.5*(data^2/par[4]^2)+data*(par[3]/par[4]^2)-0.5*(par[3]^2/par[4]^2)-log(pnorm(q=(-par[3]/par[4]),lower.tail = FALSE)))*const.calc[,2] +
                (-0.5*log(par[6]^2)-0.5*(data^2/par[6]^2)+data*(par[5]/par[6]^2)-0.5*(par[5]^2/par[6]^2)-log(pnorm(q=(-par[5]/par[6]),lower.tail = FALSE)))*const.calc[,3] +
                (-0.5*log(par[8]^2)-0.5*(data^2/par[8]^2)+data*(par[7]/par[8]^2)-0.5*(par[7]^2/par[8]^2)-log(pnorm(q=(-par[7]/par[8]),lower.tail = FALSE)))*const.calc[,4] )
    }
    
    while ((abs.e > eps) & (iter < max.iter)) 
    {
      iter = iter + 1
      
      # M-step per (delta, beta, mu1, sigma1, mu2, sigma2, mu3, sigma3) 
      const.calc = const(thetag = estimate, data=data)
      new.estimate = nlminb(start =  estimate[5:12], objective = function.maximize, lower = c(1e-7,   1e-7,   1e-7,  1e-7, 1e-7,  1e-7,   1e-7 ,  1e-7), 
                            upper = c(Inf, Inf,  Inf,  Inf,  Inf,  Inf,  Inf,  Inf), data = data, const.calc = const.calc,
                            control=list(iter.max=500,eval.max=500))$par 
      
      # M-step per (k1, k2, k3, k4) 
        new.estimate.k1 = mean( ((dgamma(data,shape=estimate[5],rate = estimate[6])*estimate[1])/p.theta(estimate,data)) )
        new.estimate.k2 = mean( (((1/pnorm(q=(-estimate[7]/estimate[8]),lower.tail = FALSE))*dnorm(data, mean = estimate[7], sd=estimate[8])*estimate[2])/p.theta(estimate,data)) )
        new.estimate.k3 = mean( (((1/pnorm(q=(-estimate[9]/estimate[10]),lower.tail = FALSE))*dnorm(data, mean = estimate[9], sd=estimate[10])*estimate[3])/p.theta(estimate,data)) )
        new.estimate.k4 = mean( (((1/pnorm(q=(-estimate[11]/estimate[12]),lower.tail = FALSE))*dnorm(data, mean = estimate[11], sd=estimate[12])*estimate[4])/p.theta(estimate,data)) )

      new.estimate = c(new.estimate.k1,new.estimate.k2,new.estimate.k3,new.estimate.k4,new.estimate) # (k1, k2, k3, k4, delta, beta, mu1, sigma1, mu2, sigma2, mu3, sigma3)
      
      abs.e = max(abs(new.estimate - estimate))
      estimate = new.estimate
    }
    
  }else{
    
    # densita' dati osservati
    p.theta = function(theta, data)
    {
      out=  as.numeric( theta[1]*(dgamma(data,shape=theta[4],rate = theta[5])*I(data>0)) + 
                            theta[2]*((1/pnorm(q=(-theta[6]/theta[7]),lower.tail = FALSE))*dnorm(data, mean = theta[6], sd=theta[7])*I(data>0)) +
                            theta[3]*((1/pnorm(q=(-theta[8]/theta[9]),lower.tail = FALSE))*dnorm(data, mean = theta[8], sd=theta[9])*I(data>0)) ) 
      out
    }
    
    # funzione che calcola costanti da utilizzare nell'M-step
    const = function(thetag,data)
    {
      mat = matrix(NA,nrow=length(data),ncol=3)
      mat[,1]= ((dgamma(data,shape=thetag[4],rate = thetag[5])*thetag[1])/p.theta(thetag,data))
      mat[,2]= (((1/pnorm(q=(-thetag[6]/thetag[7]),lower.tail = FALSE))*dnorm(data, mean = thetag[6], sd=thetag[7])*thetag[2])/p.theta(thetag,data))
      mat[,3]= (((1/pnorm(q=(-thetag[8]/thetag[9]),lower.tail = FALSE))*dnorm(data, mean = thetag[8], sd=thetag[9])*thetag[3])/p.theta(thetag,data))
      mat
    }
    
    # funzione da massimizzare per via numerica nell'M-step
    function.maximize = function(par, data, const.calc) # par  = (delta, beta, mu1, sigma1, mu2, sigma2)
    {
        -sum( (par[1]*log(par[2]) - log(gamma(par[1])) + par[1]*log(data) - par[2]*data)*const.calc[,1] +
                (-0.5*log(par[4]^2)-0.5*(data^2/par[4]^2)+data*(par[3]/par[4]^2)-0.5*(par[3]^2/par[4]^2)-log(pnorm(q=(-par[3]/par[4]),lower.tail = FALSE)))*const.calc[,2] +
                (-0.5*log(par[6]^2)-0.5*(data^2/par[6]^2)+data*(par[5]/par[6]^2)-0.5*(par[5]^2/par[6]^2)-log(pnorm(q=(-par[5]/par[6]),lower.tail = FALSE)))*const.calc[,3] )
    }
    
  
    while ((abs.e > eps) & (iter < max.iter)) 
    {
      iter = iter + 1
      
      # M-step per (delta, beta, mu1, sigma1, mu2, sigma2)
      const.calc = const(thetag = estimate, data=data)
      new.estimate = nlminb(start = estimate[4:9], objective = function.maximize, lower = c(1e-7,   1e-7,   1e-7,  1e-7, 1e-7,  1e-7), 
                            upper = c(Inf, Inf,  Inf,  Inf,  Inf,  Inf), data = data, const.calc = const.calc,
                            control=list(iter.max=500,eval.max=500))$par 
      
      # M-step per (k1, k2, k3)
        new.estimate.k1 = mean( ((dgamma(data,shape=estimate[4],rate = estimate[5])*estimate[1])/p.theta(estimate,data)) )
        new.estimate.k2 = mean( (((1/pnorm(q=(-estimate[6]/estimate[7]),lower.tail = FALSE))*dnorm(data, mean = estimate[6], sd=estimate[7])*estimate[2])/p.theta(estimate,data)) )
        new.estimate.k3 = mean( (((1/pnorm(q=(-estimate[8]/estimate[9]),lower.tail = FALSE))*dnorm(data, mean = estimate[8], sd=estimate[9])*estimate[3])/p.theta(estimate,data)) )
        
        
        new.estimate = c(new.estimate.k1,new.estimate.k2,new.estimate.k3,new.estimate) # (k1, k2, k3, delta, beta, mu1, sigma1, mu2, sigma2)
      
      abs.e = max(abs(new.estimate - estimate))
      estimate = new.estimate
    }
    
  }
  
  if (iter == max.iter) warning("Maximum number of iteration reached")
  
  list(estimate = estimate, iter = iter)
}
  
source("experimental_pdf_photon")  

# stimo i parametri per le distribuzioni delle energie di un dato numero atomico

estimate.Z26 = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,2,0.1,6,0.1,9.5,0.1),data = E.tot.26)
estimate.Z14 = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,1.5,0.1,5,0.1,7,0.1),data = E.tot.14)
estimate.Z6 = EM.mixture(start= c(1/3,1/3,1/3,0.1,0.1,1,0.1,3,0.1),data = E.tot.6)

estimate.Z26 = estimate.Z26$estimate # (k1, k2, k3, k4, delta, beta, mu1, sigma1, mu2, sigma2, mu3, sigma3)
estimate.Z14 = estimate.Z14$estimate # (k1, k2, k3, k4, delta, beta, mu1, sigma1, mu2, sigma2, mu3, sigma3)
estimate.Z6 = estimate.Z6$estimate # (k1, k2, k3, delta, beta, mu1, sigma1, mu2, sigma2)

save(estimate.Z26, file="estimate.Z26.RData")
save(estimate.Z14, file="estimate.Z14.RData")
save(estimate.Z6, file="estimate.Z6.RData")

weight.calculation = function(z)
{ # calcola i pesi k1.s, k2.s, k3.s, k4.s
  nmax = 4 # massimo dei picchi per le distribuzioni
  out = c(NA, nmax) 
  
  # pesi assunti noti delle distribuzioni mistura delle energie
  k1=c(estimate.Z6[1] , estimate.Z14[1] ,  estimate.Z26[1] )
  k2=c(estimate.Z6[2] , estimate.Z14[2] ,  estimate.Z26[2] )
  k3=c(estimate.Z6[3] , estimate.Z14[3] , estimate.Z26[3] )
  k4=c(0,   estimate.Z14[4] , estimate.Z26[4] )
  
  k1.s=  as.numeric(t(z)%*%k1)
  k2.s= as.numeric(t(z)%*%k2)
  k3.s= as.numeric(t(z)%*%k3)
  k4.s= as.numeric(t(z)%*%k4)
  
  out = c(k1.s, k2.s, k3.s, k4.s)
  out
}


fE.cond_sz = function(E, z) 
{ # distribuzione delle energie condizionata allo spessore e al numero atomico
  delta = c(estimate.Z6[4] ,estimate.Z14[5] ,estimate.Z26[5] ) 
  beta = c(estimate.Z6[5] ,estimate.Z14[6] ,estimate.Z26[6] )
  mu1=c(estimate.Z6[6] ,estimate.Z14[7] ,estimate.Z26[7] ) # primo picco
  sigma1=c(estimate.Z6[7] ,estimate.Z14[8] ,estimate.Z26[8] ) # primo picco
  mu2=c(estimate.Z6[8]  ,estimate.Z14[9] ,estimate.Z26[9] ) # secondo picco
  sigma2=c(estimate.Z6[9],estimate.Z14[10] ,estimate.Z26[10] ) # secondo picco
  mu3=c(0,estimate.Z14[11] ,estimate.Z26[11] ) # terzo picco
  sigma3=c(0,estimate.Z14[12],estimate.Z26[12]) # terzo picco
  
  Kj.s = weight.calculation(z=z)
  
  out= as.numeric(Kj.s[1] * ( dgamma(E,shape=t(z)%*%delta,rate = t(z)%*%beta)*I(E>0) ) + 
                    Kj.s[2] * ( (1/pnorm(q=(-as.numeric(t(z)%*%mu1)/as.numeric(t(z)%*%sigma1)),lower.tail = FALSE))*dnorm(E, mean = as.numeric(t(z)%*%mu1), sd=as.numeric(t(z)%*%sigma1))*I(E>0) ) +
                    Kj.s[3] * ( (1/pnorm(q=(-as.numeric(t(z)%*%mu2)/as.numeric(t(z)%*%sigma2)),lower.tail = FALSE))*dnorm(E, mean = as.numeric(t(z)%*%mu2), sd=as.numeric(t(z)%*%sigma2))*I(E>0) ) +
                    Kj.s[4] * replace( (1/pnorm(q=(-as.numeric(t(z)%*%mu3)/as.numeric(t(z)%*%sigma3)),lower.tail = FALSE))*dnorm(E, mean = as.numeric(t(z)%*%mu3), sd=as.numeric(t(z)%*%sigma3))*I(E>0), 
                                       is.na((1/pnorm(q=(-as.numeric(t(z)%*%mu3)/as.numeric(t(z)%*%sigma3)),lower.tail = FALSE))*dnorm(E, mean = as.numeric(t(z)%*%mu3), sd=as.numeric(t(z)%*%sigma3))*I(E>0)),
                                       0) ) 
  
  
  out
}
