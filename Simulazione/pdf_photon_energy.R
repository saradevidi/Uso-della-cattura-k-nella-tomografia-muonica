###### definizione delle pdf delle energie dei fotoni per ogni materiale


pdf.Z6 = function(x, E0 = 1, E1= 0.8, E2= 3, s1=0.05, s2=0.05, 
                  k1 = 0.7, k2 = 0.2, k3 = 0.1)
{
  k1 *( ((1/E0)*exp(-x/E0)) * I(x>0) ) + k2 * ( (1/pnorm(q=(-E1/s1),lower.tail = FALSE))*dnorm(x, mean = E1, sd=s1)*I(x>0) ) +
    k3 * ( (1/pnorm(q=(-E2/s2),lower.tail = FALSE))*dnorm(x, mean = E2, sd=s2)*I(x>0) )
}




pdf.Z14 = function(x, E0 = 1, E1=1.5, E2=5, E3 = 7, s1=0.05, s2=0.05, s3 = 0.05,
                   k1 = 0.7, k2 = 0.1, k3 = 0.05, k4 = 0.15)
{
  k1 * ((1/E0)*exp(-x/E0)*I(x>0)) + k2 * ( (1/pnorm(q=(-E1/s1),lower.tail = FALSE))*dnorm(x, mean = E1, sd=s1)*I(x>0) ) +
    k3 * ( (1/pnorm(q=(-E2/s2),lower.tail = FALSE))*dnorm(x, mean = E2, sd=s2)*I(x>0) ) +
    k4 * ( (1/pnorm(q=(-E3/s3),lower.tail = FALSE))*dnorm(x, mean = E3, sd=s3)*I(x>0) )
}




pdf.Z26 = function(x, E0 = 1, E1=2, E2=6, E3 = 9.5, s1=0.05, s2=0.05, s3 = 0.05,
                   k1 = 0.7, k2 = 0.05, k3 = 0.07, k4 = 0.18)
{
  k1 * ((1/E0)*exp(-x/E0)*I(x>0)) + k2 * ( (1/pnorm(q=(-E1/s1),lower.tail = FALSE))*dnorm(x, mean = E1, sd=s1)*I(x>0) ) +
    k3 * ( (1/pnorm(q=(-E2/s2),lower.tail = FALSE))*dnorm(x, mean = E2, sd=s2)*I(x>0) ) +
    k4 * ( (1/pnorm(q=(-E3/s3),lower.tail = FALSE))*dnorm(x, mean = E3, sd=s3)*I(x>0) )
}















###### generazione dalle pdf delle energie dei fotoni per ogni materiale



r.pdf.Z6 = function(n,
                    E0 = 1, E1= 0.8, E2= 3, s1=0.05, s2=0.05, 
                    k1 = 0.7, k2 = 0.2, k3 = 0.1)
{ # funzione che genera n valori casuali dalla pdf per il nr.atomico Z = 6 (pdf.Z6)
  
  out = rep(NA,n)
  
  for(i in 1:n)
  {
    # 1. genero un valore j da una multinomiale che assume i valori 1,2,3 con probabilità k1, k2, k3
    index = sample(1:3, size=1, prob= c(k1,k2,k3))
    
    # 2. genero un valore dalla densità f_j, cioè dalla densità campionata al passo 1.
    if(index == 1){
      out[i] = rexp(1, rate = 1/E0) # la prima densità mistura è una esponenziale di parametro 1/E0
    }else if(index==2){
      out[i] = (( qnorm(p=((1-pnorm(-E1/s1))*runif(1)+pnorm(-E1/s1))) )*s1) + E1 # la seconda densità mistura è gaussiana di parametri E1, s1 troncata tra 0 e +Inf
    }else{
      out[i] = (( qnorm(p=((1-pnorm(-E2/s2))*runif(1)+pnorm(-E2/s2))) )*s2) + E2 # la terza densità mistura è gaussiana di parametri E2, s2 troncata tra 0 e +Inf
    }
  }
  out
}




r.pdf.Z14 = function(n,
                     E0 = 1, E1=1.5, E2=5, E3 = 7, s1=0.05, s2=0.05, s3 = 0.05,
                     k1 = 0.7, k2 = 0.1, k3 = 0.05, k4 = 0.15)
{ # funzione che genera n valori casuali dalla pdf per il nr.atomico Z = 14 (pdf.Z14)
  
  out = rep(NA,n)
  
  for(i in 1:n)
  {
    # 1. genero un valore j da una multinomiale che assume i valori 1,2,3 con probabilità k1, k2, k3
    index = sample(1:4, size=1, prob= c(k1,k2,k3,k4))
    
    # 2. genero un valore dalla densità f_j, cioè dalla densità campionata al passo 1.
    if(index == 1){
      out[i] = rexp(1, rate = 1/E0) # la prima densità mistura è una esponenziale di parametro 1/E0
    }else if(index==2){
      out[i] = (( qnorm(p=((1-pnorm(-E1/s1))*runif(1)+pnorm(-E1/s1))) )*s1) + E1 # la seconda densità mistura è gaussiana di parametri E1, s1 troncata tra 0 e +Inf
    }else if(index==3){
      out[i] = (( qnorm(p=((1-pnorm(-E2/s2))*runif(1)+pnorm(-E2/s2))) )*s2) + E2 # la terza densità mistura è gaussiana di parametri E2, s2 troncata tra 0 e +Inf
    }else{
      out[i] = (( qnorm(p=((1-pnorm(-E3/s3))*runif(1)+pnorm(-E3/s3))) )*s3) + E3 # la quarta densità mistura è gaussiana di parametri E3, s3 troncata tra 0 e +Inf
    }
  }
  out
}




r.pdf.Z26 = function(n,
                     E0 = 1, E1=2, E2=6, E3 = 9.5, s1=0.05, s2=0.05, s3 = 0.05,
                     k1 = 0.7, k2 = 0.05, k3 = 0.07, k4 = 0.18)
{ # funzione che genera n valori casuali dalla pdf per il nr.atomico Z = 14 (pdf.Z14)
  
  out = rep(NA,n)
  
  for(i in 1:n)
  {
    # 1. genero un valore j da una multinomiale che assume i valori 1,2,3 con probabilità k1, k2, k3
    index = sample(1:4, size=1, prob= c(k1,k2,k3,k4))
    
    # 2. genero un valore dalla densità f_j, cioè dalla densità campionata al passo 1.
    if(index == 1){
      out[i] = rexp(1, rate = 1/E0) # la prima densità mistura è una esponenziale di parametro 1/E0
    }else if(index==2){
      out[i] = (( qnorm(p=((1-pnorm(-E1/s1))*runif(1)+pnorm(-E1/s1))) )*s1) + E1 # la seconda densità mistura è gaussiana di parametri E1, s1 troncata tra 0 e +Inf
    }else if(index==3){
      out[i] = (( qnorm(p=((1-pnorm(-E2/s2))*runif(1)+pnorm(-E2/s2))) )*s2) + E2 # la terza densità mistura è gaussiana di parametri E2, s2 troncata tra 0 e +Inf
    }else{
      out[i] = (( qnorm(p=((1-pnorm(-E3/s3))*runif(1)+pnorm(-E3/s3))) )*s3) + E3 # la quarta densità mistura è gaussiana di parametri E3, s3 troncata tra 0 e +Inf
    }
  }
  out
}



# matrice per i materiali e le corrispondenti funzioni
Z = c(6, 14, 26)
rfun.name = c("r.pdf.Z6", "r.pdf.Z14", "r.pdf.Z26")
pdf.fun = data.frame(Z, rfun.name)

