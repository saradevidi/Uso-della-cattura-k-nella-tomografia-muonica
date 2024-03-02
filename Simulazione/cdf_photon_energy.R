###### definizione delle cdf delle energie dei fotoni per ogni materiale



cdf.Z6 = function(q, E0 = 1, E1= 0.8, E2= 3, s1=0.05, s2=0.05, 
                  k1 = 0.7, k2 = 0.2, k3 = 0.1)
{
  if(q>=0)
  {
    k1 *( 1-exp(-q/E0) ) + k2 * ( (1/pnorm(q=(-E1/s1),lower.tail = FALSE))*( pnorm(q=((q-E1)/s1)) - pnorm(q=(-E1/s1)) ) ) +
      k3 * ( (1/pnorm(q=(-E2/s2),lower.tail = FALSE))*( pnorm(q=((q-E2)/s2)) - pnorm(q=(-E2/s2)) ) ) 
  }else{return(0)}
}

cdf.Z6 = Vectorize(cdf.Z6)



cdf.Z14 = function(q, E0 = 1, E1=1.5, E2=5, E3 = 7, s1=0.05, s2=0.05, s3 = 0.05,
                   k1 = 0.7, k2 = 0.1, k3 = 0.05, k4 = 0.15)
{
  if(q>=0)
  {
    k1 *( 1-exp(-q/E0) ) + k2 * ( (1/pnorm(q=(-E1/s1),lower.tail = FALSE))*( pnorm(q=((q-E1)/s1)) - pnorm(q=(-E1/s1)) ) ) +
      k3 * ( (1/pnorm(q=(-E2/s2),lower.tail = FALSE))*( pnorm(q=((q-E2)/s2)) - pnorm(q=(-E2/s2)) ) ) +
      k4 * ( (1/pnorm(q=(-E3/s3),lower.tail = FALSE))*( pnorm(q=((q-E3)/s3)) - pnorm(q=(-E3/s3)) ) ) 
  }else{return(0)}
}

cdf.Z14 = Vectorize(cdf.Z14)



cdf.Z26 = function(q, E0 = 1, E1=2, E2=6, E3 = 9.5, s1=0.05, s2=0.05, s3 = 0.05,
                   k1 = 0.7, k2 = 0.05, k3 = 0.07, k4 = 0.18)
{
  if(q>=0)
  {
    k1 *( 1-exp(-q/E0) ) + k2 * ( (1/pnorm(q=(-E1/s1),lower.tail = FALSE))*( pnorm(q=((q-E1)/s1)) - pnorm(q=(-E1/s1)) ) ) +
      k3 * ( (1/pnorm(q=(-E2/s2),lower.tail = FALSE))*( pnorm(q=((q-E2)/s2)) - pnorm(q=(-E2/s2)) ) ) +
      k4 * ( (1/pnorm(q=(-E3/s3),lower.tail = FALSE))*( pnorm(q=((q-E3)/s3)) - pnorm(q=(-E3/s3)) ) ) 
  }else{return(0)}
}

cdf.Z26 = Vectorize(cdf.Z26)





