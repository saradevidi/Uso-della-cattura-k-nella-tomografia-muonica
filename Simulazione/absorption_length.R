
#
# definisco le curve che descrivono il lambda in funzione dell'energia del fotone per ciascun materiale 
# considerato -> le curve che mi restituiscono il lambda per una data energia del fotone sono :
# funTotVec.Z14 (curva per il Silicio)
# funTotVec.Z6 (curva per il Carbonio)
# funTotVec.Z26 (curva per il Ferro);
# Queste curve mi servono per calcolare la funzione di sopravvivenza Survival.lambda
#

# dati per il grafico per Z = 6 (Carbonio)
dataForGraph6 = matrix(c(1.00000E-03,	2.211E+03,	
                         1.50000E-03,	7.002E+02,	
                         2.00000E-03,	3.026E+02,	
                         3.00000E-03,	9.033E+01,	
                         4.00000E-03,	3.778E+01,	
                         5.00000E-03,	1.912E+01,	
                         6.00000E-03,	1.095E+01,
                         8.00000E-03,	4.576E+00,	
                         1.00000E-02,	2.373E+00,	
                         1.50000E-02,	8.071E-01,	
                         2.00000E-02,	4.420E-01,
                         3.00000E-02,	2.562E-01,	
                         4.00000E-02,	2.076E-01,	
                         5.00000E-02,	1.871E-01,	
                         6.00000E-02,	1.753E-01,	
                         8.00000E-02,	1.610E-01,	
                         1.00000E-01,	1.514E-01,	
                         1.50000E-01,	1.347E-01,	
                         2.00000E-01,	1.229E-01,	
                         3.00000E-01,	1.066E-01,	
                         4.00000E-01,	9.546E-02,	
                         5.00000E-01,	8.715E-02,	
                         6.00000E-01,	8.058E-02,	
                         8.00000E-01,	7.076E-02,	
                         1.00000E+00,	6.361E-02,	
                         1.25000E+00,	5.690E-02,	
                         1.50000E+00,	5.179E-02,	
                         2.00000E+00,	4.442E-02,	
                         3.00000E+00,	3.562E-02,	
                         4.00000E+00,	3.047E-02,	
                         5.00000E+00,	2.708E-02,	
                         6.00000E+00,	2.469E-02,	
                         8.00000E+00,	2.154E-02,	
                         1.00000E+01,	1.959E-02), ncol=2, byrow = T)
dataForGraph6 = as.data.frame(dataForGraph6)
colnames(dataForGraph6) = c("Energy.Mev", "mu_fratto_ro")
dataForGraph6$lambda = 1/dataForGraph6$mu_fratto_ro


# dati per il grafico Z = 14 (Silicio)
dataForGraph14 = matrix(c(1.00000E-03,	1.570E+03,
                          1.50000E-03,	5.355E+02,	
                          1.83890E-03,  3.092E+02,	
                          1.83890E-03,  3.192E+03,	
                          2.00000E-03,	2.777E+03,	
                          3.00000E-03,	9.784E+02,	
                          4.00000E-03,	4.529E+02,	
                          5.00000E-03,	2.450E+02,	
                          6.00000E-03,	1.470E+02,	
                          8.00000E-03,	6.468E+01,	
                          1.00000E-02,	3.389E+01,	
                          1.50000E-02,	1.034E+01,	
                          2.00000E-02,	4.464E+00,	
                          3.00000E-02,  1.436E+00,	
                          4.00000E-02,	7.012E-01,	
                          5.00000E-02,	4.385E-01,
                          6.00000E-02,	3.207E-01,	
                          8.00000E-02,	2.228E-01,	
                          1.00000E-01,	1.835E-01,	
                          1.50000E-01,	1.448E-01,	
                          2.00000E-01,	1.275E-01,	
                          3.00000E-01,	1.082E-01,	
                          4.00000E-01,	9.614E-02,	
                          5.00000E-01,	8.748E-02,	
                          6.00000E-01,	8.077E-02,	
                          8.00000E-01,	7.082E-02,	
                          1.00000E+00,	6.361E-02,	
                          1.25000E+00,	5.688E-02,	
                          1.50000E+00,	5.183E-02,	
                          2.00000E+00,	4.480E-02,	
                          3.00000E+00,	3.678E-02,	 
                          4.00000E+00,	3.240E-02,	
                          5.00000E+00,	2.967E-02,	
                          6.00000E+00,	2.788E-02,	
                          8.00000E+00,	2.574E-02,	
                          1.00000E+01,	2.462E-02), ncol=2, byrow = T)
dataForGraph14 = as.data.frame(dataForGraph14)
colnames(dataForGraph14) = c("Energy.Mev", "mu_fratto_ro")
dataForGraph14$lambda = 1/dataForGraph14$mu_fratto_ro


# dati per il grafico Z = 26 (Ferro)
dataForGraph26 = matrix(c(1.00000E-03,	9.085E+03,	
                          1.50000E-03,	3.399E+03,	
                          2.00000E-03,	1.626E+03,	
                          3.00000E-03,	5.576E+02,	
                          4.00000E-03,	2.567E+02,	
                          5.00000E-03,	1.398E+02,	
                          6.00000E-03,	8.484E+01,	
                          7.11200E-03,	5.319E+01,	
                          7.11200E-03,	4.076E+02,	
                          8.00000E-03,	3.056E+02,	
                          1.00000E-02,	1.706E+02,	
                          1.50000E-02,	5.708E+01,	
                          2.00000E-02,	2.568E+01,	
                          3.00000E-02,	8.176E+00,	
                          4.00000E-02,	3.629E+00,	
                          5.00000E-02,	1.958E+00,	
                          6.00000E-02,	1.205E+00,	
                          8.00000E-02,	5.952E-01,	
                          1.00000E-01,	3.717E-01,	
                          1.50000E-01,	1.964E-01,	
                          2.00000E-01,	1.460E-01,	
                          3.00000E-01,	1.099E-01,	
                          4.00000E-01,	9.400E-02,	
                          5.00000E-01,	8.414E-02,	
                          6.00000E-01,	7.704E-02,	
                          8.00000E-01,	6.699E-02,	
                          1.00000E+00,	5.995E-02,
                          1.25000E+00,	5.350E-02,	
                          1.50000E+00,	4.883E-02,	
                          2.00000E+00,	4.265E-02,	
                          3.00000E+00,	3.621E-02,	
                          4.00000E+00,	3.312E-02,	
                          5.00000E+00,	3.146E-02,	
                          6.00000E+00,	3.057E-02,	
                          8.00000E+00,	2.991E-02,	
                          1.00000E+01,	2.994E-02), ncol=2, byrow = T)
dataForGraph26 = as.data.frame(dataForGraph26)
colnames(dataForGraph26) = c("Energy.Mev", "mu_fratto_ro")
dataForGraph26$lambda = 1/dataForGraph26$mu_fratto_ro





#..............................................................
# curve che interpolano i punti osservati per ciascun materiale
#..............................................................


#..................
#  Z = 26 (Ferro)
#..................


# interpolazione con loess 


lo1.Z26 = loess(dataForGraph26$lambda[1:8]~dataForGraph26$Energy.Mev[1:8], span = 0.7,
            control = loess.control(surface = "direct"))


lo2.Z26 = loess(dataForGraph26$lambda[9:36]~dataForGraph26$Energy.Mev[9:36], span = 0.2,
           control = loess.control(surface = "direct")) 



funTot.Z26 = function(x)
{
  out=((if(length(x)>1){predict(lo1.Z26,data.frame(Energy.Mev=x))}else{predict(lo1.Z26,x)}) * I( x <= dataForGraph26$Energy.Mev[8] & x > 0) ) + 
    ((if(length(x)>1){predict(lo2.Z26,data.frame(Energy.Mev=x))}else{predict(lo2.Z26,x)}) * I( x > dataForGraph26$Energy.Mev[8] & x <= dataForGraph26$Energy.Mev[36]) ) +
    (predict(lo2.Z26,dataForGraph26$Energy.Mev[36])) * I( x > dataForGraph26$Energy.Mev[36] ) 
  return(out)
}
funTotVec.Z26 = Vectorize(funTot.Z26)






#...................
#  Z = 14 (Silicio)
#...................


# interpolazione con loess e regressione lineare

lm1.Z14 = lm(dataForGraph14$lambda[1:3] ~ dataForGraph14$Energy.Mev[1:3] + I(dataForGraph14$Energy.Mev[1:3]^2))

fun1.Z14 = function(x){ as.numeric(lm1.Z14$coefficients[1]) + as.numeric(lm1.Z14$coefficients[2] * x) + as.numeric(lm1.Z14$coefficients[3] * x^2) }



lo2.Z14 = loess(dataForGraph14$lambda[4:36]~dataForGraph14$Energy.Mev[4:36], span = 0.2,
               control = loess.control(surface = "direct")) 



funTot.Z14 = function(x)
{
  out=(fun1.Z14(x) * I( x <= dataForGraph14$Energy.Mev[3] & x > 0 ) ) + 
    ((if(length(x)>1){predict(lo2.Z14,data.frame(Energy.Mev=x))}else{predict(lo2.Z14,x)}) * I( x > dataForGraph14$Energy.Mev[3] & x <= dataForGraph14$Energy.Mev[36] ) ) +
    (predict(lo2.Z14,dataForGraph14$Energy.Mev[36])) * I( x > dataForGraph14$Energy.Mev[36] ) 
  return(out)
}
funTotVec.Z14 = Vectorize(funTot.Z14)







#...................
#  Z = 6 (Carbonio)
#...................




# interpolazione con loess 

lo.Z6 = loess(dataForGraph6$lambda~dataForGraph6$Energy.Mev, span = 0.2,
              control = loess.control(surface = "direct"))


funTot.Z6 = function(x) 
{
  out= ((if(length(x)>1){predict(lo.Z6,data.frame(Energy.Mev=x))}else{predict(lo.Z6,x)}) * I(x > 0 & x <= dataForGraph6$Energy.Mev[34]) ) +
    (predict(lo.Z6,dataForGraph6$Energy.Mev[34])) * I( x > dataForGraph6$Energy.Mev[34] )
  return(out)
}
funTotVec.Z6 = Vectorize(funTot.Z6)








# dataframe che contiene i nomi delle funzioni che servono per estrarre il lambda data una certa energia e un certo materiale
function.lambda = data.frame( name = c("funTotVec.Z26", "funTotVec.Z14", "funTotVec.Z6"), Z = c(26, 14, 6) ) 




