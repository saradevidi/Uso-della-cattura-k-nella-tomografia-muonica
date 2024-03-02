##################################
# EM applicato a tutto il volume #
##################################


# carico le funzioni 
source("create_geometry.R")
source("muon_energy_loss.R")
source("pdf_photon_energy.R")
source("absorption_length.R")
source("function_generatekCapture.R")



# simulo i dati

cube1 = coord.volume(coord.ang = c(24,0,0), alt = 4, lato1 = 4, lato2 = 4)
cube2 = coord.volume.Extreme(coord.volume.matrix = cube1)
cube3 = coord.volume.Z(coord.volume.Extreme.matrix = cube2, atoms = c(6,14,26), seed_ = 123)
cube3$Z[4]=26 

coord.riv1 = coord.riv(x.coord = 3.5, coord.volume = cube1, lato2 = 4, alt = 4) 

densitaZ = matrix(c(6, 14, 26, 2.25, 2.33, 7.87), ncol=2) 

DeltaE.matrix = loss.energy(Zdensity = densitaZ, Y = 1)

kcapturePROB.matrix = cbind(DeltaE.matrix$Z,c(0.90, 0.85, 0.98))
colnames(kcapturePROB.matrix) = c("Z", "prob_cattura")
kcapturePROB.matrix = as.data.frame(kcapturePROB.matrix)  

Survival.lambda = function(spessore, densita, lambda) 
{
  exp( -(spessore * densita) / lambda)
}

Simulation = generate.kCapture(n.muon = 9000000, E_mu = c(2,40), coord.material = cube3, 
                               coord.detector = coord.riv1, DeltaE.matrix = DeltaE.matrix, 
                               kcapturePROB.matrix = kcapturePROB.matrix,
                               lato1 = 4, x = 0, seed_ = 29000)

Simulation = Simulation[(!is.na(Simulation$E_fot1) | !is.na(Simulation$E_fot2) | !is.na(Simulation$E_fot3)| !is.na(Simulation$E_fot4)) ,] # prendo solo i muoni per i quali si è visto almeno un fotone

# carico i dati ottenuti
load("data_for_inference_only_photon.RData")






# 1. Creo la funzione che mi costruisce la matrice dei dati necessaria all'EM per l'identificazione di ciascuno strato di voxel a y e z fssato

create.data = function(pixel.detector.muon = 1, # pixel del detector dei muoni entranti
                       Simulation = Simulation) # dataset ottenuto dalla funzione generate.kCapture
{
  
  # 1. Isolo tutti i muoni (con le corrispondenti energie associate dei fotoni emessi) che sono stati visti nel pixel pixel.detector.muon
  
  Simulation = Simulation[!is.na(Simulation$vox_detector_muon),]
  Simulation = Simulation[Simulation$vox_detector_muon==pixel.detector.muon,]
  
  # 2. Creo i dati osservati: prendo tutte le energie viste nel pixel 1, tutte le energie viste nel pixel 2, etc. con la corrspondente energia del muone entrante
  
  photon1 = Simulation[, c("E_fot1", "vox_detector_fot1", "E_muon")]
  photon1 = photon1[!is.na(photon1$E_fot1),]
  
  Energy.pixel.ph1 = matrix(NA, nrow=dim(Simulation)[1], ncol=(dim(coord.riv1)[1]*2))
  colnames(Energy.pixel.ph1) = c("pixel_1","pixel_2", "pixel_3", "pixel_4", "pixel_5", "pixel_6", "pixel_7","pixel_8", "pixel_9", "pixel_10", "pixel_11", "pixel_12",
                                 "pixel_13","pixel_14", "pixel_15", "pixel_16",
                                 "E_muon_pixel_1","E_muon_pixel_2","E_muon_pixel_3","E_muon_pixel_4","E_muon_pixel_5","E_muon_pixel_6","E_muon_pixel_7","E_muon_pixel_8","E_muon_pixel_9",
                                 "E_muon_pixel_10","E_muon_pixel_11","E_muon_pixel_12","E_muon_pixel_13","E_muon_pixel_14","E_muon_pixel_15","E_muon_pixel_16")
  for(i in 1:dim(coord.riv1)[1])
  {
    Energy.pixel.ph1[,i] = c(photon1$E_fot1[photon1$vox_detector_fot1==i], rep(NA,dim(Simulation)[1]-length(photon1$E_fot1[photon1$vox_detector_fot1==i])))
  }
  
  for(i in 1:dim(coord.riv1)[1])
  {
    Energy.pixel.ph1[,i+16] = c(photon1$E_muon[photon1$vox_detector_fot1==i], rep(NA,dim(Simulation)[1]-length(photon1$E_muon[photon1$vox_detector_fot1==i])))
  }
  
  photon2 = Simulation[, c("E_fot2", "vox_detector_fot2", "E_muon")]
  photon2 = photon2[!is.na(photon2$E_fot2),]
  
  Energy.pixel.ph2 = matrix(NA, nrow=dim(Simulation)[1], ncol=(dim(coord.riv1)[1]*2)) 
  colnames(Energy.pixel.ph2) = c("pixel_1","pixel_2", "pixel_3", "pixel_4", "pixel_5", "pixel_6", "pixel_7","pixel_8", "pixel_9", "pixel_10", "pixel_11", "pixel_12",
                                 "pixel_13","pixel_14", "pixel_15", "pixel_16",
                                 "E_muon_pixel_1","E_muon_pixel_2","E_muon_pixel_3","E_muon_pixel_4","E_muon_pixel_5","E_muon_pixel_6","E_muon_pixel_7","E_muon_pixel_8","E_muon_pixel_9",
                                 "E_muon_pixel_10","E_muon_pixel_11","E_muon_pixel_12","E_muon_pixel_13","E_muon_pixel_14","E_muon_pixel_15","E_muon_pixel_16")
  for(i in 1:dim(coord.riv1)[1])
  {
    Energy.pixel.ph2[,i] = c(photon2$E_fot2[photon2$vox_detector_fot2==i], rep(NA,dim(Simulation)[1]-length(photon2$E_fot2[photon2$vox_detector_fot2==i])))
  }
  
  for(i in 1:dim(coord.riv1)[1])
  {
    Energy.pixel.ph2[,i+16] = c(photon2$E_muon[photon2$vox_detector_fot2==i], rep(NA,dim(Simulation)[1]-length(photon2$E_muon[photon2$vox_detector_fot2==i])))
  }
  
  photon3 = Simulation[, c("E_fot3", "vox_detector_fot3", "E_muon")]
  photon3 = photon3[!is.na(photon3$E_fot3),]
  
  Energy.pixel.ph3 = matrix(NA, nrow=dim(Simulation)[1], ncol=(dim(coord.riv1)[1]*2))
  colnames(Energy.pixel.ph3) = c("pixel_1","pixel_2", "pixel_3", "pixel_4", "pixel_5", "pixel_6", "pixel_7","pixel_8", "pixel_9", "pixel_10", "pixel_11", "pixel_12",
                                 "pixel_13","pixel_14", "pixel_15", "pixel_16",
                                 "E_muon_pixel_1","E_muon_pixel_2","E_muon_pixel_3","E_muon_pixel_4","E_muon_pixel_5","E_muon_pixel_6","E_muon_pixel_7","E_muon_pixel_8","E_muon_pixel_9",
                                 "E_muon_pixel_10","E_muon_pixel_11","E_muon_pixel_12","E_muon_pixel_13","E_muon_pixel_14","E_muon_pixel_15","E_muon_pixel_16")
  for(i in 1:dim(coord.riv1)[1])
  {
    Energy.pixel.ph3[,i] = c(photon3$E_fot3[photon3$vox_detector_fot3==i], rep(NA,dim(Simulation)[1]-length(photon3$E_fot3[photon3$vox_detector_fot3==i])))
  }
  
  for(i in 1:dim(coord.riv1)[1])
  {
    Energy.pixel.ph3[,i+16] = c(photon3$E_muon[photon3$vox_detector_fot3==i], rep(NA,dim(Simulation)[1]-length(photon3$E_muon[photon3$vox_detector_fot3==i])))
  }
  
  photon4 = Simulation[, c("E_fot4", "vox_detector_fot4", "E_muon")]
  photon4 = photon4[!is.na(photon4$E_fot4),]
  
  Energy.pixel.ph4 = matrix(NA, nrow=dim(Simulation)[1], ncol=(dim(coord.riv1)[1]*2)) 
  colnames(Energy.pixel.ph4) = c("pixel_1","pixel_2", "pixel_3", "pixel_4", "pixel_5", "pixel_6", "pixel_7","pixel_8", "pixel_9", "pixel_10", "pixel_11", "pixel_12",
                                 "pixel_13","pixel_14", "pixel_15", "pixel_16",
                                 "E_muon_pixel_1","E_muon_pixel_2","E_muon_pixel_3","E_muon_pixel_4","E_muon_pixel_5","E_muon_pixel_6","E_muon_pixel_7","E_muon_pixel_8","E_muon_pixel_9",
                                 "E_muon_pixel_10","E_muon_pixel_11","E_muon_pixel_12","E_muon_pixel_13","E_muon_pixel_14","E_muon_pixel_15","E_muon_pixel_16")
  for(i in 1:dim(coord.riv1)[1])
  {
    Energy.pixel.ph4[,i] = c(photon4$E_fot4[photon4$vox_detector_fot4==i], rep(NA,dim(Simulation)[1]-length(photon4$E_fot4[photon4$vox_detector_fot4==i])))
  }
  
  for(i in 1:dim(coord.riv1)[1])
  {
    Energy.pixel.ph4[,i+16] = c(photon4$E_muon[photon4$vox_detector_fot4==i], rep(NA,dim(Simulation)[1]-length(photon4$E_muon[photon4$vox_detector_fot4==i])))
  }
  
  Energy.pixel = rbind(Energy.pixel.ph1, Energy.pixel.ph2, Energy.pixel.ph3, Energy.pixel.ph4)
  
  
  # matrice dei dati
  
  Energy.pixel = as.data.frame(Energy.pixel)
  
  data = c(0,0,0)
  for(pix in 1:16)
  {
    data1 = cbind(Energy.pixel[,pix][!is.na(Energy.pixel[,pix])], Energy.pixel[,16+pix][!is.na(Energy.pixel[,16+pix])], rep(pix,length(Energy.pixel[,16+pix][!is.na(Energy.pixel[,16+pix])])))
    data = rbind(data, data1)
  }
  data = data[-1,]
  
  colnames(data) = c("E", "E_muon", "pix(i)")
  data = as.data.frame(data)
  
  data
}





# stime dei parametri delle distribuzioni delle energie condizionate allo spessore e al numero atomico da usare nell'EM (ottenute in precedenza)
load("estimate.Z26.RData")
load("estimate.Z14.RData")
load("estimate.Z6.RData")

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


# funzione che calcola i valori attesi per l'E-step

Expectation = function(estimate, data)
{
  
  # calcolo theta_il1, theta_il2, theta_il3, theta_il4 per ogni energia del fotone E_il
  
  # prendo in input E_il e restituisco il theta_il1 corrispondente
  theta_il1 = I(data$E_muon>11.31)*10^-15 + 
    I((data$E_muon>4.16) & (data$E_muon<=11.31))*(1-estimate[1]-estimate[2]) + 
    I((data$E_muon>3.75 & (data$E_muon<=4.16)))*(1-estimate[1]-estimate[2] + estimate[1]) + 
    I(data$E_muon<=3.75)*(1-10^-15)
  
  # prendo in input E_il e restituisco il theta_il2 corrispondente
  theta_il2 = I(data$E_muon>22.62)*10^-15 + 
    I((data$E_muon>15.47) & (data$E_muon<=22.62))*(1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4]) + 
    I((data$E_muon>15.06 & (data$E_muon<=15.47)))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4]) + (1-estimate[1]-estimate[2])*estimate[3] + estimate[1]*(1-estimate[3]-estimate[4]) ) + 
    I((data$E_muon>11.31 & (data$E_muon<=15.06)))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4]) + (1-estimate[1]-estimate[2])*estimate[3] + estimate[1]*(1-estimate[3]-estimate[4]) + (1-estimate[1]-estimate[2])*estimate[4] + estimate[2]*(1-estimate[3]-estimate[4]) ) +
    I((data$E_muon>8.32 & (data$E_muon<=11.31)))*( estimate[2]*(1-estimate[3]-estimate[4]) + estimate[1]*(1-estimate[3]-estimate[4]) ) +
    I((data$E_muon>7.91 & (data$E_muon<=8.32)))*( estimate[2]*(1-estimate[3]-estimate[4]) + estimate[1]*(1-estimate[3]-estimate[4]) + estimate[1]*estimate[3] ) +
    I((data$E_muon>7.5 & (data$E_muon<=7.91)))*( estimate[2]*(1-estimate[3]-estimate[4]) + estimate[1]*(1-estimate[3]-estimate[4]) + estimate[1]*estimate[3] + estimate[1]*estimate[4] + estimate[2]*estimate[3] ) +
    I((data$E_muon>4.16 & (data$E_muon<=7.5)))*( estimate[2]*(1-estimate[3]-estimate[4]) + estimate[1]*(1-estimate[3]-estimate[4]) + estimate[1]*estimate[3] + estimate[1]*estimate[4] + estimate[2]*estimate[3] + estimate[2]*estimate[4] ) +
    I((data$E_muon>3.75 & (data$E_muon<=4.16)))*( estimate[2]*estimate[4] + estimate[2]*(1-estimate[3]-estimate[4]) + estimate[2]*estimate[3] ) +
    I(data$E_muon<=3.75)*10^-15
  
  # prendo in input E_il e restituisco il theta_il3 corrispondente
  theta_il3 = I(data$E_muon>33.93)*10^-15 + 
    I((data$E_muon>26.78) & (data$E_muon<=33.93))*(1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + 
    I((data$E_muon>26.37 & (data$E_muon<=26.78)))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5] + (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) ) + 
    I((data$E_muon>22.62 & (data$E_muon<=26.37)))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6] + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5] + (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) ) +
    I((data$E_muon>19.63 & (data$E_muon<=22.62)))*( (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) ) +
    I((data$E_muon>19.22 & (data$E_muon<=19.63)))*( (1-estimate[1]-estimate[2])*estimate[3]*estimate[5] + estimate[1]*(1-estimate[3]-estimate[4])*estimate[5] + estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) ) +
    I((data$E_muon>18.81 & (data$E_muon<=19.22)))*( (1-estimate[1]-estimate[2])*estimate[3]*estimate[6] + estimate[1]*(1-estimate[3]-estimate[4])*estimate[6] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[5] + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[3]*estimate[5] + estimate[1]*(1-estimate[3]-estimate[4])*estimate[5] + estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) ) +
    I((data$E_muon>15.47 & (data$E_muon<=18.81)))*( (1-estimate[1]-estimate[2])*estimate[4]*estimate[6] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6] + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[3]*estimate[6] + estimate[1]*(1-estimate[3]-estimate[4])*estimate[6] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[5] + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[3]*estimate[5] + estimate[1]*(1-estimate[3]-estimate[4])*estimate[5] + estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) ) +
    I((data$E_muon>15.06 & (data$E_muon<=15.47)))*( (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6]) + 
                                                      estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[5] + 
                                                      estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + (1-estimate[1]-estimate[2])*estimate[4]*estimate[6] +
                                                      estimate[2]*(1-estimate[3]-estimate[4])*estimate[6] + estimate[2]*estimate[4]*( 1-estimate[5]-estimate[6]) ) +
    I((data$E_muon>12.48 & (data$E_muon<=15.06)))*( estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])) +
    I((data$E_muon>12.07 & (data$E_muon<=12.48)))*( estimate[1]*estimate[3]*estimate[5] + estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])) +
    I((data$E_muon>11.66 & (data$E_muon<=12.07)))*( estimate[1]*estimate[3]*estimate[6] + estimate[1]*estimate[4]*estimate[5] + estimate[2]*estimate[3]*estimate[5] + 
                                                      estimate[1]*estimate[3]*estimate[5] + estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6]) ) +
    I((data$E_muon>11.25 & (data$E_muon<=11.66)))*( estimate[1]*estimate[3]*estimate[6] + estimate[1]*estimate[4]*estimate[5] + estimate[2]*estimate[3]*estimate[5] + 
                                                      estimate[1]*estimate[3]*estimate[5] + estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6]) +
                                                      estimate[1]*estimate[4]*estimate[6] + estimate[2]*estimate[4]*estimate[5] + estimate[2]*estimate[3]*estimate[6]) +
    I((data$E_muon>8.32 & (data$E_muon<=11.25)))*( estimate[1]*estimate[3]*estimate[6] + estimate[1]*estimate[4]*estimate[5] + estimate[2]*estimate[3]*estimate[5] + 
                                                     estimate[1]*estimate[3]*estimate[5] + estimate[1]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6]) +
                                                     estimate[1]*estimate[4]*estimate[6] + estimate[2]*estimate[4]*estimate[5] + estimate[2]*estimate[3]*estimate[6] + estimate[2]*estimate[4]*estimate[6]) +
    I((data$E_muon>7.91 & (data$E_muon<=8.32)))*( estimate[2]*estimate[4]*estimate[6] + estimate[2]*estimate[3]*estimate[6] + estimate[1]*estimate[4]*estimate[6] + estimate[2]*estimate[4]*estimate[5] + 
                                                    estimate[1]*estimate[4]*estimate[5] + estimate[2]*estimate[3]*estimate[5] + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6]) + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])) +
    I((data$E_muon>7.5 & (data$E_muon<=7.91)))*( estimate[2]*estimate[4]*estimate[6] + estimate[2]*estimate[4]*estimate[5] + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6]) ) +
    I(data$E_muon<=7.5)*10^-15
  
  # prendo in input E_il e restituisco il theta_il4 corrispondente
  theta_il4 = I(data$E_muon>45.24)*10^-15 + 
    I((data$E_muon>38.09) & (data$E_muon<=45.24))*(1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + 
    I((data$E_muon>37.68) & (data$E_muon<=38.09))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +(1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[7] + (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) + 
                                                      (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) ) + 
    I((data$E_muon>33.93) & (data$E_muon<=37.68))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +(1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[7] + (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) + 
                                                      (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + 
                                                      estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[8] ) + 
    I((data$E_muon>30.53) & (data$E_muon<=33.93))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) ) + 
    I((data$E_muon>30.12) & (data$E_muon<=30.53))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + 
                                                      estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5]*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[7] +
                                                      estimate[2]*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*estimate[7] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6])*estimate[8] + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*estimate[7] +
                                                      estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[8] + estimate[1]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8])  )  + 
    I((data$E_muon>26.78) & (data$E_muon<=30.12))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + 
                                                      estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[5]*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[7] +
                                                      estimate[2]*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*estimate[7] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[3]*(1-estimate[5]-estimate[6])*estimate[8] + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*estimate[7] +
                                                      estimate[1]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[8] + estimate[1]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*estimate[8] + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) +
                                                      estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[8] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) )  +
    I((data$E_muon>26.37) & (data$E_muon<=26.78))*( (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*estimate[8] + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) +
                                                      estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[8] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8])+
                                                      estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*estimate[7] + estimate[2]*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      estimate[1]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8]) +
                                                      estimate[2]*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) + ( 1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*estimate[7] +
                                                      (1-estimate[1]-estimate[2])*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*estimate[7] +
                                                      estimate[1]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) )  +
    I((data$E_muon>22.97) & (data$E_muon<=26.37))*( estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8])+
                                                      estimate[2]*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) + 
                                                      (1-estimate[1]-estimate[2])*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) +  estimate[1]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) )  +
    I((data$E_muon>22.56) & (data$E_muon<=22.97))*( estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8])+
                                                      estimate[2]*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) + estimate[1]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8])+
                                                      (1-estimate[1]-estimate[2])*estimate[3]*estimate[6]*estimate[8] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*estimate[7] + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*estimate[7] +
                                                      estimate[1]*(1-estimate[3]-estimate[4])*estimate[6]*estimate[8] + estimate[1]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[2]*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) +
                                                      estimate[2]*estimate[3]*(1-estimate[5]-estimate[6])*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5]*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*estimate[7] +
                                                      estimate[2]*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6])*estimate[8] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[5]*estimate[8]  )  +
    I((data$E_muon>19.22) & (data$E_muon<=22.56))*( (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*estimate[8] +
                                                      estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*estimate[8] + estimate[2]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) + 
                                                      estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) +       (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8])+
                                                      estimate[2]*estimate[3]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +  estimate[1]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5]*(1-estimate[7]-estimate[8]) +   (1-estimate[1]-estimate[2])*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) + estimate[1]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) +
                                                      (1-estimate[1]-estimate[2])*estimate[3]*estimate[6]*estimate[8] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*estimate[7] +  estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])* estimate[7] +
                                                      estimate[1]*(1-estimate[3]-estimate[4])*estimate[6]*estimate[8] + estimate[1]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[2]*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8])+
                                                      estimate[2]*estimate[3]*(1-estimate[5]-estimate[6])*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[5]*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*estimate[7] +
                                                      estimate[2]*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[1]*estimate[4]*(1-estimate[5]-estimate[6])*estimate[8] + (1-estimate[1]-estimate[2])*estimate[4]*estimate[5]*estimate[8]  ) +
    I((data$E_muon>18.81) & (data$E_muon<=19.22))*( estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*(1-estimate[7]-estimate[8]) + (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8])+
                                                      (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*estimate[7] + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*estimate[7] + estimate[1]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) + 
                                                      estimate[2]*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*estimate[7] + estimate[2]*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8])+
                                                      (1-estimate[1]-estimate[2])*estimate[4]*estimate[6]*estimate[8] + estimate[2]*(1-estimate[3]-estimate[4])*estimate[6]*estimate[8] + estimate[2]*estimate[4]*(1-estimate[5]-estimate[6])*estimate[8] +
                                                      estimate[2]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) ) +
    I((data$E_muon>15.41) & (data$E_muon<=18.81))*( estimate[1]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) + 
                                                      estimate[2]*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) + estimate[2]*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8])+ estimate[2]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) ) +
    I((data$E_muon>15) & (data$E_muon<=15.41))*( estimate[1]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) + 
                                                   estimate[2]*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) + estimate[2]*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8])+ estimate[2]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) +
                                                   estimate[1]*estimate[4]*estimate[6]*estimate[8] + estimate[2]*estimate[3]*estimate[6]*estimate[8] + estimate[2]*estimate[4]*estimate[5]*estimate[8] + estimate[2]*estimate[4]*estimate[6]*estimate[7] ) +
    I((data$E_muon>11.66) & (data$E_muon<=15))*( estimate[1]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) + estimate[2]*estimate[4]*estimate[6]*estimate[8] +
                                                   estimate[2]*estimate[4]*estimate[5]*(1-estimate[7]-estimate[8]) + estimate[2]*estimate[3]*estimate[6]*(1-estimate[7]-estimate[8])+ estimate[2]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8]) +
                                                   estimate[1]*estimate[4]*estimate[6]*estimate[8] + estimate[2]*estimate[3]*estimate[6]*estimate[8] + estimate[2]*estimate[4]*estimate[5]*estimate[8] + estimate[2]*estimate[4]*estimate[6]*estimate[7] ) +
    I((data$E_muon>11.25) & (data$E_muon<=11.66))*( estimate[2]*estimate[4]*estimate[6]*estimate[8] + estimate[2]*estimate[4]*estimate[6]*(1-estimate[7]-estimate[8])+
                                                      estimate[2]*estimate[4]*estimate[6]*estimate[7] ) +
    I(data$E_muon<=11.25)*10^-15
  
  theta = cbind(theta_il1, theta_il2, theta_il3, theta_il4)
  fun = function(x) 
  {
    # theta_ilj 
    theta_il1.s = x[1]/(sum(x))
    theta_il2.s = x[2]/(sum(x))
    theta_il3.s = x[3]/(sum(x))
    theta_il4.s = x[4]/(sum(x))
    
    c(theta_il1.s,theta_il2.s,theta_il3.s,theta_il4.s)
    
  }
  theta = as.data.frame(t(apply(theta, 1, fun)))
  
  
  # estimate[1] = p11, estimate[2] = p21, 1-estimate[1]-estimate[2] = p31
  # estimate[3] = p12, estimate[4] = p22, 1-estimate[3]-estimate[4] = p32
  # estimate[5] = p13, estimate[6] = p23, 1-estimate[5]-estimate[6] = p33
  # estimate[7] = p14, estimate[8] = p24, 1-estimate[7]-estimate[8] = p34
  # p = (p11, p21, p12, p22, p13, p23, p14, p24)
  
  Expectation_1 = ( (theta$theta_il1*(estimate[1]*fE.cond_sz.Z6.S1 + estimate[2]*fE.cond_sz.Z14.S1 + (1-estimate[1]-estimate[2])*fE.cond_sz.Z26.S1))/
                      (theta$theta_il1*(estimate[1]*fE.cond_sz.Z6.S1 + estimate[2]*fE.cond_sz.Z14.S1 + (1-estimate[1]-estimate[2])*fE.cond_sz.Z26.S1) + 
                         theta$theta_il2*(estimate[3]*fE.cond_sz.Z6.S2 + estimate[4]*fE.cond_sz.Z14.S2 + (1-estimate[3]-estimate[4])*fE.cond_sz.Z26.S2) + 
                         theta$theta_il3*(estimate[5]*fE.cond_sz.Z6.S3 + estimate[6]*fE.cond_sz.Z14.S3 + (1-estimate[5]-estimate[6])*fE.cond_sz.Z26.S3) + 
                         theta$theta_il4*(estimate[7]*fE.cond_sz.Z6.S4 + estimate[8]*fE.cond_sz.Z14.S4 + (1-estimate[7]-estimate[8])*fE.cond_sz.Z26.S4)) ) 
  
  Expectation_2 = ( (theta$theta_il2*(estimate[3]*fE.cond_sz.Z6.S2 + estimate[4]*fE.cond_sz.Z14.S2 + (1-estimate[3]-estimate[4])*fE.cond_sz.Z26.S2))/
                      (theta$theta_il1*(estimate[1]*fE.cond_sz.Z6.S1 + estimate[2]*fE.cond_sz.Z14.S1 + (1-estimate[1]-estimate[2])*fE.cond_sz.Z26.S1) + 
                         theta$theta_il2*(estimate[3]*fE.cond_sz.Z6.S2 + estimate[4]*fE.cond_sz.Z14.S2 + (1-estimate[3]-estimate[4])*fE.cond_sz.Z26.S2) + 
                         theta$theta_il3*(estimate[5]*fE.cond_sz.Z6.S3 + estimate[6]*fE.cond_sz.Z14.S3 + (1-estimate[5]-estimate[6])*fE.cond_sz.Z26.S3) + 
                         theta$theta_il4*(estimate[7]*fE.cond_sz.Z6.S4 + estimate[8]*fE.cond_sz.Z14.S4 + (1-estimate[7]-estimate[8])*fE.cond_sz.Z26.S4)) ) 
  
  Expectation_3 = ( (theta$theta_il3*(estimate[5]*fE.cond_sz.Z6.S3 + estimate[6]*fE.cond_sz.Z14.S3 + (1-estimate[5]-estimate[6])*fE.cond_sz.Z26.S3))/
                      (theta$theta_il1*(estimate[1]*fE.cond_sz.Z6.S1 + estimate[2]*fE.cond_sz.Z14.S1 + (1-estimate[1]-estimate[2])*fE.cond_sz.Z26.S1) + 
                         theta$theta_il2*(estimate[3]*fE.cond_sz.Z6.S2 + estimate[4]*fE.cond_sz.Z14.S2 + (1-estimate[3]-estimate[4])*fE.cond_sz.Z26.S2) + 
                         theta$theta_il3*(estimate[5]*fE.cond_sz.Z6.S3 + estimate[6]*fE.cond_sz.Z14.S3 + (1-estimate[5]-estimate[6])*fE.cond_sz.Z26.S3) + 
                         theta$theta_il4*(estimate[7]*fE.cond_sz.Z6.S4 + estimate[8]*fE.cond_sz.Z14.S4 + (1-estimate[7]-estimate[8])*fE.cond_sz.Z26.S4)) ) 
  out = as.data.frame(cbind(Expectation_1,Expectation_2,Expectation_3))
  return(out)
}



EM.depth.4cm = function(start, data, eps = 1e-5, max.iter = 100)
{
  library('nloptr')
  opts = list("algorithm"="NLOPT_LD_MMA",
              "xtol_rel"=1.0e-8)
  abs.e = 1
  estimate = start # p11, p21, p12, p22, p13, p23, p14, p24
  new.estimate = numeric(length(start))
  iter = 0
  
  # valore atteso della log-verosimiglianza per i dati completi: funzione da massimizzare nell'M-step
  function.maximize = function(p, data, expectation) # p = (p11, p21, p12, p22, p13, p23, p14, p24)
  {
    -sum( expectation$Expectation_1 * log( p[1]*fE.cond_sz.Z6.S1 + p[2]*fE.cond_sz.Z14.S1 + (1-p[1]-p[2])*fE.cond_sz.Z26.S1 ) +
            expectation$Expectation_2 * log( p[3]*fE.cond_sz.Z6.S2 + p[4]*fE.cond_sz.Z14.S2 + (1-p[3]-p[4])*fE.cond_sz.Z26.S2 ) +
            expectation$Expectation_3 * log( p[5]*fE.cond_sz.Z6.S3 + p[6]*fE.cond_sz.Z14.S3 + (1-p[5]-p[6])*fE.cond_sz.Z26.S3 ) +
            (1 - expectation$Expectation_1 - expectation$Expectation_2 - expectation$Expectation_3) * log( p[7]*fE.cond_sz.Z6.S4 + p[8]*fE.cond_sz.Z14.S4 + (1-p[7]-p[8])*fE.cond_sz.Z26.S4 )  
    )
  }
  
  # gradiente della funzione da massimizzare nell'M-step
  gradient = function(p, data, expectation)
  {
    c((-sum(expectation$Expectation_1*((fE.cond_sz.Z6.S1-fE.cond_sz.Z26.S1)/
                                         (p[1]*fE.cond_sz.Z6.S1+p[2]*fE.cond_sz.Z14.S1+(1-p[1]-p[2])*fE.cond_sz.Z26.S1)))),
      (-sum(expectation$Expectation_1*((fE.cond_sz.Z14.S1-fE.cond_sz.Z26.S1)/
                                         (p[1]*fE.cond_sz.Z6.S1+p[2]*fE.cond_sz.Z14.S1+(1-p[1]-p[2])*fE.cond_sz.Z26.S1)))),
      (-sum(expectation$Expectation_2*((fE.cond_sz.Z6.S2-fE.cond_sz.Z26.S2)/
                                         (p[3]*fE.cond_sz.Z6.S2+p[4]*fE.cond_sz.Z14.S2+(1-p[3]-p[4])*fE.cond_sz.Z26.S2)))),
      (-sum(expectation$Expectation_2*((fE.cond_sz.Z14.S2-fE.cond_sz.Z26.S2)/
                                         (p[3]*fE.cond_sz.Z6.S2+p[4]*fE.cond_sz.Z14.S2+(1-p[3]-p[4])*fE.cond_sz.Z26.S2)))),
      (-sum(expectation$Expectation_3*((fE.cond_sz.Z6.S3-fE.cond_sz.Z26.S3)/
                                         (p[5]*fE.cond_sz.Z6.S3+p[6]*fE.cond_sz.Z14.S3+(1-p[5]-p[6])*fE.cond_sz.Z26.S3)))),
      (-sum(expectation$Expectation_3*((fE.cond_sz.Z14.S3-fE.cond_sz.Z26.S3)/
                                         (p[5]*fE.cond_sz.Z6.S3+p[6]*fE.cond_sz.Z14.S3+(1-p[5]-p[6])*fE.cond_sz.Z26.S3)))),
      (-sum((1-expectation$Expectation_1-expectation$Expectation_2-expectation$Expectation_3)*((fE.cond_sz.Z6.S4-fE.cond_sz.Z26.S4)/
                                                                                                 (p[7]*fE.cond_sz.Z6.S4+p[8]*fE.cond_sz.Z14.S4+(1-p[7]-p[8])*fE.cond_sz.Z26.S4)))),
      (-sum((1-expectation$Expectation_1-expectation$Expectation_2-expectation$Expectation_3)*((fE.cond_sz.Z14.S4-fE.cond_sz.Z26.S4)/
                                                                                                 (p[7]*fE.cond_sz.Z6.S4+p[8]*fE.cond_sz.Z14.S4+(1-p[7]-p[8])*fE.cond_sz.Z26.S4))))
    )
  }
  
  # funzione che determina i vincoli 
  eval_g0 = function(p, data, expectation){
    c(p[1]+p[2]-1+1e-7, p[3]+p[4]-1+1e-7, p[5]+p[6]-1+1e-7, p[7]+p[8]-1+1e-7)
  }
  
  # jacobiano dei vincoli
  eval_jac_g0 = function(p, data, expectation){
    return( rbind( c( 1, 1, 0, 0, 0, 0, 0, 0),
                   c( 0, 0, 1, 1, 0, 0, 0, 0),
                   c( 0, 0, 0, 0, 1, 1, 0, 0),
                   c( 0, 0, 0, 0, 0, 0, 1, 1) ) )
  }

  while ((abs.e > eps) & (iter < max.iter)) 
  {
    iter = iter + 1
    
    # E-step
    out = Expectation(estimate = estimate, data = data) # funzione che calcola i valori attesi per l'E-step
    
    # M-step
    new.estimate = nloptr( x0=estimate, eval_f=function.maximize, eval_grad_f=gradient, eval_g_ineq = eval_g0, eval_jac_g_ineq = eval_jac_g0,
                           data=data, expectation = out,  opts=opts,
                           lb = c(rep(1e-7,8)), ub=c(rep(1-1e-7,8)) )$solution
    
    
    abs.e = max(abs(new.estimate - estimate))
    estimate = new.estimate
  }
  
  if (iter == max.iter) warning("Maximum number of iteration reached")
  p31 = 1-estimate[1]-estimate[2]
  p32 = 1-estimate[3]-estimate[4]
  p33 = 1-estimate[5]-estimate[6]
  p34 = 1-estimate[7]-estimate[8]
  estimate = c(estimate[1:2], p31, estimate[3:4], p32, estimate[5:6], p33, estimate[7:8], p34)
  list(estimate = estimate, iter = iter)
}







#-----------------------
# EM per i voxel 1,2,3,4 
#-----------------------

data = create.data(pixel.detector.muon = 1, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox1.2.3.4 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox1.2.3.4$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 5,6,7,8 
#-----------------------

data = create.data(pixel.detector.muon = 2, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox5.6.7.8 = EM.depth.4cm(start =  c(rep(1/3,8)), data=data, max.iter=100, eps = 0.001)
round(solution.vox5.6.7.8$estimate,3)[c(1,2,4,5,7,8,10,11)]


#-----------------------
# EM per i voxel 9,10,11,12 
#-----------------------
data = create.data(pixel.detector.muon = 3, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox9.10.11.12 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter = 100, eps = 0.001)
round(solution.vox9.10.11.12$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 13,14,15,16 
#-----------------------

data = create.data(pixel.detector.muon = 4, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox13.14.15.16 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=300, eps=0.001)
round(solution.vox13.14.15.16$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 17,18,19,20 
#-----------------------

data = create.data(pixel.detector.muon = 5, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox17.18.19.20 = EM.depth.4cm(start = c(rep(1/3,8)), data=data,max.iter = 100,eps=0.001)
round(solution.vox17.18.19.20$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 21,22,23,24 
#-----------------------

data = create.data(pixel.detector.muon = 6, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox21.22.23.24 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox21.22.23.24$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 25,26,27,28 
#-----------------------

data = create.data(pixel.detector.muon = 7, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox25.26.27.28 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox25.26.27.28$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 29,30,31,32
#-----------------------

cube3[c(29,30,31,32),] # voxel di cui voglio identificare il numero atomico
coord.riv1[8,] # pixel del rivelatore che corrisponde ai muoni che, entrando in questo pixel, si sono poi fermati nei voxel 29,30,31,32

data = create.data(pixel.detector.muon = 8, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox29.30.31.32 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox29.30.31.32$estimate,3)[c(1,2,4,5,7,8,10,11)]




#-----------------------
# EM per i voxel 33,34,35,36
#-----------------------

data = create.data(pixel.detector.muon = 9, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox33.34.35.36 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=200, eps=0.001)
round(solution.vox33.34.35.36$estimate,3)[c(1,2,4,5,7,8,10,11)]


#-----------------------
# EM per i voxel 37,38,39,40
#-----------------------

data = create.data(pixel.detector.muon = 10, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox37.38.39.40 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=200, eps=0.001)
round(solution.vox37.38.39.40$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 41,42,43,44
#-----------------------

cube3[c(41,42,43,44),] # voxel di cui voglio identificare il numero atomico
coord.riv1[11,] # pixel del rivelatore che corrisponde ai muoni che, entrando in questo pixel, si sono poi fermati nei voxel 41,42,43,44

data = create.data(pixel.detector.muon = 11, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox41.42.43.44 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=200, eps=0.001)
round(solution.vox41.42.43.44$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 45,46,47,48
#-----------------------

data = create.data(pixel.detector.muon = 12, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox45.46.47.48= EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=200, eps=0.001)
round(solution.vox45.46.47.48$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 49,50,51,52
#-----------------------

data = create.data(pixel.detector.muon = 13, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox49.50.51.52 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=300,eps=0.001)
round(solution.vox49.50.51.52$estimate,3)[c(1,2,4,5,7,8,10,11)]




#-----------------------
# EM per i voxel 53,54,55,56
#-----------------------

data = create.data(pixel.detector.muon = 14, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox53.54.55.56 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, eps=0.001, max.iter=100)
round(solution.vox53.54.55.56$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 57,58,59,60
#-----------------------

data = create.data(pixel.detector.muon = 15, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox57.58.59.60 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, eps=0.001,max.iter=200)
round(solution.vox57.58.59.60$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 61,62,63,64
#-----------------------

data = create.data(pixel.detector.muon = 16, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz.Z6.S2 = fE.cond_sz.Z6.S3 = fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0))
fE.cond_sz.Z14.S1 = fE.cond_sz.Z14.S2 = fE.cond_sz.Z14.S3 = fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0))
fE.cond_sz.Z26.S1 = fE.cond_sz.Z26.S2 = fE.cond_sz.Z26.S3 = fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1))

solution.vox61.62.63.64 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox61.62.63.64$estimate,3)[c(1,2,4,5,7,8,10,11)]











