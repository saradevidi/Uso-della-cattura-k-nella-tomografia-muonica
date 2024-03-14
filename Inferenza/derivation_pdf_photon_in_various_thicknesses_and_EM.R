
# carico le funzioni 
source("create_geometry.R")
source("muon_energy_loss.R")
source("pdf_photon_energy.R")
source("absorption_length.R")
source("function_generatekCapture.R")



# simulo i dati -> muoni che incidono su un volume composto com'e' stato stimato dall'EM precedentemente applicato

cube1 = coord.volume(coord.ang = c(24,0,0), alt = 4, lato1 = 4, lato2 = 4)
cube2 = coord.volume.Extreme(coord.volume.matrix = cube1)
cube3 = coord.volume.Z(coord.volume.Extreme.matrix = cube2, atoms = c(6,14,26), seed_ = 123)
cube3$Z[4]=26 
cube3$Z[20]= 26
cube3$Z[22]=26
cube3$Z[31]=26 
cube3$Z[44]= 14
cube3$Z[46]=26
cube3$Z[63]= 26
cube3$Z[34]=26
cube3$Z[39]=26

coord.riv1 = coord.riv(x.coord = 3.5, coord.volume = cube1, lato2 = 4, alt = 4) 

densitaZ = matrix(c(6, 14, 26, 2.25, 2.33, 7.87), ncol=2) 

DeltaE.matrix = loss.energy(Zdensity = densitaZ, Y = 1)

kcapturePROB.matrix = cbind(DeltaE.matrix$Z,c(0.99, 0.99, 0.99))
colnames(kcapturePROB.matrix) = c("Z", "prob_cattura")
kcapturePROB.matrix = as.data.frame(kcapturePROB.matrix)  

Survival.lambda = function(spessore, densita, lambda) 
{
  exp( -(spessore * densita) / lambda)
}

Simulation = generate.kCapture(n.muon = 8000000, E_mu = c(2,40), coord.material = cube3, 
                               coord.detector = coord.riv1, DeltaE.matrix = DeltaE.matrix, 
                               kcapturePROB.matrix = kcapturePROB.matrix,
                               lato1 = 4, x = 0, seed_ = 100000)

Simulation = Simulation[(!is.na(Simulation$E_fot1) | !is.na(Simulation$E_fot2) | !is.na(Simulation$E_fot3)| !is.na(Simulation$E_fot4)) ,] # prendo solo i muoni per i quali si è visto almeno un fotone

save(Simulation, file="data_pdf_fotoni_nei_vari_spessori.RData")

# carico i dati ottenuti
load("data_pdf_fotoni_nei_vari_spessori.RData")




# suddivido le energie dei fotoni per numero atomico e spessore di provenienza
#-----------------------------------------------------------------------------

# prendo i fotoni usciti dal carbonio nel primo spessore in x
fot.Z6.first = c(Simulation[Simulation$vox_real %in% c(13,21,37,57,61,49),9],
                 Simulation[Simulation$vox_real %in% c(13,21,37,57,61,49),10],
                 Simulation[Simulation$vox_real %in% c(13,21,37,57,61,49),11],
                 Simulation[Simulation$vox_real %in% c(13,21,37,57,61,49),12])
fot.Z6.first = fot.Z6.first[!is.na(fot.Z6.first)]
hist(fot.Z6.first,nclass = 100,prob=T)

# prendo i fotoni usciti dal silicio nel primo spessore in x
fot.Z14.first = c(Simulation[Simulation$vox_real %in% 41,9],
                  Simulation[Simulation$vox_real %in% 41,10],
                  Simulation[Simulation$vox_real %in% 41,11],
                  Simulation[Simulation$vox_real %in% 41,12])
fot.Z14.first = fot.Z14.first[!is.na(fot.Z14.first)]
hist(fot.Z14.first,nclass = 100,prob=T)

# prendo i fotoni usciti dal ferro nel primo spessore in x
fot.Z26.first = c(Simulation[Simulation$vox_real %in% c(1,5,9,17,25,29,33,45,53),9],
                  Simulation[Simulation$vox_real %in%  c(1,5,9,17,25,29,33,45,53),10],
                  Simulation[Simulation$vox_real %in%  c(1,5,9,17,25,29,33,45,53),11],
                  Simulation[Simulation$vox_real %in%  c(1,5,9,17,25,29,33,45,53),12])
fot.Z26.first = fot.Z26.first[!is.na(fot.Z26.first)]
hist(fot.Z26.first,nclass = 100,prob=T)





# prendo i fotoni usciti dal carbonio nel secondo spessore in x
fot.Z6.second = c(Simulation[Simulation$vox_real %in% c(10,54),9],
                 Simulation[Simulation$vox_real %in% c(10,54),10],
                 Simulation[Simulation$vox_real %in% c(10,54),11],
                 Simulation[Simulation$vox_real %in% c(10,54),12])
fot.Z6.second = fot.Z6.second[!is.na(fot.Z6.second)]
hist(fot.Z6.second,nclass = 100,prob=T)

# prendo i fotoni usciti dal silicio nel secondo spessore in x
fot.Z14.second = c(Simulation[Simulation$vox_real %in% c(6,14,26,30,50),9],
                  Simulation[Simulation$vox_real %in%  c(6,14,26,30,50),10],
                  Simulation[Simulation$vox_real %in%  c(6,14,26,30,50),11],
                  Simulation[Simulation$vox_real %in%  c(6,14,26,30,50),12])
fot.Z14.second = fot.Z14.second[!is.na(fot.Z14.second)]
hist(fot.Z14.second,nclass = 100,prob=T)

# prendo i fotoni usciti dal ferro nel secondo spessore in x
fot.Z26.second = c(Simulation[Simulation$vox_real %in% c(2,18,22,34,38,42,46,58,62),9],
                  Simulation[Simulation$vox_real %in%  c(2,18,22,34,38,42,46,58,62),10],
                  Simulation[Simulation$vox_real %in%  c(2,18,22,34,38,42,46,58,62),11],
                  Simulation[Simulation$vox_real %in%  c(2,18,22,34,38,42,46,58,62),12])
fot.Z26.second = fot.Z26.second[!is.na(fot.Z26.second)]
hist(fot.Z26.second,nclass = 100,prob=T)




# prendo i fotoni usciti dal carbonio nel terzo spessore in x
fot.Z6.third = c(Simulation[Simulation$vox_real %in% c(19,27,51,59),9],
                 Simulation[Simulation$vox_real %in% c(19,27,51,59),10],
                 Simulation[Simulation$vox_real %in% c(19,27,51,59),11],
                 Simulation[Simulation$vox_real %in% c(19,27,51,59),12])
fot.Z6.third = fot.Z6.third[!is.na(fot.Z6.third)]
hist(fot.Z6.third,nclass = 100,prob=T)

# prendo i fotoni usciti dal silicio nel terzo spessore in x
fot.Z14.third = c(Simulation[Simulation$vox_real %in% c(7,11,55),9],
                  Simulation[Simulation$vox_real %in% c(7,11,55),10],
                  Simulation[Simulation$vox_real %in% c(7,11,55),11],
                  Simulation[Simulation$vox_real %in% c(7,11,55),12])
fot.Z14.third = fot.Z14.third[!is.na(fot.Z14.third)]
hist(fot.Z14.third,nclass = 100,prob=T)

# prendo i fotoni usciti dal ferro nel terzo spessore in x
fot.Z26.third = c(Simulation[Simulation$vox_real %in% c(3,15,23,31,35,39,43,47,63),9],
                  Simulation[Simulation$vox_real %in%  c(3,15,23,31,35,39,43,47,63),10],
                  Simulation[Simulation$vox_real %in%  c(3,15,23,31,35,39,43,47,63),11],
                  Simulation[Simulation$vox_real %in%  c(3,15,23,31,35,39,43,47,63),12])
fot.Z26.third = fot.Z26.third[!is.na(fot.Z26.third)]
hist(fot.Z26.third,nclass = 100,prob=T)





# prendo i fotoni usciti dal carbonio nel quarto spessore in x
fot.Z6.fourth = c(Simulation[Simulation$vox_real %in% c(16,40,52,56),9],
                 Simulation[Simulation$vox_real %in% c(16,40,52,56),10],
                 Simulation[Simulation$vox_real %in% c(16,40,52,56),11],
                 Simulation[Simulation$vox_real %in% c(16,40,52,56),12])
fot.Z6.fourth = fot.Z6.fourth[!is.na(fot.Z6.fourth)]
hist(fot.Z6.fourth,nclass = 100,prob=T)

# prendo i fotoni usciti dal silicio nel quarto spessore in x
fot.Z14.fourth = c(Simulation[Simulation$vox_real %in% c(8,12,24,28,36,44,48,60),9],
                  Simulation[Simulation$vox_real %in% c(8,12,24,28,36,44,48,60),10],
                  Simulation[Simulation$vox_real %in% c(8,12,24,28,36,44,48,60),11],
                  Simulation[Simulation$vox_real %in% c(8,12,24,28,36,44,48,60),12])
fot.Z14.fourth = fot.Z14.fourth[!is.na(fot.Z14.fourth)]
hist(fot.Z14.fourth,nclass = 100,prob=T)

# prendo i fotoni usciti dal ferro nel quarto spessore in x
fot.Z26.fourth = c(Simulation[Simulation$vox_real %in% c(4,32,20,64),9],
                  Simulation[Simulation$vox_real %in%  c(4,32,20,64),10],
                  Simulation[Simulation$vox_real %in%  c(4,32,20,64),11],
                  Simulation[Simulation$vox_real %in%  c(4,32,20,64),12])
fot.Z26.fourth = fot.Z26.fourth[!is.na(fot.Z26.fourth)]
hist(fot.Z26.fourth,nclass = 100,prob=T)





# stimo i parametri delle densità f(Ei|sij=1,ztij=1)
#----------------------------------------------------

estimate.Z26.first = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,2,0.1,6,0.1,9.5,0.1),data = fot.Z26.first, eps=0.1)
estimate.Z14.first = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,1.5,0.1,5,0.1,7,0.1),data = fot.Z14.first, eps=0.1)
estimate.Z6.first = EM.mixture(start= c(1/3,1/3,1/3,0.1,0.1,1,0.1,3,0.1),data = fot.Z6.first, eps=0.1)

estimate.Z26.second = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,2,0.1,6,0.1,9.5,0.1),data = fot.Z26.second, eps=0.1)
estimate.Z14.second = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,1.5,0.1,5,0.1,7,0.1),data = fot.Z14.second, eps=0.1)
estimate.Z6.second = EM.mixture(start= c(1/3,1/3,1/3,0.1,0.1,1,0.1,3,0.1),data = fot.Z6.second, eps=0.1)

estimate.Z26.third = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,2,0.1,6,0.1,9.5,0.1),data = fot.Z26.third, eps=0.01)
estimate.Z14.third = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,1.5,0.1,5,0.1,7,0.1),data = fot.Z14.third, eps=0.01)
estimate.Z6.third = EM.mixture(start= c(1/3,1/3,1/3,0.1,0.1,1,0.1,3,0.1),data = fot.Z6.third, eps=0.01)

estimate.Z26.fourth = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,2,0.1,6,0.1,9.5,0.1),data = fot.Z26.fourth, eps=0.001)
estimate.Z14.fourth = EM.mixture(start= c(1/4,1/4,1/4,1/4,0.1,0.1,1.5,0.1,5,0.1,7,0.1),data = fot.Z14.fourth, eps=0.001)
estimate.Z6.fourth = EM.mixture(start= c(1/3,1/3,1/3,0.1,0.1,1,0.1,3,0.1),data = fot.Z6.fourth, eps=0.001)



estimate.Z26 = matrix(c(estimate.Z26.first$estimate,estimate.Z26.second$estimate,
                        estimate.Z26.third$estimate,estimate.Z26.fourth$estimate),ncol=4)
estimate.Z14 = matrix(c(estimate.Z14.first$estimate,estimate.Z14.second$estimate,
                        estimate.Z14.third$estimate,estimate.Z14.fourth$estimate),ncol=4)
estimate.Z6 = matrix(c(estimate.Z6.first$estimate,estimate.Z6.second$estimate,
                       estimate.Z6.third$estimate,estimate.Z6.fourth$estimate),ncol=4)





weight.calculation = function(z, s)
{ # calcola i pesi k1.s, k2.s, k3.s, k4.s
  nmax = 4 # massimo dei picchi per le distribuzioni
  out = c(NA, nmax) 
  
  # pesi assunti noti delle distribuzioni mistura delle energie
  k1=c(estimate.Z6[1,s] , estimate.Z14[1,s] ,  estimate.Z26[1,s] )
  k2=c(estimate.Z6[2,s] , estimate.Z14[2,s] ,  estimate.Z26[2,s] )
  k3=c(estimate.Z6[3,s] , estimate.Z14[3,s] , estimate.Z26[3,s] )
  k4=c(0,   estimate.Z14[4,s] , estimate.Z26[4,s] )
  
  k1.s=  as.numeric(t(z)%*%k1)
  k2.s= as.numeric(t(z)%*%k2)
  k3.s= as.numeric(t(z)%*%k3)
  k4.s= as.numeric(t(z)%*%k4)
  
  out = c(k1.s, k2.s, k3.s, k4.s)
  out
}


fE.cond_sz = function(E, s, z) 
{ # distribuzione delle energie condizionata allo spessore e al numero atomico
  delta = c(estimate.Z6[4,s] ,estimate.Z14[5,s] ,estimate.Z26[5,s] ) 
  beta = c(estimate.Z6[5,s] ,estimate.Z14[6,s] ,estimate.Z26[6,s] )
  mu1=c(estimate.Z6[6,s] ,estimate.Z14[7,s] ,estimate.Z26[7,s] ) # primo picco
  sigma1=c(estimate.Z6[7,s] ,estimate.Z14[8,s] ,estimate.Z26[8,s] ) # primo picco
  mu2=c(estimate.Z6[8,s]  ,estimate.Z14[9,s] ,estimate.Z26[9,s] ) # secondo picco
  sigma2=c(estimate.Z6[9,s],estimate.Z14[10,s] ,estimate.Z26[10,s] ) # secondo picco
  mu3=c(0,estimate.Z14[11,s] ,estimate.Z26[11,s] ) # terzo picco
  sigma3=c(0,estimate.Z14[12,s],estimate.Z26[12,s]) # terzo picco
  
  Kj.s = weight.calculation(z=z, s=s)
  
  out= as.numeric(Kj.s[1] * ( dgamma(E,shape=t(z)%*%delta,rate = t(z)%*%beta)*I(E>0) ) + 
                    Kj.s[2] * ( (1/pnorm(q=(-as.numeric(t(z)%*%mu1)/as.numeric(t(z)%*%sigma1)),lower.tail = FALSE))*dnorm(E, mean = as.numeric(t(z)%*%mu1), sd=as.numeric(t(z)%*%sigma1))*I(E>0) ) +
                    Kj.s[3] * ( (1/pnorm(q=(-as.numeric(t(z)%*%mu2)/as.numeric(t(z)%*%sigma2)),lower.tail = FALSE))*dnorm(E, mean = as.numeric(t(z)%*%mu2), sd=as.numeric(t(z)%*%sigma2))*I(E>0) ) +
                    Kj.s[4] * replace( (1/pnorm(q=(-as.numeric(t(z)%*%mu3)/as.numeric(t(z)%*%sigma3)),lower.tail = FALSE))*dnorm(E, mean = as.numeric(t(z)%*%mu3), sd=as.numeric(t(z)%*%sigma3))*I(E>0), 
                                       is.na((1/pnorm(q=(-as.numeric(t(z)%*%mu3)/as.numeric(t(z)%*%sigma3)),lower.tail = FALSE))*dnorm(E, mean = as.numeric(t(z)%*%mu3), sd=as.numeric(t(z)%*%sigma3))*I(E>0)),
                                       0) ) 
  
  
  out
}









#############################################
# EM applicato nuovamente a tutto il volume #
#############################################

# carico i dati del volume incognito
load("data_for_inference_only_photon.RData")


# definizione del volume incognito
source("create_geometry.R")
cube1 = coord.volume(coord.ang = c(24,0,0), alt = 4, lato1 = 4, lato2 = 4)
cube2 = coord.volume.Extreme(coord.volume.matrix = cube1)
cube3 = coord.volume.Z(coord.volume.Extreme.matrix = cube2, atoms = c(6,14,26), seed_ = 123)
cube3$Z[4]=26 
coord.riv1 = coord.riv(x.coord = 3.5, coord.volume = cube1, lato2 = 4, alt = 4) 




#-----------------------
# EM per i voxel 1,2,3,4 
#-----------------------

data = create.data(pixel.detector.muon = 1, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox1.2.3.4 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox1.2.3.4$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 5,6,7,8 
#-----------------------

data = create.data(pixel.detector.muon = 2, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox5.6.7.8 = EM.depth.4cm(start =  c(rep(1/3,8)), data=data, max.iter=100, eps = 0.001)
round(solution.vox5.6.7.8$estimate,3)[c(1,2,4,5,7,8,10,11)]


#-----------------------
# EM per i voxel 9,10,11,12 
#-----------------------
data = create.data(pixel.detector.muon = 3, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox9.10.11.12 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter = 100, eps = 0.001)
round(solution.vox9.10.11.12$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 13,14,15,16 
#-----------------------

data = create.data(pixel.detector.muon = 4, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox13.14.15.16 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=300, eps=0.001)
round(solution.vox13.14.15.16$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 17,18,19,20 
#-----------------------

data = create.data(pixel.detector.muon = 5, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox17.18.19.20 = EM.depth.4cm(start = c(rep(1/3,8)), data=data,max.iter = 100,eps=0.001)
round(solution.vox17.18.19.20$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 21,22,23,24 
#-----------------------

data = create.data(pixel.detector.muon = 6, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox21.22.23.24 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox21.22.23.24$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 25,26,27,28 
#-----------------------

data = create.data(pixel.detector.muon = 7, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox25.26.27.28 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox25.26.27.28$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 29,30,31,32
#-----------------------

cube3[c(29,30,31,32),] # voxel di cui voglio identificare il numero atomico
coord.riv1[8,] # pixel del rivelatore che corrisponde ai muoni che, entrando in questo pixel, si sono poi fermati nei voxel 29,30,31,32

data = create.data(pixel.detector.muon = 8, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox29.30.31.32 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox29.30.31.32$estimate,3)[c(1,2,4,5,7,8,10,11)]




#-----------------------
# EM per i voxel 33,34,35,36
#-----------------------

data = create.data(pixel.detector.muon = 9, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox33.34.35.36 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=200, eps=0.001)
round(solution.vox33.34.35.36$estimate,3)[c(1,2,4,5,7,8,10,11)]


#-----------------------
# EM per i voxel 37,38,39,40
#-----------------------

data = create.data(pixel.detector.muon = 10, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox37.38.39.40 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=200, eps=0.001)
round(solution.vox37.38.39.40$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 41,42,43,44
#-----------------------

cube3[c(41,42,43,44),] # voxel di cui voglio identificare il numero atomico
coord.riv1[11,] # pixel del rivelatore che corrisponde ai muoni che, entrando in questo pixel, si sono poi fermati nei voxel 41,42,43,44

data = create.data(pixel.detector.muon = 11, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox41.42.43.44 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=200, eps=0.001)
round(solution.vox41.42.43.44$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 45,46,47,48
#-----------------------

data = create.data(pixel.detector.muon = 12, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox45.46.47.48= EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=200, eps=0.001)
round(solution.vox45.46.47.48$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 49,50,51,52
#-----------------------

data = create.data(pixel.detector.muon = 13, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox49.50.51.52 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=300,eps=0.001)
round(solution.vox49.50.51.52$estimate,3)[c(1,2,4,5,7,8,10,11)]




#-----------------------
# EM per i voxel 53,54,55,56
#-----------------------

data = create.data(pixel.detector.muon = 14, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox53.54.55.56 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, eps=0.001, max.iter=100)
round(solution.vox53.54.55.56$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 57,58,59,60
#-----------------------

data = create.data(pixel.detector.muon = 15, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox57.58.59.60 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, eps=0.001,max.iter=200)
round(solution.vox57.58.59.60$estimate,3)[c(1,2,4,5,7,8,10,11)]



#-----------------------
# EM per i voxel 61,62,63,64
#-----------------------

data = create.data(pixel.detector.muon = 16, Simulation = Simulation)

fE.cond_sz.Z6.S1 = fE.cond_sz(data$E,z=c(1,0,0),s=1)
fE.cond_sz.Z6.S2 = fE.cond_sz(data$E,z=c(1,0,0),s=2)
fE.cond_sz.Z6.S3 = fE.cond_sz(data$E,z=c(1,0,0),s=3)
fE.cond_sz.Z6.S4 = fE.cond_sz(data$E,z=c(1,0,0),s=4)
fE.cond_sz.Z14.S1 =fE.cond_sz(data$E,z=c(0,1,0), s=1)
fE.cond_sz.Z14.S2 = fE.cond_sz(data$E,z=c(0,1,0), s=2)
fE.cond_sz.Z14.S3 = fE.cond_sz(data$E,z=c(0,1,0), s=3)
fE.cond_sz.Z14.S4 = fE.cond_sz(data$E,z=c(0,1,0),s=4)
fE.cond_sz.Z26.S1 =fE.cond_sz(data$E,z=c(0,0,1),s=1)
fE.cond_sz.Z26.S2 =fE.cond_sz(data$E,z=c(0,0,1),s=2)
fE.cond_sz.Z26.S3 = fE.cond_sz(data$E,z=c(0,0,1),s=3)
fE.cond_sz.Z26.S4 = fE.cond_sz(data$E,z=c(0,0,1),s=4)

solution.vox61.62.63.64 = EM.depth.4cm(start = c(rep(1/3,8)), data=data, max.iter=100,eps=0.001)
round(solution.vox61.62.63.64$estimate,3)[c(1,2,4,5,7,8,10,11)]







