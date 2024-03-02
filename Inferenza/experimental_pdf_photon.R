
# faccio incidere muoni su uno spessore di 1cm di materiale di un dato numero atomico.


# carico le funzioni 
source("create_geometry.R")
source("muon_energy_loss.R")
source("pdf_photon_energy.R")
source("absorption_length.R")
source("function_generatekCapture.R")



# FERRO (Z=26)
#_________________________


cube1 = coord.volume(coord.ang = c(24,0,0), alt = 10, lato1 = 1, lato2 = 10)
cube2 = coord.volume.Extreme(coord.volume.matrix = cube1)
cube3 = coord.volume.Z(coord.volume.Extreme.matrix = cube2, atoms = c(26,26,26), seed_ = 123)
coord.riv1 = coord.riv(x.coord = 3.5, coord.volume = cube1, lato2 = 1, alt = 10) 
densitaZ = matrix(c(6, 14, 26, 2.25, 2.33, 7.87), ncol=2) 
DeltaE.matrix = loss.energy(Zdensity = densitaZ, Y = 1)
kcapturePROB.matrix = cbind(DeltaE.matrix$Z,c(0.99, 0.99, 0.99))
colnames(kcapturePROB.matrix) = c("Z", "prob_cattura")
kcapturePROB.matrix = as.data.frame(kcapturePROB.matrix)  
Survival.lambda = function(spessore, densita, lambda) 
{
  exp( -(spessore * densita) / lambda)
}

Simulation = generate.kCapture(n.muon = 10000, E_mu = c(1,12), coord.material = cube3, 
                               coord.detector = coord.riv1, DeltaE.matrix = DeltaE.matrix, 
                               kcapturePROB.matrix = kcapturePROB.matrix,
                               lato1 = 10, x = 0, seed_ = 8888888)

Simulation.photon = Simulation[(!is.na(Simulation$E_fot1)) | (!is.na(Simulation$E_fot2)) | (!is.na(Simulation$E_fot3)) | (!is.na(Simulation$E_fot4)),]

E.tot.26 = c(Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,1], Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,2],
          Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,3], Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,4])
E.tot.26 = E.tot.26[!is.na(E.tot.26)]


library(ggplot2)
library(gridExtra)

E.Z26 = c(seq(10^-10, 10, length=length(E.tot.26)))
out.Z26 = pdf.Z26(E.Z26)
 ggplot(data.frame(x1 = E.Z26, x2 = E.tot.26, y1 = out.Z26)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "red")+ # funzione di densità teorica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 3.2),breaks = c(0,0.5,1,1.5,2,2.5,3))+
  theme_minimal()+
  labs(title="Z = 26"
   #,subtitle = "Distribuzione empirica osservata delle energie da 1cm di spessore"
   )+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))
 

 
 # CARBONIO (Z=6)
 #_________________________
 
 cube1 = coord.volume(coord.ang = c(24,0,0), alt = 10, lato1 = 1, lato2 = 10)
 cube2 = coord.volume.Extreme(coord.volume.matrix = cube1)
 cube3 = coord.volume.Z(coord.volume.Extreme.matrix = cube2, atoms = c(6,6,6), seed_ = 123)
 coord.riv1 = coord.riv(x.coord = 3.5, coord.volume = cube1, lato2 = 1, alt = 10) 
 densitaZ = matrix(c(6, 14, 26, 2.25, 2.33, 7.87), ncol=2) 
 DeltaE.matrix = loss.energy(Zdensity = densitaZ, Y = 1)
 kcapturePROB.matrix = cbind(DeltaE.matrix$Z,c(0.99, 0.99, 0.99))
 colnames(kcapturePROB.matrix) = c("Z", "prob_cattura")
 kcapturePROB.matrix = as.data.frame(kcapturePROB.matrix)  
 Survival.lambda = function(spessore, densita, lambda) 
 {
   exp( -(spessore * densita) / lambda)
 }
 
 Simulation = generate.kCapture(n.muon = 10000, E_mu = c(1,5), coord.material = cube3, 
                                coord.detector = coord.riv1, DeltaE.matrix = DeltaE.matrix, 
                                kcapturePROB.matrix = kcapturePROB.matrix,
                                lato1 = 10, x = 0, seed_ = 111)
 

 Simulation.photon = Simulation[(!is.na(Simulation$E_fot1)) | (!is.na(Simulation$E_fot2)) | (!is.na(Simulation$E_fot3)) | (!is.na(Simulation$E_fot4)),]
 

 E.tot.6 = c(Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,1], Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,2],
           Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,3], Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,4])
 E.tot.6 = E.tot.6[!is.na(E.tot.6)]

 
 library(ggplot2)
 library(gridExtra)
 
 E.Z6 = c(seq(10^-10, 10, length=length(E.tot.6)))
 out.Z6 = pdf.Z6(E.Z6)
 ggplot(data.frame(x1 = E.Z6, x2 = E.tot.6, y1 = out.Z6)) + 
   geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
   geom_line(aes(y = y1, x = x1), col= "red")+ # funzione di densità teorica
   xlab(expression(paste("E", gamma, "(MeV)")))+
   ylab("Densità")+
   scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
   scale_y_continuous(limits = c(0, 3.2),breaks = c(0,0.5,1,1.5,2,2.5,3))+
   theme_minimal()+
   labs(title="Z = 6"
   #, subtitle = "Distribuzione empirica osservata delle energie da 1cm di spessore"
   )+
   theme(
     axis.title.x = element_text(size = 20),
     axis.text.x = element_text(size = 20),
     axis.title.y = element_text(size = 20),
     axis.text.y = element_text(size = 20),
     title =  element_text(size = 20))
 

 # SILICIO (Z=14)
 #_________________________
 
 cube1 = coord.volume(coord.ang = c(24,0,0), alt = 10, lato1 = 1, lato2 = 10)
 cube2 = coord.volume.Extreme(coord.volume.matrix = cube1)
 cube3 = coord.volume.Z(coord.volume.Extreme.matrix = cube2, atoms = c(14,14,14), seed_ = 123)
 coord.riv1 = coord.riv(x.coord = 3.5, coord.volume = cube1, lato2 = 1, alt = 10) 
 densitaZ = matrix(c(6, 14, 26, 2.25, 2.33, 7.87), ncol=2) 
 DeltaE.matrix = loss.energy(Zdensity = densitaZ, Y = 1)
 kcapturePROB.matrix = cbind(DeltaE.matrix$Z,c(0.99, 0.99, 0.99))
 colnames(kcapturePROB.matrix) = c("Z", "prob_cattura")
 kcapturePROB.matrix = as.data.frame(kcapturePROB.matrix)  
 Survival.lambda = function(spessore, densita, lambda) 
 {
   exp( -(spessore * densita) / lambda)
 }
 
 Simulation = generate.kCapture(n.muon = 10000, E_mu = c(1,4), coord.material = cube3, 
                                coord.detector = coord.riv1, DeltaE.matrix = DeltaE.matrix, 
                                kcapturePROB.matrix = kcapturePROB.matrix,
                                lato1 = 10, x = 0, seed_ = 222)
 

 Simulation.photon = Simulation[(!is.na(Simulation$E_fot1)) | (!is.na(Simulation$E_fot2)) | (!is.na(Simulation$E_fot3)) | (!is.na(Simulation$E_fot4)),]

 E.tot.14 = c(Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,1], Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,2],
           Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,3], Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,4])
 E.tot.14 = E.tot.14[!is.na(E.tot.14)]
 
 library(ggplot2)
 library(gridExtra)
 
 E.Z14 = c(seq(10^-10, 10, length=length(E.tot.14)))
 out.Z14 = pdf.Z14(E.Z14)
 ggplot(data.frame(x1 = E.Z14, x2 = E.tot.14, y1 = out.Z14)) + 
   geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
   geom_line(aes(y = y1, x = x1), col= "red")+ # funzione di densità teorica
   xlab(expression(paste("E", gamma, "(MeV)")))+
   ylab("Densità")+
   scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
   scale_y_continuous(limits = c(0, 3.2),breaks = c(0,0.5,1,1.5,2,2.5,3))+
   theme_minimal()+
   labs(title="Z = 14"
        #, subtitle = "Distribuzione empirica osservata delle energie da 1cm di spessore"
        )+
   theme(
     axis.title.x = element_text(size = 20),
     axis.text.x = element_text(size = 20),
     axis.title.y = element_text(size = 20),
     axis.text.y = element_text(size = 20),
     title =  element_text(size = 20))
 
 
 
 
 
 
 