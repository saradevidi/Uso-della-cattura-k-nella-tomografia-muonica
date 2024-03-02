# carico i dati
load("data_for_inference_only_photon.RData")

# quantità di fotoni rivelati in totale
vector.energy.tot = c(Simulation$E_fot1,Simulation$E_fot2,Simulation$E_fot3,Simulation$E_fot4) 
length(vector.energy.tot[!is.na(vector.energy.tot)])
E.tot = vector.energy.tot[!is.na(vector.energy.tot)]

#******************************************************
#* Distribuzione delle energie complessivamente rivelate
#******************************************************

library(ggplot2)
library(gridExtra)

ggplot(data.frame(x2 = E.tot)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.01)+ # densità empirica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))






#*********************************
#* Composizione atomica del volume
#*********************************
library("lattice")


# carico le funzioni 
source("create_geometry.R")
source("muon_energy_loss.R")
source("pdf_photon_energy.R")
source("absorption_length.R")
source("function_generatekCapture.R")



#**
# creo il volume incognito con cui ho ottenuto i dati
#**

cube1 = coord.volume(coord.ang = c(24,0,0), alt = 4, lato1 = 4, lato2 = 4)
cube2 = coord.volume.Extreme(coord.volume.matrix = cube1)
cube3 = coord.volume.Z(coord.volume.Extreme.matrix = cube2, atoms = c(6,14,26), seed_ = 123)
cube3$Z[4]=26 
coord.riv1 = coord.riv(x.coord = 3.5, coord.volume = cube1, lato2 = 4, alt = 4) 


coord.x = c("24.5","25.5","26.5","27.5")
coord.y = c("0.5","1.5","2.5","3.5")


# LEGEND:
# dark green -> Z = 6
# light green -> Z = 14
# white -> Z = 26


square1 = rbind(cube3$Z[1:4],
                cube3$Z[5:8],
                cube3$Z[9:12],
                cube3$Z[13:16])
colnames(square1) <- coord.x                    
rownames(square1) <- coord.y  
square1.graph = levelplot(square1, main=list("z = 0.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square1.graph


square2 = rbind(cube3$Z[17:20],
                cube3$Z[21:24],
                cube3$Z[25:28],
                cube3$Z[29:32])
colnames(square2) <- coord.x                    
rownames(square2) <- coord.y  
square2.graph = levelplot(square2, main=list("z = 1.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square2.graph


square3 = rbind(cube3$Z[33:36],
                cube3$Z[37:40],
                cube3$Z[41:44],
                cube3$Z[45:48])
colnames(square3) <- coord.x                    
rownames(square3) <- coord.y  
square3.graph = levelplot(square3, main=list("z = 2.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square3.graph


square4 = rbind(cube3$Z[49:52],
                cube3$Z[53:56],
                cube3$Z[57:60],
                cube3$Z[61:64])
colnames(square4) <- coord.x                    
rownames(square4) <- coord.y  
square4.graph = levelplot(square4, main=list("z = 3.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square4.graph






#***************************************
# PDF teoriche ricavate sperimentalmente
#***************************************
#*
source("parameter_mixture_estimate")


library(ggplot2)
library(gridExtra)




# CARBONIO (Z=6)
#_________________________


E.Z6 = c(seq(10^-10, 10, length=length(E.tot.6)))
out.Z6 = fE.cond_sz(E.Z6,z=c(1,0,0)) # pdf teorica scelta
ggplot(data.frame(x1 = E.Z6, x2 = E.tot.6, y1 = out.Z6)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "red",lwd=1)+ # funzione di densità teorica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 2.2),breaks = c(0,0.5,1,1.5,2))+
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
#__________________

E.Z14 = c(seq(10^-10, 10, length=length(E.tot.14)))
out.Z14 = fE.cond_sz(E.Z14,z=c(0,1,0))
ggplot(data.frame(x1 = E.Z14, x2 = E.tot.14, y1 = out.Z14)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "red",lwd=1)+ # funzione di densità teorica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 2.2),breaks = c(0,0.5,1,1.5,2))+
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


# FERRO (Z=26)
#_________________________


E.Z26 = c(seq(10^-10, 10, length=length(E.tot.26)))
out.Z26 = fE.cond_sz(E.Z26,z=c(0,0,1))
ggplot(data.frame(x1 = E.Z26, x2 = E.tot.26, y1 = out.Z26)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "red",lwd=1)+ # funzione di densità teorica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 2.2),breaks = c(0,0.5,1,1.5,2))+
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


# grafico classificazione errata carbonio

load("data_for_inference_only_photon.RData")
E.voxj = c(Simulation[Simulation$vox_real==22, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,1], Simulation[Simulation$vox_real==22, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,2],
           Simulation[Simulation$vox_real==22, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,3], Simulation[Simulation$vox_real==22, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,4])
E.voxj = E.voxj[!is.na(E.voxj)]


E.Z = c(seq(10^-10, 10, length=length(E.voxj)))
out.Z26 = fE.cond_sz(E.Z,z=c(0,0,1))
out.Z6 = fE.cond_sz(E.Z,z=c(1,0,0))
out.Z14 = fE.cond_sz(E.Z,z=c(0,1,0))
ggplot(data.frame(x1 = E.Z, x2 = E.voxj, y1 = out.Z26, y2 = out.Z6, y3 = out.Z14)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "darkred",lwd=1)+ # funzione di densità teorica ferro
  geom_line(aes(y = y2, x = x1), col= "darkblue",lwd=1)+ # funzione di densità teorica carbonio
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 2.2),breaks = c(0,0.5,1,1.5,2))+
  theme_minimal()+
  labs(title="Voxel 22"
  )+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))



#************************************************
# Quantità di fotoni provenienti da ciascun voxel
#************************************************


# conto quanti fotoni provengono da ciascun voxel

matrix.cont = matrix(NA,64,2)
matrix.cont[,1] = c(1:64)
colnames(matrix.cont) = c("voxel", "conteggi")

for(j in 1:64)
{
  # ricavo tutte le energie provenienti dal voxel reale j
  E.voxj = c(Simulation[Simulation$vox_real==j, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,1], Simulation[Simulation$vox_real==j, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,2],
             Simulation[Simulation$vox_real==j, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,3], Simulation[Simulation$vox_real==j, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,4])
  E.voxj = E.voxj[!is.na(E.voxj)]
  
  # conto quanti sono i fotoni provenienti dal voxel reale j
  cont = length(E.voxj)
  
  # salvo la quantità di fotoni provenienti dal voxel reale j
  matrix.cont[j,2] = cont
}

matrix.cont




library("lattice")

cube3

data.square1 <- data.frame(expand.grid(x = c("24.5","25.5","26.5","27.5"), y = c("3.5","2.5","1.5","0.5")), 
                           value = c(matrix.cont[13:16,2],
                                     matrix.cont[9:12,2],
                                     matrix.cont[5:8,2],
                                     matrix.cont[1:4,2]), cont = as.character(c(matrix.cont[13:16,2],
                                                                                matrix.cont[9:12,2],
                                                                                matrix.cont[5:8,2],
                                                                                matrix.cont[1:4,2])))


square1.graph =  
  xyplot(y ~ x, data = data.square1,
         panel = function(x,y, ...) {
           ltext(x = x, y = y, labels = data.square1$value, cex = 2, font = 4,
                 fontfamily = "HersheySans",col="darkred")
         }, main=list("z = 0.5", cex=2),
                                        scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                                        xlab=list("x", cex=2),
                                        ylab=list("y", cex=2),
                                        colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                                      height=1, width=2),
                                        col.regions = gray(round( seq(1,0, length=100),2)))

square1.graph

data.square2 <- data.frame(expand.grid(x = c("24.5","25.5","26.5","27.5"), y = c("3.5","2.5","1.5","0.5")), 
                           value = c(matrix.cont[29:32,2],
                                     matrix.cont[25:28,2],
                                     matrix.cont[21:24,2],
                                     matrix.cont[17:20,2]), cont = as.character(c(matrix.cont[29:32,2],
                                                                                  matrix.cont[25:28,2],
                                                                                  matrix.cont[21:24,2],
                                                                                  matrix.cont[17:20,2])))


square2.graph = 
  xyplot(y ~ x, data = data.square2,
         panel = function(y, x, ...) {
           ltext(x = x, y = y, labels = data.square2$cont, cex = 2, font = 4,
                 fontfamily = "HersheySans",col="darkred")
         }, main=list("z = 1.5", cex=2),
         scales=list(x=list(cex=1.5),y=list(cex=1.5)),
         xlab=list("x", cex=2),
         ylab=list("y", cex=2),
         colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                       height=1, width=2),
         col.regions = gray(round( seq(1,0, length=100),2)))

square2.graph



data.square3 <- data.frame(expand.grid(x = c("24.5","25.5","26.5","27.5"), y = c("3.5","2.5","1.5","0.5")), 
                           value = c(matrix.cont[45:48,2],
                                     matrix.cont[41:44,2],
                                     matrix.cont[37:40,2],
                                     matrix.cont[33:36,2]), cont = as.character(c(matrix.cont[45:48,2],
                                                                                  matrix.cont[41:44,2],
                                                                                  matrix.cont[37:40,2],
                                                                                  matrix.cont[33:36,2])))


square3.graph =  
  xyplot(y ~ x, data = data.square3,
         panel = function(y, x, ...) {
           ltext(x = x, y = y, labels = data.square3$cont, cex = 2, font = 4,
                 fontfamily = "HersheySans",col="darkred")
         },main=list("z = 2.5", cex=2),
         scales=list(x=list(cex=1.5),y=list(cex=1.5)),
         xlab=list("x", cex=2),
         ylab=list("y", cex=2),
         colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                       height=1, width=2),
         col.regions = gray(round( seq(1,0, length=100),2)))

square3.graph



data.square4 <- data.frame(expand.grid(x = c("24.5","25.5","26.5","27.5"), y = c("3.5","2.5","1.5","0.5")), 
                           value = c(matrix.cont[61:64,2],
                                     matrix.cont[57:60,2],
                                     matrix.cont[53:56,2],
                                     matrix.cont[49:52,2]), cont = as.character( c(matrix.cont[61:64,2],
                                                                                   matrix.cont[57:60,2],
                                                                                   matrix.cont[53:56,2],
                                                                                   matrix.cont[49:52,2])))


square4.graph = 
  xyplot(y ~ x, data = data.square4,
         panel = function(y, x, ...) {
           ltext(x = x, y = y, labels = data.square4$cont, cex = 2, font = 4,
                 fontfamily = "HersheySans",col="darkred")
         }, main=list("z = 3.5", cex=2),
         scales=list(x=list(cex=1.5),y=list(cex=1.5)),
         xlab=list("x", cex=2),
         ylab=list("y", cex=2),
         colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                       height=1, width=2),
         col.regions = gray(round( seq(1,0, length=100),2)))

square4.graph


