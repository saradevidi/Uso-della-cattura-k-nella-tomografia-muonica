# codice che genera la simulazione richiamando le funzioni create
#----------------------------------------------------------------

# carico le funzioni necessarie
source("create_geometry.R")
source("muon_energy_loss.R")
source("pdf_photon_energy.R")
source("absorption_length.R")
source("function_generatekCapture.R")


#*************************
# creo il volume incognito
#*************************

# creo i centri dei voxel del volume incognito e ne scelgo la relativa locazione nel piano 3D
cube1 = coord.volume(coord.ang = c(24,0,0), alt = 10, lato1 = 10, lato2 = 10)
cube1

# aggiungo le coordinate estreme in cui ciascun voxel è compreso
cube2 = coord.volume.Extreme(coord.volume.matrix = cube1)
cube2

# aggiungo la materia di cui sono composti i voxel (scelta casualmente)
cube3 = coord.volume.Z(coord.volume.Extreme.matrix = cube2, atoms = c(6,14,26), seed_ = 123)
cube3






#****************************
# creo il piano di rivelatore
#****************************

# coordinate per il detector
coord.riv1 = coord.riv(x.coord = 3.5, coord.volume = cube1, lato2 = 10, alt = 10) 
coord.riv1






#***************************************************************************************
# calcolo la perdita di energia del muone per ogni cm attraversato nei diversi materiali
#***************************************************************************************

# definisco una matrice con 2 colonne: 
#  - nella prima inserisco i numeri atomici Z 
#  - nella seconda inserisco la corrispondente densità in gr/cm^3
densitaZ = matrix(c(6, 14, 26, 2.25, 2.33, 7.87), ncol=2) 

# calcolo la perdita di energia per ogni Z attraversato
DeltaE.matrix = loss.energy(Zdensity = densitaZ, Y = 1)
DeltaE.matrix





#************************************************************
# definisco la probabilità di cattura k per i diversi materiali
#************************************************************

kcapturePROB.matrix = cbind(DeltaE.matrix$Z,c(0.08, 0.66, 0.91))
colnames(kcapturePROB.matrix) = c("Z", "prob_cattura")
kcapturePROB.matrix = as.data.frame(kcapturePROB.matrix)  





#****************************************************************************************************************
# definisco la funzione di sopravvivenza del fotone. Tale funzione prende in input il valore lambda a seconda del
# materiale in cui il fotone si trova, restituito da :
# funTotVec.Z14 (curva per il Silicio)
# funTotVec.Z6 (curva per il Carbonio)
# funTotVec.Z26 (curva per il Ferro);
# inoltre prende in input anche la densità del materiale e lo spessore attraversato.
#****************************************************************************************************************

Survival.lambda = function(spessore, densita, lambda) 
{ # probabilità che il fotone non sia stato assorbito nello spessore percorso
  exp( -(spessore * densita) / lambda)
}




#***************************************************************
# genero N muoni dalla funzione che simula il processo k-capture
#***************************************************************

Simulation = generate.kCapture(n.muon = 3000000, E_mu = c(2,90), coord.material = cube3, 
                  coord.detector = coord.riv1, DeltaE.matrix = DeltaE.matrix, 
                  kcapturePROB.matrix = kcapturePROB.matrix,
                  lato1 = 10, x = 0, seed_ = 8888888)


