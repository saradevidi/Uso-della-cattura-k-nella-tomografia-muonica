# unità di misura su cui si basano tutte le funzioni :
# 1 voxel = 1cm x 1cm x 1cm





#+++++++++++++++++++++++
# DEFINIZIONE DEL VOLUME
#+++++++++++++++++++++++


coord.volume = function(coord.ang = c(24,0,0), alt = 10, lato1 = 10, lato2= 10) # cubo 10cm X 10cm X 10cm di default
{ # funzione che dati :
  #    - coord.ang: la coordinata del primo angolo alla base del volume da sx (davanti)
  #    - alt: altezza del volume in voxel
  #    - lato1: lato 1 del volume visto dall'alto ("base")
  #    - lato2: lato 2 del volume visto dall'alto ("altezza")
  # restituisce una matrice di coordinate per il centro di ciascun voxel
  # (supponendo che il volume venga posizionato in modo parallelo all'asse x)
  nvox = lato1*lato2 # numero di voxel in uno strato del volume -> gli strati sono definiti dall'altezza del volume
  coord.volume.out = matrix(NA, nvox, 3) # matrice che conterrà sulle righe l'id dei diversi voxel e sulle colonne le coordinate x,y,z corrispondenti
  # inizialmente coord.volume.out contiene solo i voxel del primo strato
  
  # definisco tutte le coordinate dei voxel nel primo strato
  ncambi = lato1-1
  narr = 1 # voxel 1 
  # narr tiene conto del numero del voxel di cui viene di volta in volta determinata la coordinata centrale
  for(i in 1:lato2)
  { # per tutte "le righe" del volume nel primo strato
    if(i == 1)
    { 
      coord.volume.out[narr,] = coord.ang + 0.5
      narr = narr + 1
      if(ncambi>=1)
      {
        for(j in 1:ncambi)
        {
        coord.volume.out[narr,] =  c(coord.volume.out[(narr-1),1] + 1, coord.volume.out[(narr-1),2], coord.volume.out[(narr-1),3])
        narr = narr + 1
        }
      }
    }else if(i != 1){
      coord.volume.out[narr,] = c(coord.volume.out[1,1], coord.volume.out[1,2]+(i-1), coord.volume.out[1,3])
      narr = narr + 1
      if(ncambi>=1)
      {
        for(j in 1:ncambi)
        {
          coord.volume.out[narr,] =  c(coord.volume.out[(narr-1),1] + 1, coord.volume.out[(narr-1),2], coord.volume.out[(narr-1),3])
          narr = narr + 1
        }
      }
    }
  }
  
  
  # definisco tutte le coordinate dei voxel negli strati successivi
  for(str in 1:(alt-1))
  {
    # coordinate dei voxel nello strato str
    strato.mat = cbind(coord.volume.out[1:(narr-1),1:2],coord.volume.out[1:(narr-1),3]+str) # in tutti gli strati successivi al primo ciò che cambia è solo la coordinata z
    coord.volume.out = rbind(coord.volume.out,strato.mat)
  }
  
  colnames(coord.volume.out) = c("x", "y", "z")
  coord.volume.out
}






coord.volume.Extreme = function(coord.volume.matrix)
{ # funzione che data :
  #    - coord.volume.matrix: matrice ottenuta con il comando coord.volume
  # restituisce la matrice di partenza con in aggiunta le coordinate estreme in cui ciascun voxel è compreso
  x_min = coord.volume.matrix[,colnames(coord.volume.matrix)=="x"] - 0.5
  x_max = coord.volume.matrix[,colnames(coord.volume.matrix)=="x"] + 0.5
  
  y_min = coord.volume.matrix[,colnames(coord.volume.matrix)=="y"] - 0.5
  y_max = coord.volume.matrix[,colnames(coord.volume.matrix)=="y"] + 0.5
  
  z_min = coord.volume.matrix[,colnames(coord.volume.matrix)=="z"] - 0.5
  z_max = coord.volume.matrix[,colnames(coord.volume.matrix)=="z"] + 0.5
  
  coord.volume.Extreme.out = cbind(coord.volume.matrix, x_min, x_max, y_min, y_max, z_min, z_max)
  coord.volume.Extreme.out
}




coord.volume.Z = function(coord.volume.Extreme.matrix, atoms = c(6, 14, 26), seed_ = 123) # di default tre materiali : Z = 6 (Carbonio), Z = 14 (Silicio) e Z = 26 (Ferro)
{ # funzione che data :
  #    - coord.volume.Extreme.matrix: matrice ottenuta con il comando coord.volume.Extreme
  #    - atoms: vettore numerico che riporta i numeri atomici Z di cui si vuole che il volume sia composto
  # restituisce la matrice di partenza con in aggiunta la materia di cui sono composti i singoli voxel (assegnata casualmente tra gli Z possibili indicati a priori).
  
  set.seed(seed_)
  Z = sample(atoms, dim(coord.volume.Extreme.matrix)[1], replace = T)
  
  coord.volume.Z.out = cbind(coord.volume.Extreme.matrix, Z)
  coord.volume.Z.out = as.data.frame(coord.volume.Z.out)
  coord.volume.Z.out
}











#++++++++++++++++++++++++++++++++++++
# DEFINIZIONE DEL PIANO DI RIVELATORE
#++++++++++++++++++++++++++++++++++++


coord.riv = function(x.coord = 3.5,  # coordinata x (centro del voxel) del rivelatore
                     coord.volume , # matrice di coordinate del volume ricavata con la funzione coord.volume
                     lato2 = 10, # larghezza del volume in voxel (=altezza del volume se lo immagino dall'alto)
                     alt = 10) # altezza del volume in voxel 
{ # ricava la posizione dei voxel dei rivelatori usando la posizione del volume creato 
  # (Si suppone che i rivelatori siano posti alle stesse coordinate y e z del volume e abbiano larghezza di un voxel).
  # Restituisce un data.frame di coordinate per il rivelatore in posizione x.coord .
  
  coord.riv.out = matrix(NA,(lato2*alt),3) # matrice che conterrà le coordinate dei voxel nel rivelatore
  
  x.ref = coord.volume[1,1] # leggo la prima coordinata x della matrice in coord.cube
  
  mat.ref = coord.volume[coord.volume[,1]==x.ref,][,2:3] # leggo tutte le coordinate y,z del cubo relative ad un x fissato (x.ref)
  
  coord.riv.out = cbind(rep(x.coord,dim(mat.ref)[1]),mat.ref) # modifico soltanto la coordinata x
  colnames(coord.riv.out) = c("x", "y", "z")
  
  
  # Definisco le coordinate in cui ciascun voxel è compreso 
  x_min = coord.riv.out[,colnames(coord.riv.out)=="x"] - 0.5
  x_max = coord.riv.out[,colnames(coord.riv.out)=="x"] + 0.5
  
  y_min = coord.riv.out[,colnames(coord.riv.out)=="y"] - 0.5
  y_max = coord.riv.out[,colnames(coord.riv.out)=="y"] + 0.5
  
  z_min = coord.riv.out[,colnames(coord.riv.out)=="z"] - 0.5
  z_max = coord.riv.out[,colnames(coord.riv.out)=="z"] + 0.5
  
  coord.riv.out = cbind(coord.riv.out, x_min, x_max, y_min, y_max, z_min, z_max) # ora ogni voxel ha le coordinate in cui ciascun voxel è definito
  coord.riv.out = as.data.frame(coord.riv.out)
  
  coord.riv.out
}








