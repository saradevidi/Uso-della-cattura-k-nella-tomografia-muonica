generate.kCapture = function(n.muon = 1000, # numero di muoni che voglio generare
                          E_mu = c(2,105), # limiti minimo e massimo dell'energia in MeV del fascio di muoni in entrata 
                          coord.material, # coordinate dei voxel nel volume con i corrispondenti numeri atomici Z 
                          coord.detector,  # coordinate dei voxel del rivelatore 
                          DeltaE.matrix, # dataframe con la perdita di energia del muone
                          kcapturePROB.matrix, # dataframe con le probabilità di cattura
                          lato1 = 10,# lato 1 del volume ("base")
                          x = 0, # coordinata x per il fascio di muoni in entrata
                          seed_ = 123) # set.seed per la riproducibilità
  
{
  x_partenza = x
  
  output = matrix(NA, nrow = n.muon, ncol = 37)
  colnames(output) = c("muon.exit", # indicatore che vale 1 se il muone è uscito dal volume poichè aveva ancora energia disponibile, 0 altrimenti
                       "muon.stop", # indicatore che vale 1 se il muone si è fermato nel volume poichè ha esaurito la sua energia in uno dei voxel al suo interno, 0 altrimenti
                       "kcapture", # indicatore che vale 1 se il muone è stato catturato dal nucleo dell'atomo, 0 altrimenti
                       "x_muon", "y_muon", "z_muon", # posizione del muone incidente rilevata dal detector (data la precisione del detector)
                       "vox_detector_muon",  # indicatore del voxel del detector in cui il muone incidente viene rivelato -> vox_detector_muon 
                                             # indica la riga della matrice coord.detector, cioè ad es. vox_detector_muon = 3 indica che la riga 3 della matrice
                                             # coord.detector è il voxel in cui il muone è stato rivelato. Se quel dato muone non viene visto allora è NA.
                       "E_muon", # energia in MeV del muone entrante
                       "E_fot1","E_fot2", "E_fot3", "E_fot4", # Energie dei fotoni rilevate dal detector (data la precisione del detector)
                       "x_fot1", "x_fot2","x_fot3", "x_fot4", # posizioni precise dei fotoni sul detector rilevate dal detector (data la precisione del detector)
                       "y_fot1", "y_fot2", "y_fot3", "y_fot4",
                       "z_fot1", "z_fot2", "z_fot3", "z_fot4",
                       "vox_detector_fot1", "vox_detector_fot2", "vox_detector_fot3",  "vox_detector_fot4", # indicatore del voxel del detector in cui il fotone viene rivelato -> vox_detector 
                                                                            # indica la riga della matrice coord.detector, cioè ad es. vox_detector_fot1 = 3 indica che la riga 3 della matrice
                                                                            # coord.detector è il voxel in cui il fotone è stato rivelato. Se quel dato fotone non viene rivelato allora è NA.
                       "theta_fot1", "phi_fot1", # angoli  dei fotoni reali (NON misurati dal detector)
                       "theta_fot2", "phi_fot2",
                       "theta_fot3", "phi_fot3",
                       "theta_fot4", "phi_fot4",
                       "vox_real"  # voxel reale in cui è avvenuto k-capture e/o il muone si è fermato (noto dalla simulazione) -> vox_real indica la riga della matrice coord.material,
                                   # cioè ad es. vox_real = 3 indica che la riga 3 della matrice coord.material è il voxel in cui il muone si è fermato.
  )
  
  set.seed(seed_)
  
  for(i in 1:n.muon)
  { # per ciascun muone emesso dal punto di coordinate x,y,z
    
    if(i%%1000 == 0){ # per tutti i muoni i-esimi multipli di 1000
    print(paste0("muone", sep=" ", i))} # ogni 1000 muoni stampo il punto a cui sono arrivata

    
    # definisco casualmente y e z
    y = runif(1, min = min(coord.material$y_min), max = max(coord.material$y_max)) # campiono casualmente una coordinata y per il muone in entrata, tra le coordinate min e max dell'asse y in cui il volume è definito 
    z = runif(1, min = min(coord.material$z_min), max = max(coord.material$z_max)) # campiono casualmente una coordinata z per il muone in entrata, tra le coordinate min e max dell'asse z in cui il volume è definito
    
    # il punto di partenza del muone è (x, y, z) .
    
    E_mu.new = runif(1, min = E_mu[1], max = E_mu[2]) # energia di partenza del muone (prima dell'attraversamento del tessuto)
    E_mu.forOutput = E_mu.new # mi serve per salvare l'energia del muone entrante da restituire nell'output
    
    # id.muon.exit e id.muon.stop identificano le 2 situazioni che possono verificarsi
    id.muon.exit = FALSE # variabile indicatrice che è TRUE se il muone generato esce dal volume perchè ha ancora energia disponibile, e dunque non avviene la cattura nel materiale
    id.muon.stop = FALSE # variabile indicatrice che è TRUE se il muone generato si ferma in un dato voxel, e dunque può avvenire la cattura-k
    id.kcapture = 0 # variabile indicatrice che è 1 se il muone è stato catturato dal nucleo "k-capture", 0 altrimenti
    voxel.true = NA # variabile che indica il voxel nella matrice coord.material corrispondente al voxel in cui il muone si è fermato -> se il
                    # muone non si ferma nel volume è NA
    
    direzioni.fotoni = matrix(NA, 4, 5) # matrice che conterrà sulle righe i fotoni emessi per ciascun muone, e nelle colonne gli angoli theta e phi corrispondenti (REALI), 
    # e le coordinate precise (REALI) della cattura k-capture da cui il fotone proviene (se questa avviene) -> questa matrice mi servirà per propagare il fotone fino al rivelatore e vedere se lo
    # interseca -> se lo interseca registro la posizione, e la fornirò in output insieme all'energia se quest'ultimo non viene assorbito. 
    colnames(direzioni.fotoni) = c("theta.true", "phi.true", "x.true", "y.true", "z.true") # x.true, y.true e z.true sono identici per ciascuno dei 4 fotoni poichè questi provengono dallo stesso muone.
    rownames(direzioni.fotoni) = c("fotone 1", "fotone 2", "fotone 3", "fotone 4")
    
    detector.coord.fotoni = matrix(NA, 4, 2) # matrice che conterrà sulle righe i fotoni emessi per ciascun muone, e nelle colonne le coordinate y e z del detector che corrispondono al punto
    # reale in cui il fotone interseca il detector -> questa matrice mi servirà per modificare tali coordinate secondo la precisione del detector e vedere per ciascun fotone se questi vengono poi 
    # effettivamente visti al detector.
    colnames(detector.coord.fotoni) = c("y.detect.real", "z.detect.real") 
    rownames(detector.coord.fotoni) = c("fotone 1", "fotone 2", "fotone 3", "fotone 4")
    
    energie.coord.fotoni = matrix(NA, 4, 5) # matrice che conterrà sulle righe i fotoni emessi, e nelle colonne l'energia rilevata, 
    # le posizioni dei fotoni sul detector, cioè RILEVATE dal detector e il voxel sul detector corrispondente a tale posizione.
    colnames(energie.coord.fotoni) = c("E.detect", "x.detect", "y.detect", "z.detect", "vox.detect")
    rownames(energie.coord.fotoni) = c("fotone 1", "fotone 2", "fotone 3", "fotone 4")
    
    direzioni.fotoni = as.data.frame(direzioni.fotoni)
    detector.coord.fotoni = as.data.frame(detector.coord.fotoni)
    
    x.det = coord.detector[1,1] # prendo la coordinata x del centro del rivelatore -> mi servirà per calcolare il punto di arrivo reale del fotone nelle coordinate (x,y,z) a x=x.det fissato
    
    for(j in 1:lato1) # in questo ciclo for individuo il voxel colpito dal muone a partire dalla prima coordinata x che identifica la posizione del volume, e poi proseguendo
    { # per tutte le coordinate ad ogni x successivo fino alla fine del volume (x.vox), ovvero modello la propagazione del muone nel materiale 
      # fino a quando non si ferma in un voxel di coordinate (x.vox, y, z).
      # for(j in 1:lato1) => per tutta la larghezza del volume (lato1)
      
      
      # 1. Identifico il voxel in cui entra il muone
      
      x.vox = coord.material[j,1] # centro del voxel in x
      # (x.vox, y, z) definiscono le coordinate del voxel colpito dal muone in x = x.vox poichè il muone prosegue con una traiettoria dritta
      
      # a quale voxel corrispondono le coordinate (x.vox, y, z) ? Che numero atomico Z ha ?
      
      # individuo il numero atomico del voxel colpito
      nr.atomico.vox = coord.material$Z[((coord.material$x==x.vox) & (coord.material$y_min <= y) & (coord.material$y_max >= y) & (coord.material$z_min <= z) & (coord.material$z_max >= z))]
      
      
      # 2. Definisco la perdita di energia DeltaE del muone nell'attraversare tutto il voxel (x.vox, y, z).
      
      # DeltaE.matrix$perditaEnergia.1cm[DeltaE.matrix$Z==nr.atomico.vox]  # Energia che il muone perde una volta che ha completato l'attraversamento del voxel in cui si trova (percorre 1 cm)
      E_mu.new = E_mu.new - DeltaE.matrix$loss.energy.Ycm[DeltaE.matrix$Z==nr.atomico.vox]  # energia del muone dopo che ha attraversato completamente il voxel in cui si trova
      
      
      # 3. Se il muone ha esaurito la sua energia allora il muone si ferma nel voxel identificato dalle coordinate (x.vox, y, z) => 
      #    può avvenire la cattura-k in quel voxel => simulo la cattura del muone => uscirò da questo ciclo for poichè il muone esaurisce la 
      #    sua energia quando E_mu.new <= 0 .
      if(E_mu.new <= 0)
      {
        id.muon.stop = TRUE # le coordinate in cui il muone si è fermato saranno approssimativamente (x.vox, y, z).
        # Ricavo una coordinata x più precisa in cui il muone effettivamente si ferma. Ciò lo faccio con la precisione di 1 mm, ovvero suddivido la perdita di energia DeltaE per 1 cm 
        # percorso in x nella perdita di energia per 1 mm percorso in x dividendo quest'ultima per 10. Poi, a partire dall'energia che il muone possiede a inizio voxel, sottraggo la perdita
        # di energia mm per mm, e identifico in questo modo la coordinata x avanzando di mm in mm dentro al voxel in cui sò che il muone si è fermato:
        E_iniziale = E_mu.new + DeltaE.matrix$loss.energy.Ycm[DeltaE.matrix$Z==nr.atomico.vox] # energia del muone all'inizio del voxel identificato in cui sò che il muone si ferma
        DeltaE.1mm = DeltaE.matrix$loss.energy.Ycm[DeltaE.matrix$Z==nr.atomico.vox]/10 # Energia che il muone perde attraversando 1 mm in x del voxel in cui si trova
        x = coord.material$x_min[((coord.material$x==x.vox) & (coord.material$y_min <= y) & (coord.material$y_max >= y) & (coord.material$z_min <= z) & (coord.material$z_max >= z))] # coordinata x di inizio voxel
        
        for(l in 1:10)
        {# da 1 mm a 10 mm attraversati
          E_iniziale = E_iniziale - DeltaE.1mm # energia effettiva del muone dopo l'attraversamento di 1 mm del voxel in questione
          x = x + 0.1 # coordinata x nella quale il muone avanza -> si sposta di 1 mm ogni volta
          if(E_iniziale <= 0)
          { # se il muone ha esaurito la sua energia allora si è fermato in quel mm, identificato dalla coordinata x calcolata
            break # esco da questo ciclo for => x a questo punto è la coordinata x in cui il muone si è fermato
          }
        }
        
        break
      }
    }
    
    if((E_mu.new > 0) & (j == lato1))
    { # se il muone ha ancora energia e sono arrivata alla fine del ciclo for precedente, cioè il muone ha attraversato tutto il volume,
      # allora GENERO UN ALTRO MUONE
      id.muon.exit = TRUE
    }
    
    
    # 4. Se il muone si è fermato in un voxel del volume allora simulo la cattura nel materiale, e quindi l'eventuale k-capture nelle coordinate (x, y, z)
    
    if((id.muon.stop == TRUE) & (id.muon.exit==FALSE))
    { # se il muone si è fermato in un voxel del volume
      
      # salvo il valore del voxel in cui il muone si è fermato
      voxel.true =  which(((coord.material$x==x.vox) & (coord.material$y_min <= y) & (coord.material$y_max >= y) & (coord.material$z_min <= z) & (coord.material$z_max >= z))) # riga della matrice
                                                                                                                              # coord.material corrispondente al voxel in cui il muone si è fermato
      
      
      # il muone in quel voxel decade come particella libera o viene catturato dal nucleo ?
      id.kcapture = rbinom(1,size=1,prob=kcapturePROB.matrix$prob_cattura[kcapturePROB.matrix$Z==nr.atomico.vox]) # genero da una bernoulliana con probabilità data
                                                                                                                  # dalla probabilità di cattura nel materiale Z del voxel
      
      
      if(id.kcapture == 1)
      { # se il muone viene catturato dal nucleo (id.kcapture == 1) allora simulo la cattura k-capture nel voxel individuato dalle coordinate (x, y, z),
        # se invece il muone NON è stato catturato dal nucleo allora GENERO UN ALTRO MUONE e questo if non viene svolto 
        
        # decido con un numero casuale quanti fotoni vengono emessi (decido che possono essere al massimo 4)
        n.fotoni = sample(c(1,2,3,4), size=1, prob = c(0.1, 0.1, 0.3, 0.5))
        
        # determino i valori phi_min, phi_max, th_min e th_max per quei fotoni che possono essere emessi dal voxel in cui il muone si è fermato nelle coordinate (x,y,z)
        # -> (x,y,z) nel codice corrisponde a (x0, y0,z0) nelle formule (3.2), (3.3), (3.4):
        
        #---> phi_min
        ro.phi_min = sqrt( (coord.detector[1,1]-x)^2 + (max(coord.detector$z_max) - z)^2 )
        phi_min = acos((max(coord.detector$z_max) - z)/ro.phi_min)
        
        #---> phi_max
        ro.phi_max = sqrt( (coord.detector[1,1]-x)^2 + (min(coord.detector$z_min) - z)^2 )
        phi_max = acos((min(coord.detector$z_min) - z)/ro.phi_max)
        
        #---> th_min
        th_min = atan( (max(coord.detector$y_max)-y)/(coord.detector[1,1]-x) ) + pi
        
        #---> th_max
        th_max = atan( (min(coord.detector$y_min)-y)/(coord.detector[1,1]-x) ) + pi
        
        
        
        name = pdf.fun$rfun.name[pdf.fun$Z==nr.atomico.vox] # nome della funzione che mi serve per campionare l'energia del fotone, cioè è la funzione che
                                                            # corrisponde al numero atomico del voxel in cui il fotone è stato emesso
        r.pdf = get(name) # converte la stringa di testo name in un comando
        
        
        
        for(fot in 1:n.fotoni)
        { # per tutti i fotoni emessi
          
          # campiono una energia per il fotone dalla pdf che fa riferimento al numero atomico Z del voxel in cui il fotone viene emesso
          E_foton = r.pdf(1)
          
         
          if(E_foton >= 0.001)
          { # se l'energia è >=1KeV allora mi interessa generare 2 angoli per il fotone, e dunque propagarlo fino al rivelatore.
            # Se l'energia campionata è < 1KeV allora so già che questo fotone viene subito assorbito, pertanto lo escludo da subito e passo al fotone successivo.
          
            t.indicator = 0
            while(t.indicator == 0)
            {# per ogni fotone emesso genero casualmente 2 angoli theta e phi che descrivono la sua traiettoria nello spazio 3-dimensionale (genero solo angoli 
            # che permettono al fotone di colpire effettivamente il rivelatore -> potrebbe essere che qualche volta comunque gli angoli non permettano di colpirlo, per
            # questo faccio un controllo e se gli angoli non sono ammissibili li getto via e ne genero altri 2 fintanto che non ho generato angoli ammissibili per
            # ciascuno dei miei n.fotoni)
               th.fotone = runif(1, min = th_min, max = th_max) # l'angolo theta viene campionato con uguale probabilità tra [th_min, th_max]
               phi.fotone = runif(1, min = phi_min, max = phi_max) # l'angolo phi viene campionato con uguale probabilità tra [phi_min, phi_max]
            
               ro = (x.det - x)/(sin(phi.fotone)*cos(th.fotone))
               y.det = y + (ro * sin(phi.fotone)*sin(th.fotone)) # vera coordinata y di arrivo del fotone nel rivelatore
               z.det = z + (ro * cos(phi.fotone))  # vera coordinata z di arrivo del fotone nel rivelatore
            
               # le coordinate (x.det, y.det, z.det) appartengono al rivelatore ? -> il fotone ha intersecato il rivelatore ?
               if( sum( (coord.detector$x==x.det) & (coord.detector$y_min <= y.det) & (coord.detector$y_max >= y.det) & (coord.detector$z_min <= z.det) & (coord.detector$z_max >= z.det) ) == 1 )
               {
                 t.indicator = 1 # salvo il fotone e non proseguo il ciclo while
                 direzioni.fotoni[fot,] = c(th.fotone, phi.fotone, x, y, z)
                 detector.coord.fotoni[fot,] = c(y.det, z.det)
               } # se gli angoli generati non gli permettono di intersecare il rivelatore, allora questo if non viene svolto, e genero altri 2 angoli (cioè continuo il ciclo while)
             }
          
          # A questo punto ho degli angoli theta e phi ammissibili per il fotone emesso (th.fotone e phi.fotone), cioè il fotone emesso ha angoli che incrociano il rivelatore; inoltre ho per 
          # ciascun fotone le coordinate reali in cui questi incrociano il detector (y.det e z.det).
          
          # Per ciascun fotone controllo se il fotone emesso potrebbe essere visto dal piano di rivelatore. Poi, se il fotone può essere visto
          # dal rivelatore, allora  propago il fotone fino al rivelatore calcolando di volta in volta la probabilità che questo venga assorbito (con la funzione di sopravvivenza), e, se 
          # il fotone riesce ad uscire dal volume senza essere assorbito, allora registro la posizione che viene rilevata nel piano di rivelatore -> questo lo faccio per ciascun fotone.
          # il punto di partenza del fotone è (x, y, z) = (direzioni.fotoni$x.true, direzioni.fotoni$y.true, direzioni.fotoni$z.true) con angoli 
          # th.fotone e phi.fotone = (direzioni.fotoni$theta.true, direzioni.fotoni$phi.true)
          
            
                 survival = NA 
          
          # La posizione misurata dal rivelatore non è esattamente la posizione vera, e quindi non è esattamente (x.det, y.det, z.det). Perciò
          # tali coordinate ricavate vanno leggermente modificate, e dipenderanno dalla risoluzione sulla posizione del rivelatore.
          # y.det.misurata e z.det.misurata sono quindi le coordinate effettivamente misurate nel rivelatore, pertanto devo anche verificare poi nuovamente che queste coordinate
          # y.det.misurata e z.det.misurata siano delle coordinate ammissibili, ovvero che rientrino nella geometria del rivelatore.
                 y.det.misurata = detector.coord.fotoni$y.detect.real[fot] + rnorm(1, mean = 0, sd = 0.1)
                 z.det.misurata = detector.coord.fotoni$z.detect.real[fot] + rnorm(1, mean = 0, sd = 0.1)
          
                 if( sum( (coord.detector$x==x.det) & (coord.detector$y_min <= y.det.misurata) & (coord.detector$y_max >= y.det.misurata) & (coord.detector$z_min <= z.det.misurata) & (coord.detector$z_max >= z.det.misurata) ) == 1 )
                 { # se le coordinate y.det.misurata e z.det.misurata sono delle coordinate ammissibili (=> il fotone può essere visto dal rivelatore), allora:
            
                   # salvo il valore del voxel in cui il fotone ha incrociato il rivelatore (voxel misurato, non reale)
                    voxel.detector =  which(((coord.detector$x==x.det) & (coord.detector$y_min <= y.det.misurata) & (coord.detector$y_max >= y.det.misurata) & (coord.detector$z_min <= z.det.misurata) & (coord.detector$z_max >= z.det.misurata))) # riga della matrice
                                                                                                                                                                                               # coord.detector corrispondente al voxel in cui il fotone è stato rivelato
                    survival = TRUE # inizialmente il fotone è vivo
            
            
            #    Propago il fotone nei voxel del volume e definisco di volta in volta la probabilità di sopravvivenza di quest'ultimo fino al rivelatore
            #    -> se alla fine del volume il fotone è sopravvissuto allora dichiaro che l'ho rivelato e salvo la sua posizione misurata nel detector (x.det, y.det.misurata, z.det.misurata)
            #    con la sua energia corrispondente.
            
                    tot = which(coord.material$x==x.vox)[1]   # mi dice di quanto devo retrocedere verso sinistra nell'asse x in cui il volume è definito per propagare il fotone 
            # fino al rivelatore 
            # se tot = 1 significa che il fotone si trova nella prima fila di voxel da sinistra del volume, cioè nel primo strato da sinistra  => mi è sufficiente calcolare la distanza percorsa dal
            # fotone fino alla coordinata x iniziale fissata in cui il cubo è definito (che corrisponde a x_min del voxel in cui il fotone si trova), e non ho ovviamente bisogno di
            # propagarlo ulteriormente nel materiale.
            
            # estraggo le coordinate da cui il fotone parte, e successivamente queste verranno aggiornate di volta in volta nel ciclo successivo per identificare le coordinate in cui
            # il fotone si trova man mano che si sposta lungo l'asse x fino ad arrivare a fine volume (se ci arriva, e ciò dipende dalla funzione di sopravvivenza)
                    x.fotone = direzioni.fotoni$x.true[fot]
                    y.fotone = direzioni.fotoni$y.true[fot]
                    z.fotone = direzioni.fotoni$z.true[fot]
            
            # calcolo il punto x_min del voxel da cui il fotone parte. Tale valore verrà aggiornato di volta in volta man mano che il fotone si sposta da destra verso sinistra.
            # L'aggiornamento è identificato dal valore di x_min_j2.
                    x_min_voxel.true = coord.material$x_min[voxel.true]
            
            
                    for(j2 in 1:tot)
                    { # per ogni punto di coordinata x_min dal voxel in cui il fotone si trova fino a fine volume (cioè da destra verso sinistra)
              
                      if(j2 > 1)
                      {# tutte le volte in cui il fotone si sposta di uno strato in x successivo al primo da cui parte, aggiorno x_min_j2 e le coordinate del fotone man mano che si sposta in x
                
                       # calcolo il punto di coordinate x,y,z che il fotone interseca quando entra nel voxel successivo al precedente. Con queste coordinate riesco ad ottenere poi il numero atomico Z del
                       # voxel che il fotone interseca quando entra nel voxel successivo al precedente (ciò perchè considero i numeri atomici di solo quei voxel che il fotone interseca quando entra nel voxel successivo)
                         x.fotone = x_min_j2 - 0.1 # entro di  1 mm nel voxel successivo e, fissata questa x ottengo le coordinate y e z corrispondenti
                         ro = (x.fotone - direzioni.fotoni$x.true[fot])/(sin(direzioni.fotoni$phi.true[fot])*cos(direzioni.fotoni$theta.true[fot]))
                         y.fotone = direzioni.fotoni$y.true[fot] + (ro * sin(direzioni.fotoni$phi.true[fot]) * sin(direzioni.fotoni$theta.true[fot]))
                         z.fotone = direzioni.fotoni$z.true[fot] + (ro * cos(direzioni.fotoni$phi.true[fot]))
                
                       # aggiorno x_min_j2
                         x_min_j2 = x_min_voxel.true - (j2 -1)
                
                      }else{
                         x_min_j2 = x_min_voxel.true
                      }
              
              # I. calcolo la distanza percorsa dal fotone dal punto in cui si trova (x.fotone, y.fotone, z.fotone) fino alla coordinata x_min del voxel in cui si trova di volta in volta.
              
                      if(j2 == 1)
                      {
                        Deltax = x.fotone - x_min_j2  # distanza percorsa in x dal fotone per arrivare a inizio voxel successivo (sarà sempre pari a 1cm se il fotone NON si trova nel primo voxel da cui parte)
                        d = sqrt(1 + (tan(direzioni.fotoni$phi.true[fot]*cos(direzioni.fotoni$theta.true[fot])))^2 + (tan(direzioni.fotoni$phi.true[fot]*sin(direzioni.fotoni$theta.true[fot])))^2) * Deltax # distanza percorsa dal
                                                                                                                                                                  # fotone se questo si sposta di una distanza Deltax in x.
                      }
                      if(j2 == 2)
                      { # se j2 > 1 allora Deltax sarà sempre pari a 1, pertanto d sarà sempre uguale per tutti i j2 >= 2. Dunque calcolo d solo per j2 = 2, e per j2 >=2 resterà in memoria il d qui calcolato
                         d = sqrt(1 + (tan(direzioni.fotoni$phi.true[fot]*cos(direzioni.fotoni$theta.true[fot])))^2 + (tan(direzioni.fotoni$phi.true[fot]*sin(direzioni.fotoni$theta.true[fot])))^2) # poichè Deltax = 1 
                      } 
              
              
              # II. identifico il voxel in cui il fotone si trova e dunque ne ricavo il numero atomico Z.
              
                      Z.j2 = coord.material$Z[((coord.material$x_min == x_min_j2) & (coord.material$y_min <= y.fotone) & (coord.material$y_max >= y.fotone) & (coord.material$z_min <= z.fotone) & (coord.material$z_max >= z.fotone))]
              
              
              # III. calcolo la probabilità di sopravvivenza utilizzando la funzione di sopravvivenza, e decido se il fotone sopravvive o meno.
              
                     # calcolo il lambda
                      name.fun.lambda = function.lambda$name[function.lambda$Z == Z.j2] # nome della funzione che mi serve per estrarre il lambda data l'energia del fotone e il materiale in cui il fotone si trova
                      fun.lambda = get(name.fun.lambda) # converte la stringa di testo name in un comando
                      lambda.j2 = fun.lambda(E_foton) # lambda corrispondente all'energia del fotone nel materiale in cui attualmente si trova
              
                      # estraggo la densità del materiale Z.j2
                      densita.j2 = DeltaE.matrix$density[DeltaE.matrix$Z == Z.j2]
              
                      # calcolo la probabilità di sopravvivenza
                      prob.S = Survival.lambda(spessore = d, densita = densita.j2, lambda = lambda.j2)
              
                      # decido se il fotone viene assorbito o se passa allo strato di voxel successivo
                      if(rbinom(1, 1, prob.S) == 0)
                      {# se il fotone è stato assorbito esco da questo ciclo e passo al fotone successivo
                          survival = FALSE
                          break 
                      }
                    }
            
                    if(j2 == tot & survival == TRUE)
                    {# se il fotone è sopravvissuto fino alla fine del volume 
                      
                      # modello la probabilità che il fotone venga rivelato dal tracker
                      # a seconda della sua energia:  prob=1 se E>=1 MeV, P=E se E<1 MeV
                      if(E_foton>=1)
                      {
                        probabilita = 1
                      }else{
                        probabilita = E_foton
                      }
                      
                      if(rbinom(1,1,prob = probabilita)==1)
                      {# se il fotone viene rivelato dal detector
                      # salvo le coordinate e il voxel corrispondente con cui il fotone è stato rivelato dal detector con la sua corrispondente energia rilevata
                      energie.coord.fotoni[fot,1] = abs(E_foton + rnorm(1, mean=0, sd=sqrt(0.001))) 
                      energie.coord.fotoni[fot,2:4] = c(x.det,y.det.misurata,z.det.misurata)
                      energie.coord.fotoni[fot,5] = voxel.detector
                      }
                    }
            
                 }
          }

        } # procedo col fotone successivo
      }
      
    }
    
    
    
    # assumiamo la posizione del muone che incide nel materiale con la stessa risoluzione in posizione (y,z) dei fotoni
    y.muon.misurata = y + rnorm(1, mean = 0, sd = 0.1) 
    z.muon.misurata = z + rnorm(1, mean = 0, sd = 0.1) 

    # verifico se (y,z) sono compatibili con la geometria del detector
    if( sum( (y.muon.misurata >= coord.detector$y_min) & (y.muon.misurata <= coord.detector$y_max) & (z.muon.misurata >= coord.detector$z_min) & (z.muon.misurata<=coord.detector$z_max) ) == 1 )
    { # se la coordinata del muone entrante misurata rientra nella geometria del rivelatore
      voxel.muon.misurata = which(((coord.detector$x==x.det) & (coord.detector$y_min <= y.muon.misurata) & (coord.detector$y_max >= y.muon.misurata) & (coord.detector$z_min <= z.muon.misurata) & (coord.detector$z_max >= z.muon.misurata))) # riga della matrice
                                                                                                                                                                                         # coord.detector corrispondente al voxel in cui il muone è stato rivelato
    }else{
      y.muon.misurata = NA
      z.muon.misurata = NA
      voxel.muon.misurata = NA                                                                                                                                                                          
    }

    # per ogni muone generato metto tutti i suoi output (a prescindere dal fatto che sia avvenuta k-capture o meno)
      output[i,] = c(id.muon.exit, id.muon.stop, id.kcapture, 
                   x_partenza, y.muon.misurata, z.muon.misurata, # posizione del muone incidente (punto da cui il muone parte) -> rilevata dal detector. Avrò degli NA in (y,z) se il detector non vede il muone entrante
                   voxel.muon.misurata, # voxel del detector colpito dal muone incidente -> rilevata dal detector. Avrò degli NA se il detector non vede il muone entrante
                   E_mu.forOutput, # energia del muone entrante
                   as.vector(energie.coord.fotoni[,1]), # Energie dei fotoni rilevate dal detector (se non rilevo fotoni avrò NA)
                   as.vector(energie.coord.fotoni[,2]), # posizioni dei fotoni sul detector rilevate dal detector (se non rilevo fotoni avrò NA)
                   as.vector(energie.coord.fotoni[,3]),
                   as.vector(energie.coord.fotoni[,4]),
                   as.vector(energie.coord.fotoni[,5]),
                   as.numeric(direzioni.fotoni[1,1:2]), # angoli  dei fotoni reali, NON rilevati dal detector (se ho emesso meno di 4 fotoni avrò degli NA), gli angoli vengono riportati 
                   as.numeric(direzioni.fotoni[2,1:2]), # in output anche se il fotone non viene rivelato ma k-capture è avvenuto. Gli angoli NON vengono riportati (avrò NA) se l'energia del fotone
                   as.numeric(direzioni.fotoni[3,1:2]), # emesso è < 1KeV poichè questo viene assorbito subito.
                   as.numeric(direzioni.fotoni[4,1:2]),
                   voxel.true # riga della matrice coord.material che individua il VERO voxel in cui k-capture è avvenuto e/o il muone si è fermato
                   )
    }
    as.data.frame(output)
}









