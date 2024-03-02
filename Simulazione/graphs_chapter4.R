

library("lattice")




#***************************************************
#*1. Distribuzione della posizione di stop dei muoni
#***************************************************

#------------------------------------------------------------------------------------------------
# creo i grafici che mi rappresentano la quantità di muoni che si sono fermati in ciascun voxel.
# Mostro il volume dall'alto "a fette", per ogni strato in z.
#------------------------------------------------------------------------------------------------

coord.x = c("24.5","25.5","26.5","27.5","28.5","29.5","30.5","31.5","32.5","33.5")
coord.y = c("0.5","1.5","2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5")


# grafico quantità di muoni nei primi 100 voxel (primo strato in z=0.5)
square1 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[1:10],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[11:20],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[21:30],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[31:40],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[41:50],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[51:60],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[61:70],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[71:80],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[81:90],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[91:100])
colnames(square1) <- coord.x                    
rownames(square1) <- coord.y  
square1.graph = levelplot(square1, main=list("z = 0.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square1.graph

# grafico quantità di muoni nei successivi 100 voxel (secondo strato in z=1.5)
square2 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[101:110],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[111:120],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[121:130],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[131:140],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[141:150],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[151:160],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[161:170],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[171:180],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[181:190],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[191:200])
colnames(square2) <- coord.x                    
rownames(square2) <- coord.y  
square2.graph = levelplot(square2, main=list("z = 1.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square2.graph

# grafico quantità di muoni nei successivi 100 voxel (terzo strato in z=2.5)
square3 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[201:210],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[211:220],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[221:230],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[231:240],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[241:250],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[251:260],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[261:270],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[271:280],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[281:290],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[291:300])
colnames(square3) <- coord.x                    
rownames(square3) <- coord.y  
square3.graph = levelplot(square3, main=list("z = 2.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square3.graph

# grafico quantità di muoni nei successivi 100 voxel (quarto strato in z=3.5)
square4 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[301:310],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[311:320],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[321:330],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[331:340],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[341:350],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[351:360],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[361:370],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[371:380],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[381:390],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[391:400])
colnames(square4) <- coord.x                    
rownames(square4) <- coord.y  
square4.graph = levelplot(square4, main=list("z = 3.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square4.graph

# grafico quantità di muoni nei successivi 100 voxel (quinto strato in z=4.5)
square5 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[401:410],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[411:420],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[421:430],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[431:440],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[441:450],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[451:460],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[461:470],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[471:480],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[481:490],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[491:500])
colnames(square5) <- coord.x                    
rownames(square5) <- coord.y  
square5.graph = levelplot(square4, main=list("z = 4.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square5.graph

# grafico quantità di muoni nei successivi 100 voxel (sesto strato in z=5.5)
square6 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[501:510],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[511:520],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[521:530],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[531:540],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[541:550],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[551:560],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[561:570],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[571:580],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[581:590],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[591:600])
colnames(square6) <- coord.x                    
rownames(square6) <- coord.y  
square6.graph = levelplot(square6, main=list("z = 5.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square6.graph

# grafico quantità di muoni nei successivi 100 voxel (settimo strato in z=6.5)
square7 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[601:610],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[611:620],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[621:630],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[631:640],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[641:650],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[651:660],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[661:670],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[671:680],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[681:690],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[691:700])
colnames(square7) <- coord.x                    
rownames(square7) <- coord.y  
square7.graph = levelplot(square7, main=list("z = 6.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square7.graph

# grafico quantità di muoni nei successivi 100 voxel (ottavo strato in z=7.5)
square8 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[701:710],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[711:720],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[721:730],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[731:740],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[741:750],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[751:760],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[761:770],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[771:780],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[781:790],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[791:800])
colnames(square8) <- coord.x                    
rownames(square8) <- coord.y  
square8.graph = levelplot(square8, main=list("z = 7.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square8.graph

# grafico quantità di muoni nei successivi 100 voxel (nono strato in z=8.5)
square9 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[801:810],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[811:820],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[821:830],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[831:840],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[841:850],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[851:860],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[861:870],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[871:880],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[881:890],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[891:900])
colnames(square9) <- coord.x                    
rownames(square9) <- coord.y  
square9.graph = levelplot(square9, main=list("z = 8.5", cex=2), sub="stop muoni",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square9.graph

# grafico quantità di muoni nei successivi 100 voxel (decimo strato in z=9.5)
square10 = rbind(table(Simulation$vox_real[!is.na(Simulation$vox_real)])[901:910],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[911:920],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[921:930],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[931:940],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[941:950],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[951:960],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[961:970],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[971:980],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[981:990],
                table(Simulation$vox_real[!is.na(Simulation$vox_real)])[991:1000])
colnames(square10) <- coord.x                    
rownames(square10) <- coord.y  
square10.graph = levelplot(square10, main=list("z = 9.5", cex=2), sub="stop muoni",
                           scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                           xlab=list("y", cex=2),
                           ylab=list("x", cex=2),
                           colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                         height=1, width=2),
                           col.regions = gray(round( seq(1,0, length=100),2)))
square10.graph









#*****************************************************
#*2. Distribuzioni della posizione dei fotoni rivelati
#*****************************************************

# seleziono solo quei muoni per i quali si è rivelato almeno un fotone
Simulation.photon = Simulation[(!is.na(Simulation$E_fot1)) | (!is.na(Simulation$E_fot2)) | (!is.na(Simulation$E_fot3)) | (!is.na(Simulation$E_fot4)),]

# quantità di fotoni rivelati in totale
vector.energy.tot = c(Simulation.photon$E_fot1,Simulation.photon$E_fot2,Simulation.photon$E_fot3,Simulation.photon$E_fot4) 
length(vector.energy.tot[!is.na(vector.energy.tot)])


Simulation.photon.restricted = Simulation.photon[,c(9,10,11,12,37)]


vox_real.agg = NA 

for(i in 1:nrow(Simulation.photon.restricted))
{
  nr = sum(!is.na(Simulation.photon.restricted[i,-5])) # numero di fotoni rivelati per il muone i-esimo
  if(nr>1)
  {
    vox_real. = c(rep(Simulation.photon.restricted[i,5], (nr-1)))
    vox_real.agg = c(vox_real.agg, vox_real.)
  }
}

vox_real.agg = vox_real.agg[-1]
vox_real.forTable = c(Simulation.photon.restricted$vox_real,vox_real.agg)
vox_real.agg2 = c(1:1000)[!(c(1:1000) %in% vox_real.forTable)] 
vox_real.forTable = c(vox_real.forTable,vox_real.agg2)
vox_real.forTable = as.matrix(table(vox_real.forTable),nrow=1) 
vox_real.forTable[vox_real.agg2] = 0 


coord.x = c("24.5","25.5","26.5","27.5","28.5","29.5","30.5","31.5","32.5","33.5")
coord.y = c("0.5","1.5","2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5")


# grafico quantità di fotoni nei primi 100 voxel (primo strato in z=0.5)
square1 = rbind(vox_real.forTable[1:10],
                vox_real.forTable[11:20],
                vox_real.forTable[21:30],
                vox_real.forTable[31:40],
                vox_real.forTable[41:50],
                vox_real.forTable[51:60],
                vox_real.forTable[61:70],
                vox_real.forTable[71:80],
                vox_real.forTable[81:90],
                vox_real.forTable[91:100])
colnames(square1) <- coord.x                    
rownames(square1) <- coord.y  
square1.graph = levelplot(square1, main=list("z = 0.5", cex=2), sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2))) # gives me 100 shades of gray from white (1) to black (0).
square1.graph

# grafico quantità di fotoni nei successivi 100 voxel (secondo strato in z=1.5)
square2 = rbind(vox_real.forTable[101:110],
                vox_real.forTable[111:120],
                vox_real.forTable[121:130],
                vox_real.forTable[131:140],
                vox_real.forTable[141:150],
                vox_real.forTable[151:160],
                vox_real.forTable[161:170],
                vox_real.forTable[171:180],
                vox_real.forTable[181:190],
                vox_real.forTable[191:200])
colnames(square2) <- coord.x                    
rownames(square2) <- coord.y  
square2.graph = levelplot(square2, main=list("z = 1.5", cex=2),sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square2.graph

# grafico quantità di fotoni nei successivi 100 voxel (terzo strato in z=2.5)
square3 = rbind(vox_real.forTable[201:210],
                vox_real.forTable[211:220],
                vox_real.forTable[221:230],
                vox_real.forTable[231:240],
                vox_real.forTable[241:250],
                vox_real.forTable[251:260],
                vox_real.forTable[261:270],
                vox_real.forTable[271:280],
                vox_real.forTable[281:290],
                vox_real.forTable[291:300])
colnames(square3) <- coord.x                    
rownames(square3) <- coord.y  
square3.graph = levelplot(square3, main=list("z = 2.5", cex=2),sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square3.graph

# grafico quantità di fotoni nei successivi 100 voxel (quarto strato in z=3.5)
square4 = rbind(vox_real.forTable[301:310],
                vox_real.forTable[311:320],
                vox_real.forTable[321:330],
                vox_real.forTable[331:340],
                vox_real.forTable[341:350],
                vox_real.forTable[351:360],
                vox_real.forTable[361:370],
                vox_real.forTable[371:380],
                vox_real.forTable[381:390],
                vox_real.forTable[391:400])
colnames(square4) <- coord.x                    
rownames(square4) <- coord.y  
square4.graph = levelplot(square4, main=list("z = 3.5", cex=2),sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square4.graph

# grafico quantità di fotoni nei successivi 100 voxel (quinto strato in z=4.5)
square5 = rbind(vox_real.forTable[401:410],
                vox_real.forTable[411:420],
                vox_real.forTable[421:430],
                vox_real.forTable[431:440],
                vox_real.forTable[441:450],
                vox_real.forTable[451:460],
                vox_real.forTable[461:470],
                vox_real.forTable[471:480],
                vox_real.forTable[481:490],
                vox_real.forTable[491:500])
colnames(square5) <- coord.x                    
rownames(square5) <- coord.y  
square5.graph = levelplot(square5, main=list("z = 4.5", cex=2),sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square5.graph

# grafico quantità di fotoni nei successivi 100 voxel (sesto strato in z=5.5)
square6 = rbind(vox_real.forTable[501:510],
                vox_real.forTable[511:520],
                vox_real.forTable[521:530],
                vox_real.forTable[531:540],
                vox_real.forTable[541:550],
                vox_real.forTable[551:560],
                vox_real.forTable[561:570],
                vox_real.forTable[571:580],
                vox_real.forTable[581:590],
                vox_real.forTable[591:600])
colnames(square6) <- coord.x                    
rownames(square6) <- coord.y  
square6.graph = levelplot(square6, main=list("z = 5.5", cex=2),sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square6.graph

# grafico quantità di fotoni nei successivi 100 voxel (settimo strato in z=6.5)
square7 = rbind(vox_real.forTable[601:610],
                vox_real.forTable[611:620],
                vox_real.forTable[621:630],
                vox_real.forTable[631:640],
                vox_real.forTable[641:650],
                vox_real.forTable[651:660],
                vox_real.forTable[661:670],
                vox_real.forTable[671:680],
                vox_real.forTable[681:690],
                vox_real.forTable[691:700])
colnames(square7) <- coord.x                    
rownames(square7) <- coord.y  
square7.graph = levelplot(square7, main=list("z = 6.5", cex=2),sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square7.graph

# grafico quantità di fotoni nei successivi 100 voxel (ottavo strato in z=7.5)
square8 = rbind(vox_real.forTable[701:710],
                vox_real.forTable[711:720],
                vox_real.forTable[721:730],
                vox_real.forTable[731:740],
                vox_real.forTable[741:750],
                vox_real.forTable[751:760],
                vox_real.forTable[761:770],
                vox_real.forTable[771:780],
                vox_real.forTable[781:790],
                vox_real.forTable[791:800])
colnames(square8) <- coord.x                    
rownames(square8) <- coord.y  
square8.graph = levelplot(square8, main=list("z = 7.5", cex=2),sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square8.graph

# grafico quantità di fotoni nei successivi 100 voxel (nono strato in z=8.5)
square9 = rbind(vox_real.forTable[801:810],
                vox_real.forTable[811:820],
                vox_real.forTable[821:830],
                vox_real.forTable[831:840],
                vox_real.forTable[841:850],
                vox_real.forTable[851:860],
                vox_real.forTable[861:870],
                vox_real.forTable[871:880],
                vox_real.forTable[881:890],
                vox_real.forTable[891:900])
colnames(square9) <- coord.x                    
rownames(square9) <- coord.y  
square9.graph = levelplot(square9, main=list("z = 8.5", cex=2),sub="provenienza fotoni rivelati",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = gray(round( seq(1,0, length=100),2)))
square9.graph

# grafico quantità di fotoni nei successivi 100 voxel (decimo strato in z=9.5)
square10 = rbind(vox_real.forTable[901:910],
                 vox_real.forTable[911:920],
                 vox_real.forTable[921:930],
                 vox_real.forTable[931:940],
                 vox_real.forTable[941:950],
                 vox_real.forTable[951:960],
                 vox_real.forTable[961:970],
                 vox_real.forTable[971:980],
                 vox_real.forTable[981:990],
                 vox_real.forTable[991:1000])
colnames(square10) <- coord.x                    
rownames(square10) <- coord.y  
square10.graph = levelplot(square10, main=list("z = 9.5", cex=2),sub="provenienza fotoni rivelati",
                           scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                           xlab=list("y", cex=2),
                           ylab=list("x", cex=2),
                           colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                         height=1, width=2),
                           col.regions = gray(round( seq(1,0, length=100),2)))
square10.graph




#*************************************************************************
#*3. Distribuzioni rivelabili dell'energia dei fotoni dai diversi elementi
#*************************************************************************

library(ggplot2)
library(gridExtra)

# energie dei fotoni nel voxel 821 => Z = 26
Energy.Z26.821 = Simulation.photon.restricted[Simulation.photon.restricted$vox_real==821,-5]
Energy.Z26.821 = c(Energy.Z26.821$E_fot1, Energy.Z26.821$E_fot2, Energy.Z26.821$E_fot3, Energy.Z26.821$E_fot4)
Energy.Z26.821 = Energy.Z26.821[!is.na(Energy.Z26.821)]

# energie dei fotoni nel voxel 822 => Z = 26
Energy.Z26.822 = Simulation.photon.restricted[Simulation.photon.restricted$vox_real==822,-5]
Energy.Z26.822 = c(Energy.Z26.822$E_fot1, Energy.Z26.822$E_fot2, Energy.Z26.822$E_fot3, Energy.Z26.822$E_fot4)
Energy.Z26.822 = Energy.Z26.822[!is.na(Energy.Z26.822)]

# energie dei fotoni nel voxel 823 => Z = 26
Energy.Z26.823 = Simulation.photon.restricted[Simulation.photon.restricted$vox_real==823,-5]
Energy.Z26.823 = c(Energy.Z26.823$E_fot1, Energy.Z26.823$E_fot2, Energy.Z26.823$E_fot3, Energy.Z26.823$E_fot4)
Energy.Z26.823 = Energy.Z26.823[!is.na(Energy.Z26.823)]

# energie dei fotoni nel voxel 824 => Z = 26
Energy.Z26.824 = Simulation.photon.restricted[Simulation.photon.restricted$vox_real==824,-5]
Energy.Z26.824 = c(Energy.Z26.824$E_fot1, Energy.Z26.824$E_fot2, Energy.Z26.824$E_fot3, Energy.Z26.824$E_fot4)
Energy.Z26.824 = Energy.Z26.824[!is.na(Energy.Z26.824)]


#----> densità

source("pdf_photon_energy.R")

E.Z26 = c(seq(10^-10, 10, length=2000))
out.Z26 = pdf.Z26(E.Z26)

Z26.real = ggplot(data.frame(x1 = E.Z26, y1 = out.Z26), aes(x = x1)) + 
  geom_line(aes(y = y1, x = x1)) +
  scale_x_continuous(limits = c(10^-10, 10),
                     breaks = round(seq(0, 10, length= 7),2)
  )+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03C6")+
  labs(subtitle ="Z = 26"
       # title = "Probability density function", subtitle ="Z = 26"
  )+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))





# istogrammi

E.Z26 = c(seq(10^-10, 10, length=length(Energy.Z26.821)))
out.Z26 = pdf.Z26(E.Z26)
Z26.hist.821 = ggplot(data.frame(x1 = E.Z26, x2 = Energy.Z26.821, y1 = out.Z26)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "red")+ # funzione di densità teorica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 3.2),breaks = c(0,0.5,1,1.5,2,2.5,3))+
  theme_minimal()+
  labs(title="Z = 26", subtitle = "Voxel 821")+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))


E.Z26 = c(seq(10^-10, 10, length=length(Energy.Z26.822)))
out.Z26 = pdf.Z26(E.Z26)
Z26.hist.822 = ggplot(data.frame(x1 = E.Z26, x2 = Energy.Z26.822, y1 = out.Z26)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "red")+ # funzione di densità teorica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 3.2),breaks = c(0,0.5,1,1.5,2,2.5,3))+
  theme_minimal()+
  labs(title="Z = 26", subtitle = "Voxel 822")+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))


E.Z26 = c(seq(10^-10, 10, length=length(Energy.Z26.823)))
out.Z26 = pdf.Z26(E.Z26)
Z26.hist.823 = ggplot(data.frame(x1 = E.Z26, x2 = Energy.Z26.823, y1 = out.Z26)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "red")+ # funzione di densità teorica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 3.2),breaks = c(0,0.5,1,1.5,2,2.5,3))+
  theme_minimal()+
  labs(title="Z = 26", subtitle = "Voxel 823")+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))



E.Z26 = c(seq(10^-10, 10, length=length(Energy.Z26.824)))
out.Z26 = pdf.Z26(E.Z26)
Z26.hist.824 = ggplot(data.frame(x1 = E.Z26, x2 = Energy.Z26.824, y1 = out.Z26)) + 
  geom_histogram(aes(y=after_stat(density), x=x2),binwidth=0.09)+ # densità empirica
  geom_line(aes(y = y1, x = x1), col= "red")+ # funzione di densità teorica
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Densità")+
  scale_x_continuous(limits = c(10^-10, 10), breaks = c(0,2,4,6,8,10))+
  scale_y_continuous(limits = c(0, 3.2),breaks = c(0,0.5,1,1.5,2,2.5,3))+
  theme_minimal()+
  labs(title="Z = 26", subtitle = "Voxel 824")+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))


grid.arrange(Z26.hist.821,Z26.hist.822,ncol=2)
grid.arrange(Z26.hist.823,Z26.hist.824,ncol=2)
grid.arrange(Z26.hist.821,Z26.hist.822,Z26.hist.823,Z26.hist.824,ncol=4)



#----> funzione di ripartizione

source("cdf_photon_energy.R")


# confronto funzioni di ripartizione empiriche e teoriche

E.Z26 = c(seq(10^-10, 10, length=length(Energy.Z26.821)))
out.cdf.Z26 = cdf.Z26(E.Z26) # calcolo la cdf teorica
cdf.Z26.821 = ggplot(data.frame(x1 = E.Z26, y1 = out.cdf.Z26, x2 = Energy.Z26.821)) + 
  stat_ecdf(aes(x=x2))+ # cdf empirica
  geom_line(aes(y = y1, x = x1),col="red")+ # cdf teorica
  scale_x_continuous(limits = c(10^-10, 10),
                     breaks = round(seq(0, 10, length= 5),2)
  )+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Funzione di ripartizione")+
  labs(title="Z = 26", subtitle = "Voxel 821")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    title =  element_text(size = 20))



E.Z26 = c(seq(10^-10, 10, length=length(Energy.Z26.822)))
out.cdf.Z26 = cdf.Z26(E.Z26) # calcolo la cdf teorica
cdf.Z26.822 = ggplot(data.frame(x1 = E.Z26, y1 = out.cdf.Z26, x2 = Energy.Z26.822)) + 
  stat_ecdf(aes(x=x2))+ # cdf empirica
  geom_line(aes(y = y1, x = x1),col="red")+ # cdf teorica
  scale_x_continuous(limits = c(10^-10, 10),
                     breaks = round(seq(0, 10, length= 5),2)
  )+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Funzione di ripartizione")+
  labs(title="Z = 26", subtitle = "Voxel 822")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    title =  element_text(size = 20))


E.Z26 = c(seq(10^-10, 10, length=length(Energy.Z26.823)))
out.cdf.Z26 = cdf.Z26(E.Z26) # calcolo la cdf teorica
cdf.Z26.823 = ggplot(data.frame(x1 = E.Z26, y1 = out.cdf.Z26, x2 = Energy.Z26.823)) + 
  stat_ecdf(aes(x=x2))+ # cdf empirica
  geom_line(aes(y = y1, x = x1),col="red")+ # cdf teorica
  scale_x_continuous(limits = c(10^-10, 10),
                     breaks = round(seq(0, 10, length= 5),2)
  )+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Funzione di ripartizione")+
  labs(title="Z = 26", subtitle = "Voxel 823")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    title =  element_text(size = 20))


E.Z26 = c(seq(10^-10, 10, length=length(Energy.Z26.824)))
out.cdf.Z26 = cdf.Z26(E.Z26) # calcolo la cdf teorica
cdf.Z26.824 = ggplot(data.frame(x1 = E.Z26, y1 = out.cdf.Z26, x2 = Energy.Z26.824)) + 
  stat_ecdf(aes(x=x2))+ # cdf empirica
  geom_line(aes(y = y1, x = x1),col="red")+ # cdf teorica
  scale_x_continuous(limits = c(10^-10, 10),
                     breaks = round(seq(0, 10, length= 5),2)
  )+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("Funzione di ripartizione")+
  labs(title="Z = 26", subtitle = "Voxel 824")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    title =  element_text(size = 20))

grid.arrange(cdf.Z26.821,cdf.Z26.822,cdf.Z26.823,cdf.Z26.824,ncol=4)



#*********************************
#* Composizione atomica del volume
#*********************************

coord.x = c("24.5","25.5","26.5","27.5","28.5","29.5","30.5","31.5","32.5","33.5")
coord.y = c("0.5","1.5","2.5","3.5","4.5","5.5","6.5","7.5","8.5","9.5")


# LEGEND:
# dark green -> Z = 6
# light green -> Z = 14
# white -> Z = 26


# grafico quantità di fotoni nei primi 100 voxel (primo strato in z=0.5)
square1 = rbind(cube3$Z[1:10],
                cube3$Z[11:20],
                cube3$Z[21:30],
                cube3$Z[31:40],
                cube3$Z[41:50],
                cube3$Z[51:60],
                cube3$Z[61:70],
                cube3$Z[71:80],
                cube3$Z[81:90],
                cube3$Z[91:100])
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

# grafico quantità di fotoni nei successivi 100 voxel (secondo strato in z=1.5)
square2 = rbind(cube3$Z[101:110],
                cube3$Z[111:120],
                cube3$Z[121:130],
                cube3$Z[131:140],
                cube3$Z[141:150],
                cube3$Z[151:160],
                cube3$Z[161:170],
                cube3$Z[171:180],
                cube3$Z[181:190],
                cube3$Z[191:200])
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

# grafico quantità di fotoni nei successivi 100 voxel (terzo strato in z=2.5)
square3 = rbind(cube3$Z[201:210],
                cube3$Z[211:220],
                cube3$Z[221:230],
                cube3$Z[231:240],
                cube3$Z[241:250],
                cube3$Z[251:260],
                cube3$Z[261:270],
                cube3$Z[271:280],
                cube3$Z[281:290],
                cube3$Z[291:300])
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

# grafico quantità di fotoni nei successivi 100 voxel (quarto strato in z=3.5)
square4 = rbind(cube3$Z[301:310],
                cube3$Z[311:320],
                cube3$Z[321:330],
                cube3$Z[331:340],
                cube3$Z[341:350],
                cube3$Z[351:360],
                cube3$Z[361:370],
                cube3$Z[371:380],
                cube3$Z[381:390],
                cube3$Z[391:400])
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

# grafico quantità di fotoni nei successivi 100 voxel (quinto strato in z=4.5)
square5 = rbind(cube3$Z[401:410],
                cube3$Z[411:420],
                cube3$Z[421:430],
                cube3$Z[431:440],
                cube3$Z[441:450],
                cube3$Z[451:460],
                cube3$Z[461:470],
                cube3$Z[471:480],
                cube3$Z[481:490],
                cube3$Z[491:500])
colnames(square5) <- coord.x                    
rownames(square5) <- coord.y  
square5.graph = levelplot(square5, main=list("z = 4.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square5.graph

# grafico quantità di fotoni nei successivi 100 voxel (sesto strato in z=5.5)
square6 = rbind(cube3$Z[501:510],
                cube3$Z[511:520],
                cube3$Z[521:530],
                cube3$Z[531:540],
                cube3$Z[541:550],
                cube3$Z[551:560],
                cube3$Z[561:570],
                cube3$Z[571:580],
                cube3$Z[581:590],
                cube3$Z[591:600])
colnames(square6) <- coord.x                    
rownames(square6) <- coord.y  
square6.graph = levelplot(square6, main=list("z = 5.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square6.graph

# grafico quantità di fotoni nei successivi 100 voxel (settimo strato in z=6.5)
square7 = rbind(cube3$Z[601:610],
                cube3$Z[611:620],
                cube3$Z[621:630],
                cube3$Z[631:640],
                cube3$Z[641:650],
                cube3$Z[651:660],
                cube3$Z[661:670],
                cube3$Z[671:680],
                cube3$Z[681:690],
                cube3$Z[691:700])
colnames(square7) <- coord.x                    
rownames(square7) <- coord.y  
square7.graph = levelplot(square7, main=list("z = 6.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square7.graph

# grafico quantità di fotoni nei successivi 100 voxel (ottavo strato in z=7.5)
square8 = rbind(cube3$Z[701:710],
                cube3$Z[711:720],
                cube3$Z[721:730],
                cube3$Z[731:740],
                cube3$Z[741:750],
                cube3$Z[751:760],
                cube3$Z[761:770],
                cube3$Z[771:780],
                cube3$Z[781:790],
                cube3$Z[791:800])
colnames(square8) <- coord.x                    
rownames(square8) <- coord.y  
square8.graph = levelplot(square8, main=list("z = 7.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square8.graph

# grafico quantità di fotoni nei successivi 100 voxel (nono strato in z=8.5)
square9 = rbind(cube3$Z[801:810],
                cube3$Z[811:820],
                cube3$Z[821:830],
                cube3$Z[831:840],
                cube3$Z[841:850],
                cube3$Z[851:860],
                cube3$Z[861:870],
                cube3$Z[871:880],
                cube3$Z[881:890],
                cube3$Z[891:900])
colnames(square9) <- coord.x                    
rownames(square9) <- coord.y  
square9.graph = levelplot(square9, main=list("z = 8.5", cex=2),sub="composizione atomica",
                          scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                          xlab=list("y", cex=2),
                          ylab=list("x", cex=2),
                          colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                        height=1, width=2),
                          col.regions = terrain.colors(n=100,alpha = 0.4)) 
square9.graph

# grafico quantità di fotoni nei successivi 100 voxel (decimo strato in z=9.5)
square10 = rbind(cube3$Z[901:910],
                 cube3$Z[911:920],
                 cube3$Z[921:930],
                 cube3$Z[931:940],
                 cube3$Z[941:950],
                 cube3$Z[951:960],
                 cube3$Z[961:970],
                 cube3$Z[971:980],
                 cube3$Z[981:990],
                 cube3$Z[991:1000])
colnames(square10) <- coord.x                    
rownames(square10) <- coord.y  
square10.graph = levelplot(square10, main=list("z = 9.5", cex=2),sub="composizione atomica",
                           scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                           xlab=list("y", cex=2),
                           ylab=list("x", cex=2),
                           colorkey=list(labels=list(cex=1.5, font=2, col="black"),
                                         height=1, width=2),
                           col.regions = terrain.colors(n=100,alpha = 0.4)) 
square10.graph





#******************************************************
#* Distribuzione delle energie complessivamente rivelate
#******************************************************

E.tot = c(Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,1], Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,2],
          Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,3], Simulation.photon[, c("E_fot1", "E_fot2", "E_fot3", "E_fot4")][,4])
E.tot = E.tot[!is.na(E.tot)]


# istogramma di tutte le energie rivelate complessivamente in tutti i pixel

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
