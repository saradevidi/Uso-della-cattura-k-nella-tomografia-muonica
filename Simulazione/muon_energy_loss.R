
#+++++++++++++++++++++++++++++++++++++++++++++++++
#+ CALCOLO LA PERDITA MEDIA DI ENERGIA DEL MUONE
#+++++++++++++++++++++++++++++++++++++++++++++++++


loss.energy = function(Zdensity, Y = 1)
{ # funzione che dati:
  #   - Zdensity: matrice con 2 colonne (nella prima inserisco i numeri atomici Z, nella seconda inserisco la corrispondente densità in gr/cm^3)
  #   - Y: distanza attraversata in cm
  # calcola la perdita di energia del muone nell'attraversare Y cm di materiali di differenti numeri atomici Z.
  perditaEnergia.Y = ( 2.35 - 0.28 * log(Zdensity[,1]) ) * Zdensity[,2] * Y
  out = cbind(Zdensity, perditaEnergia.Y)
  colnames(out) = c("Z", "density", "loss.energy.Ycm")
  out = as.data.frame(out)
  out
}




