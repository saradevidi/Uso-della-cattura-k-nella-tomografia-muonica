#--------------------------------------------
# grafici per le pdf delle energie dei fotoni
#--------------------------------------------

library(ggplot2)
library(gridExtra)


source("pdf_photon_energy.R")
source("absorption_length.R")

E.Z6 = c(seq(10^-10, 10, length=50000))

out.Z6 = pdf.Z6(E.Z6)
out.Z14 = pdf.Z14(E.Z14)
out.Z26 = pdf.Z26(E.Z26)

ggplot(data.frame(x1 = E.Z6, y1 = out.Z6, y2 = out.Z14, y3 = out.Z26), aes(x = x1)) + 
  geom_line(aes(y = out.Z6, x = x1), color="black", linewidth=1) +
  geom_line(aes(y = out.Z14, x = x1), color="black", linewidth=1)+
  geom_line(aes(y = out.Z26, x = x1), color="black", linewidth=1)+
  scale_x_continuous(limits = c(10^-10, 10),
                     breaks = round(seq(0, 10, length= 7),2)
  )+
  scale_y_continuous(limits = c(10^-10, 2.3))+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03C6")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))




#--------------------------------
# grafici per l'absorption length
#--------------------------------


## z=26

uno = seq(dataForGraph26$Energy.Mev[1], dataForGraph26$Energy.Mev[8], length=50)
due = c(seq(dataForGraph26$Energy.Mev[9]+0.0001, dataForGraph26$Energy.Mev[19], length=100),
        seq(dataForGraph26$Energy.Mev[19], dataForGraph26$Energy.Mev[36], length=100))

y1 = funTotVec.Z26(c(uno))
y2 = funTotVec.Z26(c(due))

Z26 = ggplot(data.frame(x1 = c(uno), x2 = c(due), y1 = y1, y2 = y2), aes(x =  c(x1, x2))) + 
  geom_line(aes(y = y1, x = x1)) +
  geom_line(aes(y = y2, x = x2))+
  geom_point(data = dataForGraph26, aes(x = Energy.Mev, y = lambda),alpha = 1/2, colour = "darkgray", size = 2)+
  scale_x_continuous(limits = c(dataForGraph26$Energy.Mev[1], dataForGraph26$Energy.Mev[36]),
                     breaks = round(seq(dataForGraph26$Energy.Mev[1], dataForGraph26$Energy.Mev[36], length= 5),3)
  )+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03BB")+
  labs(subtitle ="Z = 26")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))


# zoom
Z26Zoom = ggplot(data.frame(x1 = c(uno), x2 = c(due), y1 = y1, y2 = y2), aes(x =  c(x1, x2))) + 
  geom_line(aes(y = y1, x = x1)) +
  geom_line(aes(y = y2, x = x2))+
  geom_point(data = dataForGraph26, aes(x = Energy.Mev, y = lambda),alpha = 1/2, colour = "darkgray", size = 2)+
  scale_x_continuous(limits = c(dataForGraph26$Energy.Mev[1], dataForGraph26$Energy.Mev[36])
                     ,breaks = round(seq(dataForGraph26$Energy.Mev[1],  0.025, length= 4),3)
  )+
  coord_cartesian(xlim=c(0.001,0.025), ylim=c(0,0.05))+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03BB")+
  labs(title="Zoom")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))

grid.arrange(Z26,Z26Zoom,ncol=2)





## Z=14

uno.Z14 = seq(dataForGraph14$Energy.Mev[1], dataForGraph14$Energy.Mev[3], length=50)
due.Z14 = c(seq(dataForGraph14$Energy.Mev[4]+0.0001, dataForGraph14$Energy.Mev[19], length=100),
            seq(dataForGraph14$Energy.Mev[19], dataForGraph14$Energy.Mev[36], length=100))



y1.Z14 = funTotVec.Z14(c(uno.Z14))
y2.Z14 = funTotVec.Z14(c(due.Z14))

Z14 = ggplot(data.frame(x1 = c(uno.Z14), x2 = c(due.Z14), y1 = y1.Z14, y2 = y2.Z14), aes(x = c(x1, x2))) + 
  geom_line(aes(y = y1, x = x1)) +
  geom_line(aes(y = y2, x = x2))+
  geom_point(data = dataForGraph14, aes(x = Energy.Mev, y = lambda),alpha = 1/2, colour = "darkgray", size = 2)+
  scale_x_continuous(limits = c(dataForGraph14$Energy.Mev[1], dataForGraph14$Energy.Mev[36]),
                     breaks = round(seq(dataForGraph26$Energy.Mev[1], dataForGraph26$Energy.Mev[36], length= 5),3)
  )+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03BB")+
  labs(subtitle ="Z = 14")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))


# zoom
Z14Zoom = ggplot(data.frame(x1 = c(uno.Z14), x2 = c(due.Z14), y1 = y1.Z14, y2 = y2.Z14), aes(x =  c(x1, x2))) + 
  geom_line(aes(y = y1, x = x1)) +
  geom_line(aes(y = y2, x = x2))+
  geom_point(data = dataForGraph14, aes(x = Energy.Mev, y = lambda),alpha = 1/2, colour = "darkgray", size = 2)+
  scale_x_continuous(limits = c(dataForGraph14$Energy.Mev[1], dataForGraph14$Energy.Mev[36])
                     ,breaks = round(seq(dataForGraph26$Energy.Mev[1],  0.025, length= 4),3)
  )+
  coord_cartesian(xlim=c(0.001,0.025), ylim=c(0,0.05))+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03BB")+
  labs(title="Zoom")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))

grid.arrange(Z14,Z14Zoom,ncol=2)



## Z=6

uno.Z6 = c(seq(dataForGraph6$Energy.Mev[1], dataForGraph6$Energy.Mev[15], length=200),
           seq(dataForGraph6$Energy.Mev[15], dataForGraph6$Energy.Mev[34], length=100))



y1.Z6 = funTotVec.Z6(c(uno.Z6))


Z6 = ggplot(data.frame(x1 = c(uno.Z6), y1 = y1.Z6), aes(x = c(x1))) + 
  geom_line(aes(y = y1, x = x1)) +
  geom_point(data = dataForGraph6, aes(x = Energy.Mev, y = lambda),alpha = 1/2, colour = "darkgray", size = 2)+
  scale_x_continuous(limits = c(dataForGraph6$Energy.Mev[1], dataForGraph6$Energy.Mev[34]),
                     breaks = round(seq(dataForGraph26$Energy.Mev[1], dataForGraph26$Energy.Mev[36], length= 5),3)
  )+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03BB")+
  labs(subtitle ="Z = 6")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))


# zoom
Z6Zoom =  ggplot(data.frame(x1 = c(uno.Z6), y1 = y1.Z6), aes(x = c(x1))) + 
  geom_line(aes(y = y1, x = x1)) +
  geom_point(data = dataForGraph6, aes(x = Energy.Mev, y = lambda),alpha = 1/2, colour = "darkgray", size = 2)+
  scale_x_continuous(limits = c(dataForGraph6$Energy.Mev[1], dataForGraph6$Energy.Mev[34])
                     ,breaks = round(seq(dataForGraph26$Energy.Mev[1],  0.025, length= 4),3)
  )+
  coord_cartesian(xlim=c(0.001,0.025), ylim=c(0,0.05))+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03BB")+
  labs(title="Zoom")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20)) 

grid.arrange(Z6,Z6Zoom,ncol=2)





### Confronto funzioni stimate e funzioni originali ottenute unendo i punti dei valori 
### osservati

graph.original =  ggplot(data.frame(x1 = c(dataForGraph6$Energy.Mev,rep(dataForGraph6$Energy.Mev[34],2)), x2 = dataForGraph14$Energy.Mev, x3 = dataForGraph26$Energy.Mev,
                                    y1 = c(dataForGraph6$lambda,rep(dataForGraph6$lambda[34],2)), y2 = dataForGraph14$lambda, y3 = dataForGraph26$lambda), aes(x = c(x1, x2, x3))) + 
  geom_line(aes(y = y1, x = x1), color="darkgreen", linewidth=0.8) + # Z = 6
  geom_line(aes(y = y2, x = x2), color="darkblue", linewidth=0.8) + # Z = 14
  geom_line(aes(y = y3, x = x3), color="darkred", linewidth=0.8) + # Z = 26
  scale_x_continuous(limits = c(dataForGraph6$Energy.Mev[1], dataForGraph6$Energy.Mev[34]),
                     breaks = round(seq(dataForGraph26$Energy.Mev[1], dataForGraph6$Energy.Mev[34], length= 6),3)
  )+
  coord_cartesian(xlim=c(0.001,10), ylim=c(0,50))+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03BB")+
  theme_minimal()+
  labs(title="Originale")+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))


# Z=26
uno = seq(0, dataForGraph26$Energy.Mev[8], length=50)
due = c(seq(dataForGraph26$Energy.Mev[9]+0.0001, dataForGraph26$Energy.Mev[19], length=100),
        seq(dataForGraph26$Energy.Mev[19], 11, length=100))
y1 = funTotVec.Z26(c(uno))
y2 = funTotVec.Z26(c(due))
## Z=14
uno.Z14 = seq(0, dataForGraph14$Energy.Mev[3], length=50)
due.Z14 = c(seq(dataForGraph14$Energy.Mev[4]+0.0001, dataForGraph14$Energy.Mev[19], length=100),
            seq(dataForGraph14$Energy.Mev[19], 11, length=100))
y1.Z14 = funTotVec.Z14(c(uno.Z14))
y2.Z14 = funTotVec.Z14(c(due.Z14))
## Z=6
uno.Z6 = c(seq(0, dataForGraph6$Energy.Mev[15], length=200),
           seq(dataForGraph6$Energy.Mev[15], 11, length=100))
y1.Z6 = funTotVec.Z6(c(uno.Z6))

graph.function =  ggplot(data.frame(x1Z6 = uno.Z6, x1Z14 = c(uno.Z14, rep(uno.Z14[50],250)), x2Z14 = c(due.Z14, rep(due.Z14[200],100)), x1Z26 = c(uno, rep(uno[50],250)), x2Z26= c(due, rep(due[200], 100)),
                                    y1Z6 = y1.Z6,  y1Z14 = c(y1.Z14, rep(y1.Z14[50],250)),   y2Z14 = c(y2.Z14, rep(y2.Z14[200],100)),   y1Z26 = c(y1, rep(y1[50], 250)),  y2Z26 = c(y2, rep(y2[200], 100))), 
                         aes(x = c(x1Z6))) + 
  geom_line(aes(y = y1Z6, x = x1Z6), color="darkgreen", linewidth=0.8) + # Z = 6
  geom_line(aes(y = y1Z14, x = x1Z14), color="darkblue", linewidth=0.8) + # Z = 14
  geom_line(aes(y = y2Z14, x = x2Z14), color="darkblue", linewidth=0.8) +
  geom_line(aes(y = y1Z26, x = x1Z26), color="darkred", linewidth=0.8) + # Z = 26
  geom_line(aes(y = y2Z26, x = x2Z26), color="darkred", linewidth=0.8) +
  scale_x_continuous(limits = c(0, 11),
                     breaks = round(seq(dataForGraph26$Energy.Mev[1], dataForGraph6$Energy.Mev[34], length= 6),3)
  )+
  coord_cartesian(xlim=c(0.001,10), ylim=c(0,50))+
  xlab(expression(paste("E", gamma, "(MeV)")))+
  ylab("\u03BB")+
  theme_minimal()+
  labs(title="Modellazione")+
  theme(
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    title =  element_text(size = 20))


grid.arrange(graph.original,graph.function,ncol=2)






