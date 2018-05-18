############ Gráfica calidad de seequencias Steptanthus después del demultiplexeo ##############

## Cargar libreria
library(ggplot2)

## Leer base de datos
dat <- read.csv("../strep_demulti_fastqs/2_fastqc_qualitycharts/1_summary.csv")

## Graficar
head (dat)
dat2<-table(dat$Wrong)
class(dat2)

dat2 <- as.data.frame(dat2)
dat2 <- dat2[-c(1),]

##Graficar
ggplot(data=dat2, aes(x=Var1, y= Freq)) + geom_bar(stat="identity")
myplot<-ggplot(data=dat2, aes(x=Var1, y= Freq))
myplot + geom_bar(stat="identity")

#aes cambiar color
x <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) + geom_bar(stat="identity")
ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1))

##quitar fondo gris
y <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1)) + 
  theme_bw()

##cambiar nombre de ejes
z <- y + ylab("Frequency") + xlab("Problem")
z  

##quitar los nombres que tiene para cada barra
a <- z + theme(axis.text.x=element_blank())
a

sessionInfo()
