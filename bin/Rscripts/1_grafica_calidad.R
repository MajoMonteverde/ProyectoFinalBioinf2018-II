############ Quality plot of Fastqc results of Streptanthus sequences after demultiplexing ##############

## Load library
library(ggplot2)

## Read csv summary
dat <- read.csv("../strep_demulti_fastqs/2_fastqc_qualitycharts/1_summary.csv")

## Plot
head (dat)
dat2<-table(dat$Wrong)
class(dat2)

dat2 <- as.data.frame(dat2)
dat2 <- dat2[-c(1),]

## Plot
ggplot(data=dat2, aes(x=Var1, y= Freq)) + geom_bar(stat="identity")
myplot<-ggplot(data=dat2, aes(x=Var1, y= Freq))
myplot + geom_bar(stat="identity")

# Change color
x <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) + geom_bar(stat="identity")
ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1))

## Remove grey background
y <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1)) + 
  theme_bw()

## Change axis names
z <- y + ylab("Frequency") + xlab("Problem")
z  

## Rmeove names for each bar
a <- z + theme(axis.text.x=element_blank())
a
