############ Quality plot of Fastqc results of Streptanthus sequences after trimming ##############

## Load library
library(ggplot2)

## Read data 
dat <- read.csv("../strep_demulti_fastqs/4_fastqc_trim_qualitycharts/1_summary_trim.csv")

## Plot
head (dat)
dat2<-table(dat$Wrong)
class(dat2)

dat2 <- as.data.frame(dat2)
dat2 <- dat2[-c(1),]

##Plot
ggplot(data=dat2, aes(x=Var1, y= Freq)) + geom_bar(stat="identity")
myplot<-ggplot(data=dat2, aes(x=Var1, y= Freq))
myplot + geom_bar(stat="identity")

#Change color
x <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) + geom_bar(stat="identity")
ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1))

##Remove grey background
y <- ggplot(data=dat2, aes(x=Var1, y= Freq, color=Var1)) +
  geom_bar(stat="identity", aes(color=Var1, fill=Var1)) + 
  theme_bw()

##Change axis names
z <- y + ylab("Frequency") + xlab("Problem")
z  

##Remove names of each bar
a <- z + theme(axis.text.x=element_blank())
a


