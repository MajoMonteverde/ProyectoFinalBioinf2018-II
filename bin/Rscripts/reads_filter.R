##################### RAW READS --> FILTER #########################


## Load library
library(ggplot2)


## Load data
dat <- read.csv("reads_filter.csv")
head(dat)


##Plot
x <- ggplot(data=dat, aes(x=ids, y=percentage)) + geom_point(show.legend = F)

##Change axis names
z <- x + ylab("Percentage") + xlab("Individual")
z  


##Remove grey background
b <-  z + theme(axis.text.x=element_blank()) + theme_bw()
b


##Remove names of each bar
a <- b + theme(axis.text.x=element_blank())
a



##################### Coverage #####################

## load data
dat2 <- read.csv("cobertura.csv")
head(dat2)

##Plot
x <- ggplot(data=dat2, aes(x=sample, y=sample_coverage)) + geom_point(show.legend = F)
x

##Changes axis names
z <- x + ylab("Sample Coverage") + xlab("Individual")
z  


## Remove grey background
b <-  z + theme(axis.text.x=element_blank()) + theme_bw()
b


## Remove name of each bar
a <- b + theme(axis.text.x=element_blank())
a
