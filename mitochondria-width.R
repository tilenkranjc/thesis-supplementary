require(plyr)
require(ggplot2)
require(grid)

# load data
data<-data.frame()
for (i in list.files("width",pattern="p1.*.tab",full.names = TRUE)) {
  data<-rbind(data,cbind(read.table(i),data.frame(trt="siPRAF1")))
}

for (i in list.files("width",pattern="neg.*.tab",full.names = TRUE)) {
  data<-rbind(data,cbind(read.table(i),data.frame(trt="siNEG")))
}

# Scale in meters per pixel, acquired from image metadata
scale<-1.2698395083648987E-07
# in micrometers per pixel
scale1<-scale*10E6


means1<-ddply(data,.(trt),summarise,mean=mean(MinFeret),
             median=median(MinFeret),
             sd=sd(MinFeret),
             se=sd(MinFeret)/sqrt(length(MinFeret)))

t.test(MinFeret~trt,data=data)

# Remove outliers with boxplot
# split on treatment to remove outliers with boxplot
a<-split(data,data$trt)

# create a new dataframe without outliers
data1<-data.frame()

for (j in names(a)) {
  d<-a[[j]]
  b<-boxplot(d$MinFeret)$stats[c(1,5),]
  data1<-rbind(data1,d[d["MinFeret"]>b[1] & d["MinFeret"]<b[2],])
}

# calculate new means
means<-ddply(data1,.(trt),summarise,mean=mean(MinFeret),
             median=median(MinFeret),
             sd=sd(MinFeret),
             se=sd(MinFeret)/sqrt(length(MinFeret)))

# and calculate means in micrometers
means.um<-ddply(data1,.(trt),summarise,mean=mean(MinFeret)*scale1,
             median=median(MinFeret)*scale1,
             sd=sd(MinFeret)*scale1,
             se=sd(MinFeret)*scale1/sqrt(length(MinFeret)))
# correct factor names
means.um$trt<-factor(means.um$trt,levels=c("siNEG","siPRAF1"),ordered=TRUE)

# perform a t-test
t.test(MinFeret~trt,data=data1)

# get star offset and create labels per group
# Symbol  Meaning
# ns  P > 0.05
# *  P ≤ 0.05
# **  P ≤ 0.01
# ***  P ≤ 0.001
# ****   P ≤ 0.0001

# set limits of plot and star significance
max.point<-max(means.um$mean)
limits<-c(0,max.point*1.2)
offset<-limits[2]*0.05
means.um$stars<-c("****","")

text_points<-9
text_mm<-9*0.35

# plot
plot<-ggplot(means.um, aes(x=factor(trt), y=mean))+geom_bar(stat="identity",fill='grey30', width=0.5)+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position='dodge',colour='black', width=0.1) +
  geom_text(aes(y=mean+offset,label=stars),size=text_mm)
plot + theme_bw() + labs(title="Thickness of mitochondria\nin PRAF1 depleted cells",
                         x="Treatment", y=expression(paste("Mean thickness in ",mu,"m"))) + 
  scale_y_continuous(limits=limits, expand = c(0,0))+
  theme(text=element_text(size=text_points),plot.margin=unit(c(1,1,0,0),"mm"))

# save
ggsave(file="rnai-width.pdf",width=80,height=100,units="mm")
ggsave(file="rnai-width.tiff",width=80,height=100,units="mm")
