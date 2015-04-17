require(plyr)
require(ggplot2)
require(sda)
require(mclust)
require(grid)
require(reshape2)

# Import data of automatically segmented cells
data.cells<-read.csv('data/MyExpt_FilteredCells.csv', header=TRUE, sep=',')
# Import data of expanded nuclei segmentation
data.expcells<-read.csv('data/MyExpt_FilteredExpCells.csv', header=TRUE, sep=',')
# select only texture features
fts.cells<-c(6,76:102,105:157)
fts.expcells<-c(6,28:54,57:109)
# separate labels and numerical data
labels.cells<-data.cells$Metadata_gene
labels.expcells<-data.expcells$Metadata_gene
data1.cells<-data.cells[,fts.cells]
data1.expcells<-data.expcells[,fts.expcells]

# join y1,y2 and y3 from expcells to the other dataset
selection.cells<-c("p1c","p2c","p3c","y4c","y5c","y6c","y7c")
selection.expcells<-c("y1c","y2c","y3c")

# join them in the same dataset
data<-rbind(data1.cells[data1.cells$Metadata_gene %in% selection.cells,],
      data1.expcells[data1.expcells$Metadata_gene %in% selection.expcells,])

labels<-data$Metadata_gene
# mean number of cells per gene
mean(count(labels)$freq)

# helper functions to analyse data
do_pca<-function(x,labels){
  # scale to zero mean unit variance
  d<-scale(x)
  
  # do PCA
  d.pca<-prcomp(d.selected)
  # convert to dataframe
  d1<-as.data.frame(d.pca$x)
  # add labels
  d1$labels<-labels
  return(d1)
}

# calculate mean points
means<-function(x,labels,lda){
  # first do pca
  d<-do_pca(x,labels,lda)
  # then calculate means
  m<-ddply(d,.(labels),summarise,PC1=mean(PC1),PC2=mean(PC2),PC3=mean(PC3),PC4=mean(PC4))
  return(m)
}

# plot PCA plot
plot_pca<-function(x,labels,fts){
  text_pt<-9
  text_mm<-9*0.35
  
  p<-ggplot(do_pca(x,labels), aes_string(x=fts[1], y=fts[2]))+
    geom_point(aes(colour=labels), shape=16, size=1)+
    geom_point(data=means(x,labels), aes_string(x=fts[1], y=fts[2], fill="labels"), size=2, shape=23,colour="black")
  p+theme_bw()+coord_fixed()+theme(text=element_text(size=text_pt),legend.position="none")
}

# Make a sequence of 1 to 4
pcs<-t(combn(seq(1,4),2))

# Plot all combinations of 4 PCAs and save them
for(i in seq(1,nrow(pcs))){
  
  plot_pca(data[,-1],labels,fts=paste("PC",pcs[i,],sep=""))
  ggsave(paste("PC",paste(pcs[i,],collapse=""),".pdf",sep=""),width=75,height=75,units = "mm")

}


# Do model bsed clustering
fit<-Mclust(data[,-1])

# join cluster predictions with original classes
classes<-data.frame(class=fit$classification,label=data[,1])
counts<-count(classes,.(label,class))

lbls<-levels(labels)

ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# plot piecharts
for(i in seq(1,10)) {
  
  lbl<-lbls[i]

  # plots
  text_pt<-9
  text_mm<-9*0.35
  
  pie <- ggplot(counts[counts$label==lbl,], aes(x = "",y=freq, fill = factor(class))) +
    geom_bar(width = 1,stat="identity")
  pie + theme_bw() + coord_polar(theta = "y") + scale_fill_manual(name="Class",
                                                                  values = ggplotColours(n=3)) +
    labs(x=NULL,y=NULL) +
    theme(text=element_text(size=text_pt),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position="none")
  
  ggsave(paste(lbl,"-piechart.pdf",sep=""),width=30,height=30,units = "mm")

}

# Plot the BIC plot
data.bic<-melt(fit$BIC)

text_pt<-9
text_mm<-9*0.35

bic<-ggplot(data.bic, aes(x=Var1, y=value, fill=Var2))+ geom_line(aes(colour=Var2))+
  geom_point(aes(colour=Var2), shape=22, size=2,colour="black")
bic+theme_bw() +theme(text=element_text(size=text_pt), plot.margin=unit(c(1,0,0,0),"mm")) +
  labs(x="Number of clusters", y="BIC")+scale_color_discrete(name="Model")+scale_fill_discrete(name="Model")+
  scale_x_continuous(limits=c(1,9),breaks=seq(0, 9, 1))

ggsave("bic.pdf",width=150,height=100,unit="mm")
ggsave("bic.tiff",width=150,height=100,unit="mm")
