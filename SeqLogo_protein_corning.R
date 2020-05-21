library(ggpubr)
#install.packages("ggseqlogo")
library(tidyverse)
library(ggseqlogo)

CP<-read.table("~/Descargas/C_P_mirnas.txt",header=T)
names(CP) = gsub(pattern = ".tailed_reads", replacement = "", x = names(CP))

names(CP)


m<-list()
m

for( i in c("R-380-3","R-10b-5p","R-369-3","R-10a-5","R-451-5","R-127-3","R-423-5","R-100-5","R-148a-5","R-320a-3")) {

    CP%>%
    filter(str_detect(CP$miR.name,i))-> TCV
    
  #####sumamos aquellos miRNAs que tienen a en primerlugar y ...
  
  for (a in c("A","T","C","G")) {
    TCV[paste("II",a, sep = "_")] <- rowSums(TCV[,str_detect(names(TCV),as.character(paste(".",a,sep = "")))])  
    
  }
  
    for (a in c("A","T","C","G")) {
    TCV[paste("I",a, sep = "_")] <- rowSums(TCV[,startsWith(names(TCV),a)])  
  }

  
TCV_s<-TCV[,c(3,27:34)]
  
z<-colSums(TCV_s)
  
  TCVx<-matrix(c(z[2:9]), nrow=4,ncol=2,byrow = F)
  row.names(TCVx)<-c("A","U","C","G")
  
  
  x<-z[1]-colSums(TCVx)
  
  TCVx<-rbind(TCVx,x)
  row.names(TCVx)<-c("A","U","C","G","NA")
  
  TCVx<-TCVx[,c(2,1)]
  
  ggseqlogo(TCVx,method="prob" )+
    labs(title = i)->m[[i]]
}

#library(ggpubr)
ggarrange(plotlist=m, 
          labels = c("A", "B","C","D","E","F","G","H","i","J"))

######################################333
########################################3
####################
CP%>%
  filter(str_detect(CP$miR.name,"R-380-3"))-> TCV

names(TCV)

for (a in c("A","T","C","G")) {
  TCV[paste("II",a, sep = "_")] <- rowSums(TCV[,str_detect(names(TCV),as.character(paste(".",a,sep = "")))])  
  
}

for (a in c("A","T","C","G")) {
  TCV[paste("I",a, sep = "_")] <- rowSums(TCV[,startsWith(names(TCV),a)])  
}


TCV_s<-TCV[,c(3,27:34)]
TCV_s
z<-colSums(TCV_s)
z
TCVx<-matrix(c(z[2:9]), nrow=4,ncol=2,byrow = F)
row.names(TCVx)<-c("A","U","C","G")
TCVx

x<-z[1]-colSums(TCVx)

TCVx<-rbind(TCVx,x)
row.names(TCVx)<-c("A","U","C","G","NA")


ggseqlogo(TCVx,method="prob" )+
  labs(title = i)->m[[i]]