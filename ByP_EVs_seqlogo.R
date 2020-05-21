BV<-read.table("ByP_EVs_mirnas.txt",header=T)
names(BV) = gsub(pattern = ".tailed_reads", replacement = "", x = names(BV))

names(BV)


m<-list()
m

for( i in c("R-451-5p","R-191-5p","R-10b-5p","R-486-5p","R-126-5p","R-99b-5","R-127-3","R-25-3p","R-143-3p","R-26a-1-5p")) {
  BV%>%
    filter(str_detect(BV$miR.name,i))-> TCV
  
  
  
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
ggarrange(plotlist=m, 
          labels = c("A", "B","C","D","E","F","G","H","i","J"))
