###network
rm(list = ls())
library(tidyverse)
library(RColorBrewer)
library(scales)
library(igraph)

setwd("C:/cellphoneDB")
pvalues<-read.table("means.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
pvalues=pvalues[,12:ncol(pvalues)]
statdf=as.data.frame(colSums(pvalues))
colnames(statdf)=c("number")

statdf$indexb=str_replace(rownames(statdf),"^.*[|]","")
statdf$indexa=str_replace(rownames(statdf),"[|]..*$","")
rankname=sort(unique(statdf$indexa)) 

A=c()
B=c()
C=c()
remaining=rankname
for (i in rankname[-6]) {
  remaining=setdiff(remaining,i)
  for (j in remaining) {
    count=statdf[statdf$indexa == i & statdf$indexb == j,"number"]+
      statdf[statdf$indexb == i & statdf$indexa == j,"number"]
    A=append(A,i)
    B=append(B,j)
    C=append(C,count)
  }
}

statdf2=data.frame(indexa=A,indexb=B,number=C)
statdf2=statdf2 %>% rbind(statdf[statdf$indexa==statdf$indexb,c("indexa","indexb","number")])
statdf2=statdf2[statdf2$number > 0,] #过滤掉值为0的观测

#Set the color of nodes and connections
color1=c("#8DD3C7", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")
names(color1)=rankname
color2=colorRampPalette(brewer.pal(9, "Reds")[3:7])(20) 
names(color2)=1:20 


net <- graph_from_data_frame(statdf2[,c("indexa","indexb","number")])
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = order(membership(group)))

E(net)$width <- E(net)$number / 10 
E(net)$color <- color2[as.character(ifelse(E(net)$number > 20,20,E(net)$number))] 
E(net)$label = E(net)$number
E(net)$label.color <- "black" 
V(net)$label.color <- "black" 
V(net)$color <- color1[names(V(net))] 


loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

#pdf("interaction.num.3.pdf",width = 6,height = 6)
plot(net,
     edge.arrow.size = 0,
     edge.curved = 0, 
     vertex.frame.color = "black", 
     layout = coords,
     vertex.label.cex = 1,
     vertex.size = 30) 
#dev.off()

###gene/cell model
a<-read.table("GSE148476.txt",header = T,sep = "\t")
a<-a[,-1]
info<-data.frame(sample=colnames(a),RN=t(a[1,]))
info<-info[-1,]
colnames(info)[2]<-"RN"
gene<-read.table("F:/gene.txt",header = F,sep = "\t")
a1<-a[a$Symbol%in%gene$V1,]

b<-data.frame(lapply(a1[,2:ncol(a1)],as.numeric))
b<-cbind(Sample=a1$Symbol,b)
#result<-aggregate(b,by = list(b$Sample),FUN = mean)
#b<-na.omit(b)
value_1<-info[which(info$RN=="Responder"),]
value_1<-value_1$sample
value_2<-info[which(info$RN=="nonResponder"),]
value_2<-value_2$sample
n1 <- sample(c(1:2),30,replace = T,prob = c(0.5, 0.5))
sample1 <- value_1[n1==1]
data1<-b[,sample1]
rownames(data1)<-b[,1]
sample2 <- value_1[n1==2]
data2<-b[,sample2]
rownames(data2)<-b[,1]

n2 <- sample(c(1:2),32,replace = TRUE,prob = c(0.5, 0.5))
sample3 <- value_2[n2==1]
data3<-b[,sample3]
rownames(data1)<-b[,1]
sample4 <- value_2[n2==2]
data4<-b[,sample4]
rownames(data2)<-b[,1]


data11<-cbind(data1,data3)
tra<-t(data11)
tra<-cbind(sample=rownames(tra),tra)
tra<-merge(tra,info,by="sample")

data22<-cbind(data2,data4)
val<-t(data22)
val<-cbind(sample=rownames(val),val)
val<-merge(val,info,by="sample")

a<-read.table("train-gene.txt",header = T,sep = "\t")
b<-a[,-1]
b<-data.frame(lapply(b[,1:37],as.numeric))
b<-cbind(RN=a$RN,b)
b$RN<-factor(b$RN,levels = c("nonResponder","Responder"))

#write.table(b,"exp.txt",quote = F,sep = "\t",row.names = F)
x <- as.matrix(b[,2:ncol(b)])
d1<-c()
for (i in 1:dim(x)[2]) {
  fit.ful<-glm(RN ~ x[,i], family = binomial(), 
               data = b,control=list(maxit=100))
  c<-summary(fit.ful)
  d<-c$coefficients
  d<-d[2,]
  d1<-rbind(d1,d)
}
rownames(d1)<-colnames(x)
d1<-data.frame(d1)
d1<-d1[which(d1$Pr...z..<0.05),]

value<-row.names(d1)
b1<-b[,value]
b1<-cbind(RN=b$RN,b1)
x <- as.matrix(b1[,2:ncol(b1)])
uni.ful<-glm(RN ~., family = binomial(),data = b1,control=list(maxit=100))
c<-summary(uni.ful)
d<-c$coefficients
val<-read.table("validation-gene.txt",header = T,sep = "\t")
val<-val[,-1]
val_1<-data.frame(lapply(val[,1:37],as.numeric))
val_1<-cbind(RN=val$RN,val_1)
val_1$RN<-factor(val_1$RN,levels = c("nonResponder","Responder"))
val_2<-val_1[,value]
val_2<-cbind(RN=val_1$RN,val_2)
val_2$RN<-factor(val_2$RN,levels = c("nonResponder","Responder"))

library(pROC)
#pr1 <- predict(uni.ful,newdata= val_2,type = c("response"))
m<-data.frame(coef=coef(uni.ful))
m_1<-as.matrix(m[-1,])
rownames(m_1)<-rownames(m)[2:nrow(m)] 
k<-colnames(val_2)[2:ncol(val_2)]
m_1<-m_1[k,]
val_3<-as.matrix(val_2[,2:ncol(val_2)])
y=as.vector(val_3%*%m_1)

#Training set
#pr1<-predict(uni.ful,newdata= val_2,type=c("response"))
pr2<-predict(uni.ful,data=b1,type=c("response"))
##ROC
roccurve <- roc(b1$RN ~ pr2)
value_1<-auc(roccurve)

roccurve1 <- roc(val_2$RN~y)
value_2<-auc(roccurve1)

value_1
value_2
d

plot.roc(roccurve,xlim = c(1,0),ylim=c(0,1),
         identity.lwd=2,identity.col="black",identity.lty=2,
         print.auc=TRUE,col="red")
plot.roc(roccurve1,xlim = c(1,0),ylim=c(0,1),
         identity.lwd=2,identity.col="black",identity.lty=2,
         print.auc=TRUE,col="red")

#modelFit <- train(RN~.,data=b1,method="glm")
#predictions <- predict(modelFit,newdata=b1)
#confusionMatrix(b1$RN,predictions)


completepath <- list.files(reviewpath, pattern = "*.txt$", full.names = TRUE)
review <- lapply(completepath, read.table)
docname <- list.files(reviewpath, pattern = "*.txt$")
docname<-data.frame(docname)
library(tidyverse)
docname<-separate(docname,docname,c("name"),sep = "[.]")
names(review)<-as.vector(docname$name)

library(do)
total<-c()
for (i in 1:length(review)) {
  val<-data.frame(review[[i]])
  colnames(val)<-val[1,]
  val<-val[-1,]
  val$Symbol<-Replace(val$Symbol,from="-",to=".")
  
  #if(sum(val$Symbol%in%value)==length(value)){
  val<-val[,-1]
  a_1<-val[,2:ncol(val)]
  a_1<-a_1[2:nrow(val),]
  b<-data.frame(lapply(a_1,as.numeric))
  b<-scale(b)
  a_2<-data.frame(cbind(Symbol=val$Symbol[2:nrow(val)],b))
  a_2<-rbind(val[1,],a_2)
  
  a_3<-a_2[a_2$Symbol%in%value,]
  a_4<-rbind(a_2[1,],a_3) 
  a_4<-a_4[-1,]
  a_5<-data.frame(lapply(a_4[,2:ncol(a_4)], as.numeric))
  rownames(a_5)<-a_4$Symbol
  anno<-data.frame(a_2[1,])
  anno<-anno[,-1]
  anno<-data.frame(t(anno))
  colnames(anno)<-"RN"
  anno_colors =list(RN=c(Responder="#364C8E", nonResponder="#CB3A1F"))
  library(ggplot2)
  pheatmap::pheatmap(a_5,show_colnames = F,
                     cluster_rows=F,cluster_cols = F,
                   scale = "row",angle_col = 0,
                     annotation_colors =anno_colors,
                     annotation_col =anno,
                     color = colorRampPalette(colors = c("#2F2D58","white","#E03322"))(100),
                     border_color = "grey",
                     breaks=c(seq(-3,-0.0000001,length.out=50),0,
                              seq(0.0000001,3,length.out=50)))
  
  
  val_1<-val[val$Symbol%in%value,]
  val_1<-rbind(val[1,],val_1) 
  val_1<-t(val_1)
  colnames(val_1)<-val_1[1,]
  val_1<-data.frame(val_1[-1,])
  val_2<-data.frame(lapply(val_1[,],as.numeric))
  val_2<-cbind(RN=val_1$NA.,val_2)
  val_2<-val_2[,-2]
  val_2$RN<-factor(val_2$RN,levels = c("nonResponder","Responder"))
  
  library(pROC)
  library(caret)
  # pr1 <- predict(uni.ful,newdata= val_2,type = c("response"))
  
  
  m<-data.frame(coef=coef(uni.ful))
  m_1<-as.matrix(m[-1,])
  rownames(m_1)<-rownames(m)[2:nrow(m)] 
  
  k<-colnames(val_2)[2:ncol(val_2)]
  m_1<-m_1[k,]
  val_3<-as.matrix(val_2[,2:ncol(val_2)])
  y=as.vector(val_3%*%m_1)
  
  ##Validation set
  #pr2<-predict(uni.ful,data=b1,type=c("response"))
  ##做ROC
  #roccurve <- roc(b1$RN ~ pr2)
  #value_1<-auc(roccurve)
  
  roccurve1 <- roc(val_2$RN~y)
  value_2<-auc(roccurve1)
  #}
  #else{
  #value_2<-NA
  
  total<-c(total,value_2)
  plot.roc(roccurve1,xlim = c(1,0),ylim=c(0,1),
           identity.lwd=2,identity.col="black",identity.lty=2,
           print.auc=TRUE,col="red")
}
result<-data.frame(GSE=docname$name,roc=total)
#write.table(result,"result-gene.txt",quote = F,sep = "\t",row.names = F)

plot.roc(roccurve,xlim = c(1,0),ylim=c(0,1))
plot.roc(roccurve1,xlim = c(1,0),ylim=c(0,1))


###heatmap
pheatmap::pheatmap(martix,show_colnames = T,
                   cluster_rows=F,cluster_cols = F,
                   border_color = "grey",
                   color = colorRampPalette(colors = c("white","#E03322"))(100),
                   breaks=c(seq(0,8.99999999,length.out=50),9,
                            seq(9.0000001,50,length.out=50)))

