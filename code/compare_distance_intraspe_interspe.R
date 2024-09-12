embeddings_pythoncalib <- read.csv("data/embeddings_pythoncalib.csv",sep=";")

rownames(embeddings_pythoncalib) = paste0(embeddings_pythoncalib$id,embeddings_pythoncalib$view)

embeddings_pythoncalib =embeddings_pythoncalib[,c(-2:-1)]

distpheno = as.matrix(dist(embeddings_pythoncalib))
distpheno = melt(distpheno)

distpheno = distpheno[as.numeric(distpheno$Var1)<as.numeric(distpheno$Var2),,F]
distpheno = distpheno[distpheno$value>0,,F]


data_photos = read.csv("./data/data_photos_UV_pythoncalib.csv",sep=";")

distpheno$sp1 = data_photos[match(gsub("D|V","",distpheno$Var1),data_photos$id),]$sp
distpheno$sp2 = data_photos[match(gsub("D|V","",distpheno$Var2),data_photos$id),]$sp

distpheno$sex1 = data_photos[match(gsub("D|V","",distpheno$Var2),data_photos$id),]$sex
distpheno$sex2 = data_photos[match(gsub("D|V","",distpheno$Var2),data_photos$id),]$sex

distpheno$intraspe = as.character(distpheno$sp1)==as.character(distpheno$sp2)
distpheno$view1 = str_sub(distpheno$Var1,-1,-1)
distpheno$view2 = str_sub(distpheno$Var2,-1,-1)
distpheno$D= (distpheno$view1=="D") & (distpheno$view2=="D")
distpheno$V= (distpheno$view1=="V") & (distpheno$view2=="V")
distpheno$Fem= (distpheno$sex1=="F") & (distpheno$sex2=="F")
distpheno$M= (distpheno$sex1=="M") & (distpheno$sex2=="M")

distpheno = distpheno[complete.cases(distpheno),,]
distphenoFD = distpheno[(distpheno$D==T) & (distpheno$Fem==T),,F]
distphenoFV = distpheno[(distpheno$V==T) & (distpheno$Fem==T),,F]
distphenoMD = distpheno[(distpheno$D==T) & (distpheno$M==T),,F]
distphenoMV = distpheno[(distpheno$V==T) & (distpheno$M==T),,F]

library(patchwork)

p1 = ggplot(distphenoFD,aes(x=intraspe,y=value))+
  geom_violin(fill="darkblue")+
  stat_compare_means()+theme_classic()

p2 = ggplot(distphenoFV,aes(x=intraspe,y=value))+
  geom_violin(fill="lightblue")+
  stat_compare_means()+theme_classic()


p3 = ggplot(distphenoMD,aes(x=intraspe,y=value))+
  geom_violin(fill="darkgreen")+
  stat_compare_means()+theme_classic()

p4 = ggplot(distphenoMV,aes(x=intraspe,y=value))+
  geom_violin(fill="lightgreen")+
  stat_compare_means()+theme_classic()

(p1|p2)/(p3|p4)
