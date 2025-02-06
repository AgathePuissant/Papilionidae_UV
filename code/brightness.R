#------------------ Prepare libraries and functions ----------------------------

library(diverge)
library(reshape2)
library(ggpubr)
library(ape)
library(ggpubr)
library(phytools)
library(RColorBrewer)

source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_phenotype.R")

#------------------ Get the specimens data (UV) -------------------------------------

pc = "./data/pca_embeddings_UV.csv"
lvl = "sp"
if (lvl == "form"){
  adp=T
}else{
  adp = F
}


list_get_phenotype = get_phenotype(c("F"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_FD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("F"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FV <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))

list_match <- match_tree(meanphen_match = meanphen, data_match = data_FV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FV <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_MD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_MD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_MD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_MV <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))

list_match <- match_tree(meanphen_match = meanphen, data_match = data_MV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_MV <- list_match[[2]]

meanphen_grayscale = rbind(meanphen_FD,meanphen_FV,meanphen_MD,meanphen_MV)
data_grayscale = rbind(data_FD,data_FV,data_MD,data_MV)


meanphen_FD_grayscale = meanphen_FD
meanphen_FV_grayscale = meanphen_FV
meanphen_MD_grayscale = meanphen_MD
meanphen_MV_grayscale = meanphen_MV



#------------- Merge the data with mean brightness info ------------------------

meanb <- read.csv("./data/mean_brightness.csv", sep=";")


data_grayscale <- merge(data_grayscale,meanb, by = c("id","view"))


#------------- Prepare data for pairwise comparisons ---------------------------

my_comparisons <- expand.grid(unique(paste0(data_grayscale$sex,data_grayscale$view)),unique(paste0(data_grayscale$sex,data_grayscale$view)))
my_comparisons <- my_comparisons[as.numeric(my_comparisons[,1])>as.numeric(my_comparisons[,2]),]
my_comparisons <- split(as.matrix(my_comparisons), 1:nrow(my_comparisons))
my_comparisons <- lapply(my_comparisons,function(x) {return(as.vector(x))})
my_comparisons <- my_comparisons[c(1:2,5:6)]


data_grayscale$group = paste0(data_grayscale$sex,data_grayscale$view)


stat.test <- compare_means(
  meanb ~ group , data = data_grayscale,
  method = "t.test"
)

stat.test <- stat.test[!((stat.test$group1=="MD" & stat.test$group2=="FV") | (stat.test$group1=="FV" & stat.test$group2=="MD")),]
stat.test <- stat.test[!((stat.test$group1=="MV" & stat.test$group2=="FD") | (stat.test$group1=="FD" & stat.test$group2=="MV")),]

#------------- Compare Male Dorsal and Male Ventral  ---------------------------

x1= data_grayscale[data_grayscale$view=="D" & data_grayscale$sex=="M",,drop=F]$meanb
names(x1) <- data_grayscale[data_grayscale$view=="D" & data_grayscale$sex=="M",,drop=F]$tipsgenre

x2= data_grayscale[data_grayscale$view=="V" & data_grayscale$sex=="M",,drop=F]$meanb
names(x2) <- data_grayscale[data_grayscale$view=="V" & data_grayscale$sex=="M",,drop=F]$tipsgenre

matching = match.phylo.data(subtree,x1)
subtree_analysis = matching$phy
# x1 = matching$data
# x2 = x2[names(x1)]

x2 = x2[names(x2) %in% subtree_analysis$tip.label]
x1 = x1[names(x1) %in% subtree_analysis$tip.label]

phyl.pairedttest(subtree_analysis, x1,x2)
pval = phyl.pairedttest(subtree_analysis, x1,x2)$P.dbar*4

d <- data.frame(tipsgenre = c(names(x1),names(x2)), meanb = c(x1,x2), view = c(rep("x1",length(x1)),rep("x2",length(x2))))

ggplot(d, aes(view, meanb)) + 
  geom_violin() +
  geom_point() +
  geom_line(aes(group = tipsgenre)) +
  theme_classic()+
  annotate("text",1,0.4, label=as.character(round(pval,3)))

stat.test[(stat.test$group1=="MD" & stat.test$group2=="MV") | (stat.test$group1=="MV" & stat.test$group2=="MD"),]$p.adj=pval


#------------- Compare Female Dorsal and Female Ventral  -----------------------

x1= data_grayscale[data_grayscale$view=="D" & data_grayscale$sex=="F",,drop=F]$meanb
names(x1) <- data_grayscale[data_grayscale$view=="D" & data_grayscale$sex=="F",,drop=F]$tipsgenre

x2= data_grayscale[data_grayscale$view=="V" & data_grayscale$sex=="F",,drop=F]$meanb
names(x2) <- data_grayscale[data_grayscale$view=="V" & data_grayscale$sex=="F",,drop=F]$tipsgenre

matching = match.phylo.data(subtree,x2)
subtree_analysis = matching$phy
# x2 = matching$data
# x1 = x1[names(x2)]

x2 = x2[names(x2) %in% subtree_analysis$tip.label]
x1 = x1[names(x1) %in% subtree_analysis$tip.label]

phyl.pairedttest(subtree_analysis, x1,x2)
pval = phyl.pairedttest(subtree_analysis, x1,x2)$P.dbar*4

d <- data.frame(tipsgenre = c(names(x1),names(x2)), meanb = c(x1,x2), view = c(rep("x1",length(x1)),rep("x2",length(x2))))

ggplot(d, aes(view, meanb)) + 
  geom_violin() +
  geom_point() +
  geom_line(aes(group = tipsgenre)) +
  theme_classic()+
  annotate("text",1,0.4, label=as.character(round(pval,3)))

stat.test[(stat.test$group1=="FD" & stat.test$group2=="FV") | (stat.test$group1=="FV" & stat.test$group2=="FD"),]$p.adj=pval

#------------- Compare Female Dorsal and Male Dorsal  --------------------------

x1= data_grayscale[data_grayscale$view=="D" & data_grayscale$sex=="F",,drop=F]$meanb
names(x1) <- data_grayscale[data_grayscale$view=="D" & data_grayscale$sex=="F",,drop=F]$tipsgenre

x2= data_grayscale[data_grayscale$view=="D" & data_grayscale$sex=="M",,drop=F]$meanb
names(x2) <- data_grayscale[data_grayscale$view=="D" & data_grayscale$sex=="M",,drop=F]$tipsgenre

matching = match.phylo.data(subtree,x2)
subtree_analysis = matching$phy
# x2 = matching$data
# x1 = x1[names(x2)]
# x1=x1[!is.na(x1)]
# x2 = x2[names(x1)]
# matching = match.phylo.data(subtree,x2)
# subtree_analysis = matching$phy

x2 = x2[names(x2) %in% subtree_analysis$tip.label]
x1 = x1[names(x1) %in% subtree_analysis$tip.label]

phyl.pairedttest(subtree_analysis, x1,x2)
pval = phyl.pairedttest(subtree_analysis, x1,x2)$P.dbar*4

d <- data.frame(tipsgenre = c(names(x1),names(x2)), meanb = c(x1,x2), view = c(rep("x1",length(x1)),rep("x2",length(x2))))

ggplot(d, aes(view, meanb)) + 
  geom_violin() +
  geom_point() +
  geom_line(aes(group = tipsgenre)) +
  theme_classic()+
  annotate("text",1,0.4, label=as.character(round(pval,3)))

stat.test[(stat.test$group1=="MD" & stat.test$group2=="FD") | (stat.test$group1=="FD" & stat.test$group2=="MD"),]$p.adj=pval

#------------- Compare Female Ventral and Male Ventral  ------------------------

x1= data_grayscale[data_grayscale$view=="V" & data_grayscale$sex=="F",,drop=F]$meanb
names(x1) <- data_grayscale[data_grayscale$view=="V" & data_grayscale$sex=="F",,drop=F]$tipsgenre

x2= data_grayscale[data_grayscale$view=="V" & data_grayscale$sex=="M",,drop=F]$meanb
names(x2) <- data_grayscale[data_grayscale$view=="V" & data_grayscale$sex=="M",,drop=F]$tipsgenre

matching = match.phylo.data(subtree,x2)
subtree_analysis = matching$phy
# x2 = matching$data
# x1 = x1[names(x2)]
# x1=x1[!is.na(x1)]
# x2 = x2[names(x1)]
# matching = match.phylo.data(subtree,x2)
# subtree_analysis = matching$phy

x2 = x2[names(x2) %in% subtree_analysis$tip.label]
x1 = x1[names(x1) %in% subtree_analysis$tip.label]

phyl.pairedttest(subtree_analysis, x1,x2)
pval = phyl.pairedttest(subtree_analysis, x1,x2)$P.dbar*4

d <- data.frame(tipsgenre = c(names(x1),names(x2)), meanb = c(x1,x2), view = c(rep("x1",length(x1)),rep("x2",length(x2))))

ggplot(d, aes(view, meanb)) + 
  geom_violin() +
  geom_point() +
  geom_line(aes(group = tipsgenre)) +
  theme_classic()+
  annotate("text",1,0.4, label=as.character(round(pval,3)))

stat.test[(stat.test$group1=="FV" & stat.test$group2=="MV") | (stat.test$group1=="MV" & stat.test$group2=="FV"),]$p.adj=pval

#------------- Create summary plot of all comparisons  -------------------------

stat.test$y.position = c(80/255,70/255,70/255,90/255)

stat.test$p.signif <- "n.s."
stat.test$p.signif[stat.test[,5]<0.1] <- "."
stat.test$p.signif[stat.test[,5]<0.05] <- "*"
stat.test$p.signif[stat.test[,5]<0.01] <- "**"
stat.test$p.signif[stat.test[,5]<0.001] <- "***"


data_grayscale$group <- factor(data_grayscale$group, levels = c("FD", "MD", "FV","MV"))


library(RColorBrewer)
color_paired = brewer.pal(n=4,"Paired")
color_paired=color_paired[c(1,3,2,4)]

ggplot(data=data_grayscale, aes(x=group,y=meanb))+
  geom_violin(aes(fill = group))+
  stat_summary(fun.data=mean_sdl, geom="pointrange", alpha =0.3)+
  stat_pvalue_manual(stat.test, label="p.signif",y.position = "y.position")+
  theme_classic()+
  xlab("Sex and wing side")+
  ylab("Mean UV brightness of wings")+
  scale_fill_manual(values=color_paired)+ 
  theme(legend.position = "none",
        text = element_text(size = 15))

#-------------------------------------------------------------------------------




