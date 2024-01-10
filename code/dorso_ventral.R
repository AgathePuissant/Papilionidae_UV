#------------------ Prepare libraries and functions ----------------------------

library(diverge)
library(reshape2)
library(ggpubr)
library(ape)
library(ggpubr)
library(phytools)
library(RColorBrewer)
library(geomorph)
library(caper)

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



list_get_phenotype = get_phenotype(c("F"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV.csv",path_coords = pc)
meanphen <- list_get_phenotype[[1]]
data_FD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_FD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("F"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV.csv",path_coords = pc)
meanphen <- list_get_phenotype[[1]]
data_FV <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))

list_match <- match_tree(meanphen_match = meanphen, data_match = data_FV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FV <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV.csv",path_coords = pc)
meanphen <- list_get_phenotype[[1]]
data_MD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_MD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_MD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV.csv",path_coords = pc)
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

subtree_analysis = subtree #Save the tree to use in comparisons


pc="./data/pca_embeddings_match.csv"
#------------------ Get the specimens data (visible) -------------------------------------

pc="./data/pca_embeddings_visible.csv"

list_get_phenotype = get_phenotype(c("F"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_visible.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_FD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("F"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_visible.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FV <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))

list_match <- match_tree(meanphen_match = meanphen, data_match = data_FV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FV <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_visible.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_MD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_MD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_MD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_visible.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_MV <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))

list_match <- match_tree(meanphen_match = meanphen, data_match = data_MV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_MV <- list_match[[2]]

meanphen = rbind(meanphen_FD,meanphen_FV,meanphen_MD,meanphen_MV)
data = rbind(data_FD,data_FV,data_MD,data_MV)

meanphen_FD = meanphen_FD[rownames(meanphen_FD_grayscale),]
meanphen_FV = meanphen_FV[rownames(meanphen_FV_grayscale),]
meanphen_MD = meanphen_MD[rownames(meanphen_MD_grayscale),]
meanphen_MV = meanphen_MV[rownames(meanphen_MV_grayscale),]


#------------------ Compare rates of evolution ---------------------------------
sex="F" #Change this parameter to test for males or females

if (sex=="M"){
  meanphen_D_vis = meanphen_MD
  meanphen_V_vis = meanphen_MV
  meanphen_D_grayscale = meanphen_MD_grayscale
  meanphen_V_grayscale = meanphen_MV_grayscale
  meanphen_D_grayscale = meanphen_D_grayscale[rownames(meanphen_V_grayscale),,drop=F]
  meanphen_V_grayscale = meanphen_V_grayscale[rownames(meanphen_D_grayscale),,drop=F]
  meanphen_D_grayscale = meanphen_D_grayscale[complete.cases(meanphen_D_grayscale),,drop=F]
  meanphen_V_grayscale = meanphen_V_grayscale[complete.cases(meanphen_V_grayscale),,drop=F]
}else
{
  meanphen_D_vis = meanphen_FD
  meanphen_V_vis = meanphen_FV
  meanphen_D_grayscale = meanphen_FD_grayscale
  meanphen_V_grayscale = meanphen_FV_grayscale
}


#Rate ratio UV
res.comp_grayscale = compare.multi.evol.rates(as.matrix(cbind(meanphen_D_grayscale,meanphen_V_grayscale)),subtree_analysis,gp=c(rep("D",length(meanphen_D_grayscale)),rep("V",length(meanphen_V_grayscale))))
res.comp_grayscale

#Rate ratio visible
res.comp_vis = compare.multi.evol.rates(as.matrix(cbind(meanphen_D_vis,meanphen_V_vis)),subtree_analysis,gp=c(rep("D",length(meanphen_D_vis)),rep("V",length(meanphen_V_vis))))
res.comp_vis

#------------------ PGLS dorso-ventral -----------------------------------------

#Compute dorso-ventral distance
DV_visible = as.data.frame(diag(as.matrix(dist(meanphen_D_vis, meanphen_V_vis))))
colnames(DV_visible)<-"DV_vis"

DV_grayscale = as.data.frame(diag(as.matrix(dist(meanphen_D_grayscale, meanphen_V_grayscale))))
colnames(DV_grayscale)<-"DV_grayscale"
rownames(DV_grayscale) <- rownames(meanphen_D_grayscale)

#Create dataframe for PGLS
df_cap = data.frame(rownames(DV_visible),DV_visible$DV_vis,DV_grayscale$DV_grayscale)
colnames(df_cap) <- c("Species","DV_visible","DV_grayscale")

#Perform and plot PGLS
comp.data<-comparative.data(subtree_analysis, df_cap, names.col="Species", vcv.dim=2, warn.dropped=TRUE)
model<-pgls(DV_grayscale~DV_visible, data=comp.data, lambda="ML")
summary(model)


ggplot(df_cap, aes(x=DV_visible,y=DV_grayscale))+
  geom_point()+
  geom_abline(aes(intercept = coef(model)[1], slope = coef(model)[2]), linetype="dashed")+
  theme_classic()
