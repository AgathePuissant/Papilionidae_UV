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
library(dplyr)

source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_phenotype.R")

#------------------ Get the specimens data (UV) -------------------------------------

pc = "./data/pca_embeddings_UV_match_all.csv"
lvl = "sp"
if (lvl == "form"){
  adp=T
}else{
  adp = F
}

path_data_UV = "./data/data_photos_UV_pythoncalib.csv"

list_get_phenotype = get_phenotype(c("F"),c("D"), mode = 'mean', level = lvl, path_data_photos = path_data_UV,path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_FD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("F"),c("V"), mode = 'mean', level = lvl, path_data_photos = path_data_UV,path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FV <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))

list_match <- match_tree(meanphen_match = meanphen, data_match = data_FV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FV <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("D"), mode = 'mean', level = lvl, path_data_photos = path_data_UV,path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_MD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_MD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_MD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("V"), mode = 'mean', level = lvl, path_data_photos = path_data_UV,path_coords = pc, reduce_dataset=T)
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

# list_sp <- read.csv("~/GitHub/Papilionidae_UV/data/list_sp.csv", sep=";")
# meanphen_MD_grayscale = meanphen_MD_grayscale[rownames(meanphen_MD_grayscale) %in% list_sp$x,,F]
# meanphen_MV_grayscale = meanphen_MV_grayscale[rownames(meanphen_MV_grayscale) %in% list_sp$x,,F]
# meanphen_FD_grayscale = meanphen_FD_grayscale[rownames(meanphen_FD_grayscale) %in% list_sp$x,,F]
# meanphen_FV_grayscale = meanphen_FV_grayscale[rownames(meanphen_FV_grayscale) %in% list_sp$x,,F]

subtree = match.phylo.data(subtree,meanphen_MD_grayscale)$phy

#------------------ Compute pairwise distance in UV ----------------------------

# sis=extract_sisters(subtree)
sis=read.table("./data/sis_list.csv", head=T,sep=";",row.names = 1)
colnames(sis) <- c("sp1","sp2")

MD = create_distpheno(meanphen_MD,"M","D", level=lvl)
MD$sex="M"
MD$view="D"
MV = create_distpheno(meanphen_MV,"M","V", level=lvl)
MV$sex="M"
MV$view="V"
FD = create_distpheno(meanphen_FD,"F","D", level=lvl)
FD$sex="F"
FD$view="D"
FV = create_distpheno(meanphen_FV,"F","V", level=lvl)
FV$sex="F"
FV$view="V"

M = rbind(MD,MV)
Fem = rbind(FD,FV)

MF_grayscale=rbind(M,Fem)
MF_grayscale$type = "grayscale"

#------------------ Get the specimens data (visible) -------------------------------------

pc="./data/pca_embeddings_match_all.csv"

path_data_visible = "./data/data_photos_visible.csv"

list_get_phenotype = get_phenotype(c("F"),c("D"), mode = 'mean', level = lvl, path_data_photos = path_data_visible,path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_FD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("F"),c("V"), mode = 'mean', level = lvl, path_data_photos = path_data_visible,path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FV <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))

list_match <- match_tree(meanphen_match = meanphen, data_match = data_FV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FV <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("D"), mode = 'mean', level = lvl, path_data_photos = path_data_visible,path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_MD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_MD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_MD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("V"), mode = 'mean', level = lvl, path_data_photos = path_data_visible,path_coords = pc, reduce_dataset=T)
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

# meanphen_FD = meanphen_FD[rownames(meanphen_FD) %in% list_sp$x,,F]
# meanphen_FV = meanphen_FV[rownames(meanphen_FV) %in% list_sp$x,,F]
# meanphen_MD = meanphen_MD[rownames(meanphen_MD) %in% list_sp$x,,F]
# meanphen_MV = meanphen_MV[rownames(meanphen_MV) %in% list_sp$x,,F]



#------------------ Compute pairwise distance in visible -----------------------

MD = create_distpheno(meanphen_MD,"M","D", level=lvl)
MD$sex="M"
MD$view="D"
MV = create_distpheno(meanphen_MV,"M","V", level=lvl)
MV$sex="M"
MV$view="V"
FD = create_distpheno(meanphen_FD,"F","D", level=lvl)
FD$sex="F"
FD$view="D"
FV = create_distpheno(meanphen_FV,"F","V", level=lvl)
FV$sex="F"
FV$view="V"

M = rbind(MD,MV)
Fem = rbind(FD,FV)

MF_visible=rbind(M,Fem)
MF_visible$type = "visible"

#------------------- GLM -------------------------------------------------------

#Data preparation
MF_visible = MF_visible[!is.na(MF_visible$Var1),]
MF_grayscale = MF_grayscale[!is.na(MF_grayscale$Var1),]
MF_visible = MF_visible[match(paste(MF_grayscale$Var1,MF_grayscale$Var2,MF_grayscale$sex,MF_grayscale$view),paste(MF_visible$Var1,MF_visible$Var2,MF_visible$sex,MF_visible$view)),]
MF_grayscale$visdist = MF_visible[match(paste(MF_grayscale$Var1,MF_grayscale$Var2,MF_grayscale$sex,MF_grayscale$view),paste(MF_visible$Var1,MF_visible$Var2,MF_visible$sex,MF_visible$view)),]$value
MF_grayscale = MF_grayscale[MF_grayscale$Var1!=MF_grayscale$Var2,]

list_new=unique(paste0(MF_grayscale$Var1,"-",MF_grayscale$Var2))
setdiff(rownames(sis),list_new)

#GLM
model<-glm(value~(visdist+overlap+distphylo):sex:view, data=MF_grayscale)

summary(model)

#-------------------------------------------------------------------------------