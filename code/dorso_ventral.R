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

pc = "./data/pca_embeddings_UV_match_all.csv"
lvl = "sp"
if (lvl == "form"){
  adp=T
}else{
  adp = F
}


list_get_phenotype = get_phenotype(c("F"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV_pythoncalib.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_FD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("F"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV_pythoncalib.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_FV <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))

list_match <- match_tree(meanphen_match = meanphen, data_match = data_FV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_FV <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV_pythoncalib.csv",path_coords = pc, reduce_dataset=T)
meanphen <- list_get_phenotype[[1]]
data_MD <- list_get_phenotype[[2]]
sp_data <- list_get_phenotype[[4]]
rm(list=c("list_get_phenotype"))


list_match <- match_tree(meanphen_match = meanphen, data_match = data_MD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
subtree <- list_match[[1]]
meanphen_MD <- list_match[[2]]

list_get_phenotype = get_phenotype(c("M"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV_pythoncalib.csv",path_coords = pc, reduce_dataset=T)
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


list_sp <- read.csv("~/GitHub/Papilionidae_UV/data/list_sp.csv", sep=";")
# meanphen_MD_grayscale = meanphen_MD_grayscale[rownames(meanphen_MD_grayscale) %in% list_sp$x,,F]
# meanphen_MV_grayscale = meanphen_MV_grayscale[rownames(meanphen_MV_grayscale) %in% list_sp$x,,F]
# meanphen_FD_grayscale = meanphen_FD_grayscale[rownames(meanphen_FD_grayscale) %in% list_sp$x,,F]
# meanphen_FV_grayscale = meanphen_FV_grayscale[rownames(meanphen_FV_grayscale) %in% list_sp$x,,F]



#------------------ Get the specimens data (visible) -------------------------------------

pc="./data/pca_embeddings_match_all.csv"

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
# 
# meanphen_FD = meanphen_FD[rownames(meanphen_FD_grayscale),]
# meanphen_FV = meanphen_FV[rownames(meanphen_FV_grayscale),]
# meanphen_MD = meanphen_MD[rownames(meanphen_MD_grayscale),]
# meanphen_MV = meanphen_MV[rownames(meanphen_MV_grayscale),]


#------------------ Compare rates of evolution ---------------------------------
sex="M" #Change this parameter to test for males or females

if (sex=="M"){
  meanphen_D_vis = meanphen_MD
  meanphen_V_vis = meanphen_MV
  meanphen_D_grayscale = meanphen_MD_grayscale
  meanphen_V_grayscale = meanphen_MV_grayscale
  meanphen_D_grayscale = meanphen_D_grayscale[rownames(meanphen_V_grayscale),,drop=F]
  meanphen_V_grayscale = meanphen_V_grayscale[rownames(meanphen_D_grayscale),,drop=F]
  meanphen_D_grayscale = meanphen_D_grayscale[complete.cases(meanphen_D_grayscale),,drop=F]
  meanphen_V_grayscale = meanphen_V_grayscale[complete.cases(meanphen_V_grayscale),,drop=F]
  meanphen_D_vis = meanphen_D_vis[rownames(meanphen_V_vis),,drop=F]
  meanphen_V_vis = meanphen_V_vis[rownames(meanphen_D_vis),,drop=F]
  meanphen_D_vis = meanphen_D_vis[complete.cases(meanphen_D_vis),,drop=F]
  meanphen_V_vis = meanphen_V_vis[complete.cases(meanphen_V_vis),,drop=F]
  
  }else
{
  meanphen_D_vis = meanphen_FD
  meanphen_V_vis = meanphen_FV
  meanphen_D_grayscale = meanphen_FD_grayscale
  meanphen_V_grayscale = meanphen_FV_grayscale
  meanphen_D_grayscale = meanphen_D_grayscale[rownames(meanphen_V_grayscale),,drop=F]
  meanphen_V_grayscale = meanphen_V_grayscale[rownames(meanphen_D_grayscale),,drop=F]
  meanphen_D_grayscale = meanphen_D_grayscale[complete.cases(meanphen_D_grayscale),,drop=F]
  meanphen_V_grayscale = meanphen_V_grayscale[complete.cases(meanphen_V_grayscale),,drop=F]
  meanphen_D_vis = meanphen_D_vis[rownames(meanphen_V_vis),,drop=F]
  meanphen_V_vis = meanphen_V_vis[rownames(meanphen_D_vis),,drop=F]
  meanphen_D_vis = meanphen_D_vis[complete.cases(meanphen_D_vis),,drop=F]
  meanphen_V_vis = meanphen_V_vis[complete.cases(meanphen_V_vis),,drop=F]
}

subtree_analysis_1 = match.phylo.data(subtree,meanphen_D_grayscale)$phy
meanphen_D_grayscale = match.phylo.data(subtree,meanphen_D_grayscale)$data
meanphen_V_grayscale = match.phylo.data(subtree,meanphen_V_grayscale)$data
dim(meanphen_V_grayscale)

#Rate ratio UV
res.comp_grayscale = compare.multi.evol.rates(as.matrix(cbind(meanphen_D_grayscale,meanphen_V_grayscale)),subtree_analysis_1,gp=c(rep("D",length(meanphen_D_grayscale)),rep("V",length(meanphen_V_grayscale))))
res.comp_grayscale

subtree_analysis_2 = match.phylo.data(subtree,meanphen_D_vis)$phy
#Rate ratio visible
res.comp_vis = compare.multi.evol.rates(as.matrix(cbind(meanphen_D_vis,meanphen_V_vis)),subtree_analysis_2,gp=c(rep("D",length(meanphen_D_vis)),rep("V",length(meanphen_V_vis))))
res.comp_vis


#---------------Prepare-------------------

meanphen_D_vis2 = meanphen_D_vis
meanphen_D_vis2$spname = rownames(meanphen_D_vis2)
meanphen_V_vis2 = meanphen_V_vis
meanphen_V_vis2$spname = rownames(meanphen_V_vis2)

meanphen_D_grayscale2 = meanphen_D_grayscale
meanphen_D_grayscale2$spname = rownames(meanphen_D_grayscale2)
meanphen_V_grayscale2 = meanphen_V_grayscale
meanphen_V_grayscale2$spname = rownames(meanphen_V_grayscale2)


MF_visible2 = MF_visible
MF_visible2$spname=MF_visible2$Var1
MF_visible2 = MF_visible2[,c(7,12)]
MF_visible2 = MF_visible2[!duplicated(MF_visible2),]
MF_visible3 = MF_visible
MF_visible3$spname=MF_visible3$Var2
MF_visible3 = MF_visible3[,c(7,12)]
MF_visible3 = MF_visible3[!duplicated(MF_visible3),]

meanphen_D_vis2$overlap = NA
# First join with MF_visible2
meanphen_D_vis2 <- meanphen_D_vis2 %>%
  left_join(MF_visible2, by = "spname", suffix = c("", ".MF2")) %>%
  mutate(overlap = coalesce(overlap, overlap.MF2)) 
meanphen_D_vis2$overlap.MF2 = NULL
# Second join with MF_visible3
meanphen_D_vis2 <- meanphen_D_vis2 %>%
  left_join(MF_visible3, by = "spname", suffix = c("", ".MF3")) %>%
  mutate(overlap = coalesce(overlap, overlap.MF3))
meanphen_D_vis2$overlap.MF3 = NULL


meanphen_V_vis2$overlap = NA
# First join with MF_visible2
meanphen_V_vis2 <- meanphen_V_vis2 %>%
  left_join(MF_visible2, by = "spname", suffix = c("", ".MF2")) %>%
  mutate(overlap = coalesce(overlap, overlap.MF2)) 
meanphen_V_vis2$overlap.MF2 = NULL
# Second join with MF_visible3
meanphen_V_vis2 <- meanphen_V_vis2 %>%
  left_join(MF_visible3, by = "spname", suffix = c("", ".MF3")) %>%
  mutate(overlap = coalesce(overlap, overlap.MF3))
meanphen_V_vis2$overlap.MF3 = NULL


meanphen_D_grayscale2$overlap = NA
# First join with MF_visible2
meanphen_D_grayscale2 <- meanphen_D_grayscale2 %>%
  left_join(MF_visible2, by = "spname", suffix = c("", ".MF2")) %>%
  mutate(overlap = coalesce(overlap, overlap.MF2)) 
meanphen_D_grayscale2$overlap.MF2 = NULL
# Second join with MF_visible3
meanphen_D_grayscale2 <- meanphen_D_grayscale2 %>%
  left_join(MF_visible3, by = "spname", suffix = c("", ".MF3")) %>%
  mutate(overlap = coalesce(overlap, overlap.MF3))
meanphen_D_grayscale2$overlap.MF3 = NULL


meanphen_V_grayscale2$overlap = NA
# First join with MF_visible2
meanphen_V_grayscale2 <- meanphen_V_grayscale2 %>%
  left_join(MF_visible2, by = "spname", suffix = c("", ".MF2")) %>%
  mutate(overlap = coalesce(overlap, overlap.MF2)) 
meanphen_V_grayscale2$overlap.MF2 = NULL
# Second join with MF_visible3
meanphen_V_grayscale2 <- meanphen_V_grayscale2 %>%
  left_join(MF_visible3, by = "spname", suffix = c("", ".MF3")) %>%
  mutate(overlap = coalesce(overlap, overlap.MF3))
meanphen_V_grayscale2$overlap.MF3 = NULL

#---------------Overlap select-------------------

overlap_bin = meanphen_D_vis2[meanphen_D_vis2$overlap<=0,]$spname

meanphen_D_vis2 = meanphen_D_vis2[meanphen_D_vis2$spname %in% overlap_bin,]
rownames(meanphen_D_vis2)=meanphen_D_vis2$spname
meanphen_D_vis2$spname=NULL
meanphen_D_vis2$overlap=NULL


meanphen_V_vis2 = meanphen_V_vis2[meanphen_V_vis2$spname %in% overlap_bin,]
rownames(meanphen_V_vis2)=meanphen_V_vis2$spname
meanphen_V_vis2$spname=NULL
meanphen_V_vis2$overlap=NULL

meanphen_D_grayscale2 = meanphen_D_grayscale2[meanphen_D_grayscale2$spname %in% overlap_bin,]
rownames(meanphen_D_grayscale2)=meanphen_D_grayscale2$spname
meanphen_D_grayscale2$spname=NULL
meanphen_D_grayscale2$overlap=NULL

meanphen_V_grayscale2 = meanphen_V_grayscale2[meanphen_V_grayscale2$spname %in% overlap_bin,]
rownames(meanphen_V_grayscale2)=meanphen_V_grayscale2$spname
meanphen_V_grayscale2$spname=NULL
meanphen_V_grayscale2$overlap=NULL

#---------------rate ratio symp-------------------

subtree_analysis_1 = match.phylo.data(subtree,meanphen_D_grayscale2)$phy
meanphen_D_grayscale2 = match.phylo.data(subtree,meanphen_D_grayscale2)$data
meanphen_V_grayscale2 = match.phylo.data(subtree,meanphen_V_grayscale2)$data
dim(meanphen_V_grayscale2)

#Rate ratio UV
res.comp_grayscale = compare.multi.evol.rates(as.matrix(cbind(meanphen_D_grayscale2,meanphen_V_grayscale2)),subtree_analysis_1,gp=c(rep("D",length(meanphen_D_grayscale2)),rep("V",length(meanphen_V_grayscale2))))
res.comp_grayscale

subtree_analysis_2 = match.phylo.data(subtree,meanphen_D_vis2)$phy
#Rate ratio visible
res.comp_vis = compare.multi.evol.rates(as.matrix(cbind(meanphen_D_vis2,meanphen_V_vis2)),subtree_analysis_2,gp=c(rep("D",length(meanphen_D_vis2)),rep("V",length(meanphen_V_vis2))))
res.comp_vis


# compare.multi.evol.rates(as.matrix(cbind(meanphen_V_vis,meanphen_V_grayscale)),subtree_analysis_2,gp=c(rep("V_vis",length(meanphen_V_vis)),rep("V_grayscale",length(meanphen_V_grayscale))),Subset=F)

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


