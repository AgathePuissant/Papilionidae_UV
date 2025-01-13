
#------------------ Prepare libraries and functions ----------------------------

library(reshape2)
library(ape)
library(phytools)
library(caper)
library(mvMORPH)

source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_phenotype.R")

#------------------ Get the specimens data (UV) -------------------------------------

pc = "./data/pca_embeddings_pythoncalib_retrained.csv"
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

subtree_analysis = subtree #Save the tree to use in comparisons

#-----------------------------------------------------------

subtree = match.phylo.data(subtree, meanphen_MV_grayscale)$phy
bm = mvgls(as.matrix(meanphen_MV_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_MV_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_MV_grayscale)~1, tree=subtree, model="EB")
lam_MV = mvgls(as.matrix(meanphen_MV_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_MV)

lam_MV


bm = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="EB")
lam_MD = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_MD)

lam_MD


#-----------------------------------


subtree = match.phylo.data(subtree, meanphen_FV_grayscale)$phy
bm = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="EB")
lam_FV = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_FV)

lam_FV


bm = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="EB")
lam_FD = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_FD)

lam_FD



lam_MD$param ; lam_MV$param
lam_FD$param ; lam_FV$param

