
#------------------ Prepare libraries and functions ----------------------------

library(reshape2)
library(ape)
library(phytools)
library(caper)
library(mvMORPH)
library(motmot)

source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_phenotype.R")

#------------------ Get the specimens data (UV) -------------------------------------

pc = "./data/pca_embeddings_UV.csv"
# pc = "./data/pca_embeddings_retrained.csv"
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

#---------------------------- Rate ratio Males -------------------------------

subtree = match.phylo.data(subtree, meanphen_MV_grayscale)$phy
bm = mvgls(as.matrix(meanphen_MV_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_MV_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_MV_grayscale)~1, tree=subtree, model="EB")
lam_MV = mvgls(as.matrix(meanphen_MV_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_MV)

lam_MV

subtree_lam1 = transformPhylo(subtree, model="lambda", lambda = lam_MV$param)

resbm = mvBM(subtree_lam1, as.matrix(meanphen_MV_grayscale), method="pic")


bm = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="EB")
lam_MD = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_MD)

lam_MD

subtree_lam2 = transformPhylo(subtree, model="lambda", lambda = lam_MD$param)

resbm2 = mvBM(subtree_lam2, as.matrix(meanphen_MD_grayscale), method="pic")




meanphen_Mtot = cbind(meanphen_MD_grayscale,meanphen_MV_grayscale)


bm = mvgls(as.matrix(meanphen_Mtot)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_Mtot)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_Mtot)~1, tree=subtree, model="EB")
lam_Mtot = mvgls(as.matrix(meanphen_Mtot)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_Mtot)

lam_Mtot


subtree_lam3 = transformPhylo(subtree, model="lambda", lambda = lam_Mtot$param)

resbm3 = mvBM(subtree_lam3, as.matrix(meanphen_Mtot), method="pic", param = list(constraint="equal"))

simu_tot = mvSIM(subtree_lam3,nsim=100,model="BM",param=resbm3)

null_ratio = rep(NA, 100)
for (i in c(1:100)){
  print(i)
  dorsal = simu_tot[[i]][,1:11]
  ventral = simu_tot[[i]][,12:22]
  resbm_dorsal = mvBM(subtree_lam2, dorsal, method="pic", echo = F, diagnostic = F)
  resbm_ventral = mvBM(subtree_lam1, ventral, method="pic", echo = F, diagnostic = F)
  null_ratio[i] = mean(diag(resbm_ventral$sigma))/mean(diag(resbm_dorsal$sigma))
  null_ratio[i] = mean(sqrt(diag(resbm_ventral$sigma)))/mean(sqrt(diag(resbm_dorsal$sigma)))
}

real_ratio = mean(diag(resbm$sigma))/mean(diag(resbm2$sigma))
real_ratio
p=sum((null_ratio>real_ratio) & (!is.na(null_ratio)))/sum(!is.na(null_ratio))
p

#--------------- Rate ratio Females -------------


subtree = match.phylo.data(subtree, meanphen_FV_grayscale)$phy
bm = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="EB")
lam_FV = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_FV)

lam_FV

subtree_lam1 = transformPhylo(subtree, model="lambda", lambda = lam_FV$param)

resbm = mvBM(subtree_lam1, as.matrix(meanphen_FV_grayscale), method="pic")


bm = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="EB")
lam_FD = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_FD)

lam_FD

subtree_lam2 = transformPhylo(subtree, model="lambda", lambda = lam_FD$param)

resbm2 = mvBM(subtree_lam2, as.matrix(meanphen_FD_grayscale), method="pic")

meanphen_Ftot = cbind(meanphen_FD_grayscale,meanphen_FV_grayscale)


bm = mvgls(as.matrix(meanphen_Ftot)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_Ftot)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_Ftot)~1, tree=subtree, model="EB")
lam_Ftot = mvgls(as.matrix(meanphen_Ftot)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_Ftot)

lam_Ftot



library(motmot)
subtree_lam3 = transformPhylo(subtree, model="lambda", lambda = lam_Ftot$param)

resbm3 = mvBM(subtree, as.matrix(meanphen_Ftot), method="pic", param = list(constraint=constraints))
resbm3

resbm3 = mvBM(subtree_lam3, as.matrix(meanphen_Ftot), method="pic", param = list(constraint="equal"))

simu_tot = mvSIM(subtree_lam3,nsim=100,model="BM",param=resbm3)

null_ratio = rep(NA, 100)
for (i in c(1:100)){
  print(i)
  dorsal = simu_tot[[i]][,1:11]
  ventral = simu_tot[[i]][,12:22]
  resbm_dorsal = mvBM(subtree_lam2, dorsal, method="pic", echo = F, diagnostic = F)
  resbm_ventral = mvBM(subtree_lam1, ventral, method="pic", echo = F, diagnostic = F)
  null_ratio[i] = mean(diag(resbm_ventral$sigma))/mean(diag(resbm_dorsal$sigma))
  null_ratio[i] = mean(sqrt(diag(resbm_ventral$sigma)))/mean(sqrt(diag(resbm_dorsal$sigma)))
}

real_ratio = mean(diag(resbm$sigma))/mean(diag(resbm2$sigma))
real_ratio
p=sum((null_ratio>real_ratio) & (!is.na(null_ratio)))/sum(!is.na(null_ratio))
p


#------------- Lambda values -------------------------


lam_MD$param ; lam_MV$param
lam_FD$param ; lam_FV$param


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






#------------------ Rate ratio Males -----------------------------

subtree = match.phylo.data(subtree, meanphen_MV)$phy
bm = mvgls(as.matrix(meanphen_MV)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_MV)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_MV)~1, tree=subtree, model="EB")
lam_MV = mvgls(as.matrix(meanphen_MV)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_MV)

lam_MV

subtree_lam1 = transformPhylo(subtree, model="lambda", lambda = lam_MV$param)

resbm = mvBM(subtree_lam1, as.matrix(meanphen_MV), method="pic")


bm = mvgls(as.matrix(meanphen_MD)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_MD)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_MD)~1, tree=subtree, model="EB")
lam_MD = mvgls(as.matrix(meanphen_MD)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_MD)

lam_MD

library(motmot)
subtree_lam2 = transformPhylo(subtree, model="lambda", lambda = lam_MD$param)

resbm2 = mvBM(subtree_lam2, as.matrix(meanphen_MD), method="pic")


meanphen_Mtot = cbind(meanphen_MD,meanphen_MV)


bm = mvgls(as.matrix(meanphen_Mtot)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_Mtot)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_Mtot)~1, tree=subtree, model="EB")
lam_Mtot = mvgls(as.matrix(meanphen_Mtot)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_Mtot)

lam_Mtot


subtree_lam3 = transformPhylo(subtree, model="lambda", lambda = lam_Mtot$param)

resbm3 = mvBM(subtree_lam3, as.matrix(meanphen_Mtot), method="pic", param = list(constraint="equal"))
resbm3
mean(sqrt(diag(resbm3$sigma)))


simu_tot = mvSIM(subtree_lam3,nsim=100,model="BM",param=resbm3)

null_ratio = rep(NA, 100)
for (i in c(1:100)){
  print(i)
  dorsal = simu_tot[[i]][,1:11]
  ventral = simu_tot[[i]][,12:22]
  resbm_dorsal = mvBM(subtree_lam2, dorsal, method="pic", echo = F, diagnostic = F)
  resbm_ventral = mvBM(subtree_lam1, ventral, method="pic", echo = F, diagnostic = F)
  null_ratio[i] = mean(diag(resbm_ventral$sigma))/mean(diag(resbm_dorsal$sigma))
  null_ratio[i] = mean(sqrt(diag(resbm_ventral$sigma)))/mean(sqrt(diag(resbm_dorsal$sigma)))
}

real_ratio = mean(diag(resbm$sigma))/mean(diag(resbm2$sigma))
real_ratio
p=sum((null_ratio>real_ratio) & (!is.na(null_ratio)))/sum(!is.na(null_ratio))
p



#--------------- Rate ratio Females -------------


subtree = match.phylo.data(subtree, meanphen_FV)$phy
bm = mvgls(as.matrix(meanphen_FV)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_FV)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_FV)~1, tree=subtree, model="EB")
lam_FV = mvgls(as.matrix(meanphen_FV)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_FV)

lam_FV

subtree_lam1 = transformPhylo(subtree, model="lambda", lambda = lam_FV$param)

resbm = mvBM(subtree_lam1, as.matrix(meanphen_FV), method="pic")


bm = mvgls(as.matrix(meanphen_FD)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_FD)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_FD)~1, tree=subtree, model="EB")
lam_FD = mvgls(as.matrix(meanphen_FD)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_FD)

lam_FD

subtree_lam2 = transformPhylo(subtree, model="lambda", lambda = lam_FD$param)

resbm2 = mvBM(subtree_lam2, as.matrix(meanphen_FD), method="pic")


meanphen_Ftot = cbind(meanphen_FD,meanphen_FV)


bm = mvgls(as.matrix(meanphen_Ftot)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_Ftot)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_Ftot)~1, tree=subtree, model="EB")
lam_Ftot = mvgls(as.matrix(meanphen_Ftot)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_Ftot)

lam_Ftot


subtree_lam3 = transformPhylo(subtree, model="lambda", lambda = lam_Ftot$param)

resbm3 = mvBM(subtree_lam3, as.matrix(meanphen_Ftot), method="pic", param = list(constraint="equal"))

simu_tot = mvSIM(subtree_lam3,nsim=100,model="BM",param=resbm3)

null_ratio = rep(NA, 100)
for (i in c(1:100)){
  print(i)
  dorsal = simu_tot[[i]][,1:11]
  ventral = simu_tot[[i]][,12:22]
  resbm_dorsal = mvBM(subtree_lam2, dorsal, method="pic", echo = F, diagnostic = F)
  resbm_ventral = mvBM(subtree_lam1, ventral, method="pic", echo = F, diagnostic = F)
  null_ratio[i] = mean(diag(resbm_ventral$sigma))/mean(diag(resbm_dorsal$sigma))
  null_ratio[i] = mean(sqrt(diag(resbm_ventral$sigma)))/mean(sqrt(diag(resbm_dorsal$sigma)))
}


real_ratio = mean(diag(resbm$sigma))/mean(diag(resbm2$sigma))
real_ratio
p=sum((null_ratio>real_ratio) & (!is.na(null_ratio)))/sum(!is.na(null_ratio))
p

