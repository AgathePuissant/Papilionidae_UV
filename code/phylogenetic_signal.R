
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

library(motmot)
subtree_lam1 = transformPhylo(subtree, model="lambda", lambda = lam_MV$param)

resbm = mvBM(subtree_lam1, as.matrix(meanphen_MV_grayscale), method="pic")
resbm
mean(diag(resbm$sigma))


bm = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="EB")
lam_MD = mvgls(as.matrix(meanphen_MD_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_MD)

lam_MD

library(motmot)
subtree_lam2 = transformPhylo(subtree, model="lambda", lambda = lam_MD$param)

resbm2 = mvBM(subtree_lam2, as.matrix(meanphen_MD_grayscale), method="pic")
resbm2
mean(diag(resbm2$sigma))



mean(diag(resbm$sigma))/mean(diag(resbm2$sigma))

#-----------------------------------

# Set up libraries
library(geiger) # For phylogenetic operations
library(motmot) # For evolutionary rate calculations
library(mvMORPH) # For multivariate Brownian Motion and OU models

# Define parameters
n_perm <- 100  # Number of permutations
observed_ratio <- mean(diag(resbm$sigma)) / mean(diag(resbm2$sigma))  # Observed ratio

# Initialize storage for null ratios
null_ratios <- numeric(n_perm)



# Permutation loop
for (i in 1:n_perm) {
  print(i)
  # Permute phenotypic data
  permuted_meanphen_MV <- meanphen_MV_grayscale[sample(1:nrow(meanphen_MV_grayscale)), ]
  rownames(permuted_meanphen_MV) = rownames(meanphen_MV_grayscale)
  permuted_meanphen_MD <- meanphen_MD_grayscale[sample(1:nrow(meanphen_MD_grayscale)), ]
  rownames(permuted_meanphen_MD) = rownames(meanphen_MD_grayscale)
  
  # Recalculate lambda-transformed trees
  perm_tree_mv <- transformPhylo(subtree, model="lambda", lambda = lam_MV$param)
  perm_tree_md <- transformPhylo(subtree, model="lambda", lambda = lam_MD$param)
  
  # Recalculate evolutionary rates
  perm_resbm_mv <- mvBM(perm_tree_mv, as.matrix(permuted_meanphen_MV), method="pic")
  perm_resbm_md <- mvBM(perm_tree_md, as.matrix(permuted_meanphen_MD), method="pic")
  
  # Calculate ratio of rates and store
  null_ratios[i] <- mean(diag(perm_resbm_mv$sigma)) / mean(diag(perm_resbm_md$sigma))
}

# Compare observed ratio to null distribution
p_value <- sum(null_ratios[null_ratios>0] >= observed_ratio) / n_perm

# Output results
cat("Observed ratio:", observed_ratio, "\n")
cat("P-value:", p_value, "\n")

# Visualize null distribution
hist(null_ratios, main="Null Distribution of Ratios", xlab="Rate Ratio", col="gray")
abline(v = observed_ratio, col = "red", lwd = 2, lty = 2)  # Add observed ratio line


#----------------------------


subtree = match.phylo.data(subtree, meanphen_FV_grayscale)$phy
bm = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="EB")
lam_FV = mvgls(as.matrix(meanphen_FV_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_FV)

lam_FV

library(motmot)
subtree_lam1 = transformPhylo(subtree, model="lambda", lambda = lam_FV$param)

resbm = mvBM(subtree_lam1, as.matrix(meanphen_FV_grayscale), method="pic")
resbm
mean(diag(resbm$sigma))


bm = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="BM")
ou = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="OU")
eb = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="EB")
lam_FD = mvgls(as.matrix(meanphen_FD_grayscale)~1, tree=subtree, model="lambda")


GIC(bm) ; GIC(ou) ; GIC(eb) ; GIC(lam_FD)

lam_FD

library(motmot)
subtree_lam2 = transformPhylo(subtree, model="lambda", lambda = lam_FD$param)

resbm2 = mvBM(subtree_lam2, as.matrix(meanphen_FD_grayscale), method="pic")
resbm2
mean(diag(resbm2$sigma))



mean(diag(resbm$sigma))/mean(diag(resbm2$sigma))


lam_MD$param ; lam_MV$param
lam_FD$param ; lam_FV$param

