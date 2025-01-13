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

pc = "./data/pca_embeddings_pythoncalib_retrained.csv"
lvl = "sp"
if (lvl == "form"){
  adp=T
}else{
  adp = F
}

path_data_UV = "./data/data_photos_UV.csv"

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

subtree = match.phylo.data(subtree,meanphen_MD_grayscale)$phy


#------------------ Get the specimens data (visible) -------------------------------------

pc="./data/pca_embeddings_visible.csv"

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




#-------------------CCA--------------------------------------------
meanphen=meanphen[rownames(meanphen_grayscale),]
# Step 1: Standardize the data
X_scaled <- scale(meanphen)
Y_scaled <- scale(meanphen_grayscale)

# Step 2: Perform CCA
cca_result <- cancor(X_scaled, Y_scaled)

# # Step 3: Display results
# cat("Canonical Correlations:\n")
# print(cca_result$cor)

# Step 4: Project data into the canonical space
X_proj <- X_scaled %*% cca_result$xcoef  # Transformed X
Y_proj <- Y_scaled %*% cca_result$ycoef  # Transformed Y

# Step 5: Compare distances in the aligned space
dist_X <- dist(X_proj)  # Pairwise distances in aligned space 1
dist_Y <- dist(Y_proj)  # Pairwise distances in aligned space 2

# # Example: Compute correlation between distance matrices
# distance_correlation <- cor(as.vector(dist_X), as.vector(dist_Y))
# cat("Distance Correlation between aligned spaces:", distance_correlation, "\n")

meanphen = X_proj
meanphen_grayscale = Y_proj

# Split the dataset into 4 based on the suffix pattern (no suffix, then 1, 2, 3)
split_datasets <- function(dataset) {
  # Identify rows with no suffix or numeric suffix
  row_suffixes <- ifelse(grepl("\\d$", rownames(dataset)),
                         sub(".*(\\d)$", "\\1", rownames(dataset)),
                         "0")

  # Split the dataset based on the suffix
  split_data <- split(as.data.frame(dataset), row_suffixes)

  # Process each subset: ensure it's a data frame and clean rownames/colnames
  result <- lapply(split_data, function(df) {
    df <- as.matrix(df)  # Ensure the subset is a matrix
    # Remove the numeric suffix from rownames
    rownames(df) <- gsub("\\d$", "", rownames(df))
    # Rename columns to coord0, coord1, ...
    colnames(df) <- paste0("coord", seq_len(ncol(df)) - 1)
    return(df)
  })

  return(result)
}

# Applying the function to your dataset
datasets <- split_datasets(meanphen)

# Access the individual datasets
meanphen_FD <- datasets[["0"]]  # Rows with no suffix
meanphen_FV <- datasets[["1"]]
meanphen_MD <- datasets[["2"]]
meanphen_MV <- datasets[["3"]]

# Applying the function to your dataset
datasets <- split_datasets(meanphen_grayscale)

# Access the individual datasets
meanphen_FD_grayscale <- datasets[["0"]]  # Rows with no suffix
meanphen_FV_grayscale <- datasets[["1"]]
meanphen_MD_grayscale <- datasets[["2"]]
meanphen_MV_grayscale <- datasets[["3"]]



# #------------Procrustes-----------------
# meanphen=meanphen[rownames(meanphen_grayscale),]
# # Step 1: Standardize the data
# X_scaled <- scale(meanphen)
# Y_scaled <- scale(meanphen_grayscale)
# 
# # Step 2: Perform CCA
# proc_result <- procrustes(X_scaled, Y_scaled, scale=F)
# 
# # # Step 3: Display results
# # cat("Canonical Correlations:\n")
# # print(cca_result$cor)
# 
# # Step 4: Project data into the canonical space
# X_proj <- proc_result$X  # Transformed X
# Y_proj <- proc_result$Y  # Transformed Y
# 
# # Step 5: Compare distances in the aligned space
# dist_X <- dist(X_proj)  # Pairwise distances in aligned space 1
# dist_Y <- dist(Y_proj)  # Pairwise distances in aligned space 2
# 
# # # Example: Compute correlation between distance matrices
# # distance_correlation <- cor(as.vector(dist_X), as.vector(dist_Y))
# # cat("Distance Correlation between aligned spaces:", distance_correlation, "\n")
# 
# meanphen = X_proj
# meanphen_grayscale = Y_proj
# 
# # Split the dataset into 4 based on the suffix pattern (no suffix, then 1, 2, 3)
# split_datasets <- function(dataset) {
#   # Identify rows with no suffix or numeric suffix
#   row_suffixes <- ifelse(grepl("\\d$", rownames(dataset)),
#                          sub(".*(\\d)$", "\\1", rownames(dataset)),
#                          "0")
# 
#   # Split the dataset based on the suffix
#   split_data <- split(as.data.frame(dataset), row_suffixes)
# 
#   # Process each subset: ensure it's a data frame and clean rownames/colnames
#   result <- lapply(split_data, function(df) {
#     df <- as.matrix(df)  # Ensure the subset is a matrix
#     # Remove the numeric suffix from rownames
#     rownames(df) <- gsub("\\d$", "", rownames(df))
#     # Rename columns to coord0, coord1, ...
#     colnames(df) <- paste0("coord", seq_len(ncol(df)) - 1)
#     return(df)
#   })
# 
#   return(result)
# }
# 
# # Applying the function to your dataset
# datasets <- split_datasets(meanphen)
# 
# # Access the individual datasets
# meanphen_FD <- datasets[["0"]]  # Rows with no suffix
# meanphen_FV <- datasets[["1"]]
# meanphen_MD <- datasets[["2"]]
# meanphen_MV <- datasets[["3"]]
# 
# # Applying the function to your dataset
# datasets <- split_datasets(meanphen_grayscale)
# 
# # Access the individual datasets
# meanphen_FD_grayscale <- datasets[["0"]]  # Rows with no suffix
# meanphen_FV_grayscale <- datasets[["1"]]
# meanphen_MD_grayscale <- datasets[["2"]]
# meanphen_MV_grayscale <- datasets[["3"]]

#------------------ Compute pairwise distance in UV ----------------------------

sis=extract_sisters(subtree)

MD = create_distpheno(meanphen_MD_grayscale,"M","D", level=lvl)
MD$sex="M"
MD$view="D"
MV = create_distpheno(meanphen_MV_grayscale,"M","V", level=lvl)
MV$sex="M"
MV$view="V"
FD = create_distpheno(meanphen_FD_grayscale,"F","D", level=lvl)
FD$sex="F"
FD$view="D"
FV = create_distpheno(meanphen_FV_grayscale,"F","V", level=lvl)
FV$sex="F"
FV$view="V"

M = rbind(MD,MV)
Fem = rbind(FD,FV)

MF_grayscale=rbind(M,Fem)
MF_grayscale$type = "grayscale"

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

# MF_grayscale$value <- scale(MF_grayscale$value)
# MF_grayscale$visdist <- scale(MF_grayscale$visdist)


# MF_grayscale = MF_grayscale[MF_grayscale$Var1!="Protesilaus_molops",,F]

#GLM
model<-glm(value~(visdist+overlap+distphylo):sex:view, data=MF_grayscale)

summary(model)


#-------------------------------------------------------------------------------

ggplot(MF_grayscale[(MF_grayscale$view=="V") & (MF_grayscale$sex=="M"),],aes(x=visdist, y=value, label=Var1))+
  geom_point()+
  geom_text(size=2)


# id_toretain = paste0(data_UV[((data_UV$tipsgenre %in% sis$sp1) | (data_UV$tipsgenre %in% sis$sp2)) & (data_UV$sex == "M"),]$id, "V.jpg")


