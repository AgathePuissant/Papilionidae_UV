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
meanphen_FD <- as.data.frame(datasets[["0"]])  # Rows with no suffix
meanphen_FV <- as.data.frame(datasets[["1"]])
meanphen_MD <- as.data.frame(datasets[["2"]])
meanphen_MV <- as.data.frame(datasets[["3"]])

# Applying the function to your dataset
datasets <- split_datasets(meanphen_grayscale)

# Access the individual datasets
meanphen_FD_grayscale <- as.data.frame(datasets[["0"]])  # Rows with no suffix
meanphen_FV_grayscale <- as.data.frame(datasets[["1"]])
meanphen_MD_grayscale <- as.data.frame(datasets[["2"]])
meanphen_MV_grayscale <- as.data.frame(datasets[["3"]])



# #------------Procrustes-----------------
# meanphen=meanphen[rownames(meanphen_grayscale),]
# # Step 1: Standardize the data
# X_scaled <- scale(meanphen)
# Y_scaled <- scale(meanphen_grayscale)
# 
# # Step 2: Perform CCA
# proc_result <- procrustes(X_scaled, Y_scaled)
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
# meanphen_FD <- as.data.frame(datasets[["0"]])  # Rows with no suffix
# meanphen_FV <- as.data.frame(datasets[["1"]])
# meanphen_MD <- as.data.frame(datasets[["2"]])
# meanphen_MV <- as.data.frame(datasets[["3"]])
# 
# # Applying the function to your dataset
# datasets <- split_datasets(meanphen_grayscale)
# 
# # Access the individual datasets
# meanphen_FD_grayscale <- as.data.frame(datasets[["0"]])  # Rows with no suffix
# meanphen_FV_grayscale <- as.data.frame(datasets[["1"]])
# meanphen_MD_grayscale <- as.data.frame(datasets[["2"]])
# meanphen_MV_grayscale <- as.data.frame(datasets[["3"]])

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
  meanphen_D_vis = meanphen_D_vis[rownames(meanphen_V_vis),,drop=F]
  meanphen_V_vis = meanphen_V_vis[rownames(meanphen_D_vis),,drop=F]
  meanphen_D_vis = meanphen_D_vis[complete.cases(meanphen_D_vis),,drop=F]
  meanphen_V_vis = meanphen_V_vis[complete.cases(meanphen_V_vis),,drop=F]
  
  }else{
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



#------------------ PGLS dorso-ventral -----------------------------------------

#Compute dorso-ventral distance
DV_visible = as.data.frame(diag(as.matrix(dist(meanphen_D_vis, meanphen_V_vis))))
colnames(DV_visible)<-"DV_vis"

DV_grayscale = as.data.frame(diag(as.matrix(dist(meanphen_D_grayscale, meanphen_V_grayscale))))
colnames(DV_grayscale)<-"DV_grayscale"
rownames(DV_grayscale) <- rownames(meanphen_D_grayscale)

DV_visible = DV_visible[rownames(DV_grayscale),,F]

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


