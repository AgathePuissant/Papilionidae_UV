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



#------------------ Compute pairwise distance in UV ----------------------------

sis=extract_sisters(subtree)

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

#-------------------------------------------------------


library(readr)
data_photos <- read_delim("data/data_photos_UV.csv", 
                                  delim = ";", escape_double = FALSE, trim_ws = TRUE)
data_photos$tipsgenre = paste0(data_photos$genus,"_",data_photos$sp)
meanb <- read.csv("./data/meanb.csv", sep=";")
meanbdf <- merge(data_photos,meanb, by = c("id"), all.y=T)
meanbdf = meanbdf[!is.na(meanbdf$tipsgenre),,F]
meanbdf =meanbdf %>% 
  group_by(tipsgenre,view,sex) %>%
  summarise(avg = mean(meanb))
meanbdf2 <- meanbdf[meanbdf$sex=="M" & meanbdf$view=="D",]

male_D = meanbdf2$avg
male_D <- as.data.frame(male_D)
male_D$id <- meanbdf2$tipsgenre

meanbdf2 <- meanbdf[meanbdf$sex=="M" & meanbdf$view=="V",]

male_V = meanbdf2$avg
male_V <- as.data.frame(male_V)
male_V$id <- meanbdf2$tipsgenre

meanbdf2 <- meanbdf[meanbdf$sex=="F" & meanbdf$view=="D",]

female_D = meanbdf2$avg
female_D <- as.data.frame(female_D)
female_D$id <- meanbdf2$tipsgenre

meanbdf2 <- meanbdf[meanbdf$sex=="F" & meanbdf$view=="V",]

female_V = meanbdf2$avg
female_V <- as.data.frame(female_V)
female_V$id <- meanbdf2$tipsgenre

library(tidyverse)

#put all data frames into list
df_list_M <- list(male_D,male_V)
df_list_F <- list(female_D,female_V)

#merge all data frames in list
brightall_M <- df_list_M %>% reduce(full_join, by="id") 
brightall_M = brightall_M[complete.cases(brightall_M),]
rownames(brightall_M)<-brightall_M$id
brightall_M$id<-NULL
colnames(brightall_M) <- c("MD","MV")

brightall_F <- df_list_F %>% reduce(full_join, by="id") 
brightall_F = brightall_F[complete.cases(brightall_F),]
rownames(brightall_F)<-brightall_F$id
brightall_F$id<-NULL
colnames(brightall_F) <- c("FD","FV")


# brightall_M <- log(1+brightall_M)
# brightall_F <- log(1+brightall_F)
brightall_M$spnames <- rownames(brightall_M)
brightall_F$spnames <- rownames(brightall_F)
keep = intersect(rownames(brightall_F),rownames(brightall_M))
brightall_F = brightall_F[(rownames(brightall_F) %in% keep),,F]
brightall_M = brightall_M[(rownames(brightall_M) %in% keep),,F]

brightness <- cbind(brightall_F,brightall_M)
brightness2 <- melt(brightness)
brightness2 <- brightness2[,-1]

library(ggtree)
library(pals)  
library(phytools)
library(phangorn)


overlap = read.csv(file.path(getwd(),"data/jaccard.csv"),sep=";",row.names=1,dec=",")
sis$overlap = overlap[match(paste(sis$sp1,sis$sp2),paste(overlap$Var1,overlap$Var2)),3]
sis <- sis[!is.na(sis$overlap),,drop=F]
sis$mrca <- apply(sis, 1, function(x) findMRCA(subtree, c(x[1],x[2])))

desc = apply(sis, 1, function(x) getDescendants(subtree, x[4]))
desc = lapply(desc, function(x) if(length(x)<2){ return(Descendants(as.numeric(x),subtree))} else{return(x)})
desc=matrix(unlist(desc),ncol=2,byrow=T)

sis$descendants1 <- desc[,1]
sis$descendants2 <- desc[,2]
sympatry2 <- data.frame(nodenb = c(sis$descendants1,sis$descendants2), val = c(sis$overlap,sis$overlap))
sympatry2 <- sympatry2[complete.cases(sympatry2),,drop=F]
# sympatry2$symp <- sympatry[sympatry2$spnames,]


source("./code/basis_functions/plot_genus.R")
p3 <- plot_genres(subtree, layout_fan=F, palette_grey=T,sizebranch=0.1)
p3

p3 <- p3 + 
  new_scale_fill()+
  geom_highlight(data=sympatry2, mapping=aes(node=nodenb, fill=val), alpha = 0.7, extendto=63)+ 
  scale_fill_gradient(low="#B2F5BE", high = "#124E1D")+ xlim(NA, 250)
p3

library(ggtreeExtra)



p4 <- p3 +
  
  new_scale_fill()+
  
  geom_fruit(
    data = brightall_F,
    geom = geom_bar,
    stat = "identity",
    mapping = aes(x = -FD, y = spnames, fill = FD), 
    orientation = "y",
    offset = 0.7,
    axis.params = list(axis = "xy", line.size = 0.5, nbreak = 5, line.colour = "black", text.size = 2, limits=c(-1,0)),
    width = 0.5,
    grid.params = list()
  ) +
  
  geom_fruit(
    data = brightall_F,
    geom = geom_bar,
    stat = "identity",
    mapping = aes(x = FV, y = spnames, fill = FV), 
    orientation = "y",
    offset = -0.19,
    axis.params = list(axis = "xy", line.size = 0.5, nbreak = 5, line.colour = "black", text.size = 2,limits=c(0,1)),
    width = 0.5,
    grid.params = list()
  ) +
  
  geom_fruit(
    data = brightall_M,
    geom = geom_bar,
    stat = "identity",
    mapping = aes(x = -MD, y = spnames, fill = MD), 
    orientation = "y",
    offset = 0.6,
    axis.params = list(axis = "xy", line.size = 0.5, nbreak = 5, line.colour = "black", text.size = 2,limits=c(-1,0)),
    width = 0.5,
    grid.params = list()
  ) +
  
  geom_fruit(
    data = brightall_M,
    geom = geom_bar,
    stat = "identity",
    mapping = aes(x = MV, y = spnames, fill = MV), 
    orientation = "y",
    offset = -0.19,
    axis.params = list(axis = "xy", line.size = 0.5, nbreak = 5, line.colour = "black", text.size = 2, limits=c(0,1)),
    width = 0.5,
    grid.params = list()
  ) +
  
  # Uniform color scale for all bars
  scale_fill_viridis_c(option = "magma") 

# +
#   
#   # Set limits for both positive and negative values
#   scale_x_continuous(
#     limits = c(-1,1),  # Adjust these limits as per your data range
#     breaks = seq(-1, 1, by = 0.1),
#     expand = c(0, 0)
#   )

p4
