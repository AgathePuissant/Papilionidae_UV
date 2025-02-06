#------------------ Prepare libraries and functions ----------------------------

library(diverge)
library(reshape2)
library(ggpubr)
library(ape)
library(ggpubr)
library(phytools)
library(RColorBrewer)
library(lme4)
library(AICcmodavg)
library(multcomp)

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

ages = read.csv("./data/collection_date.csv",sep=";")


data_grayscale <- merge(data_grayscale,meanb, by = c("id","view"))
data_grayscale$collection_date = ages[match(data_grayscale$id,ages$id),]$collection_date

#------------------ Effect of specimen age --------------------------------


data_lme = data_grayscale

data_lme$specimen_age = 2025 - data_lme$collection_date

data_lme = data_lme[!is.na(data_lme$collection_date),,drop=F]

model_mixed <- lmer(meanb~specimen_age+(1|tipsgenre),
                    data=data_lme,
                    REML = F)

summary(model_mixed) 

null_model = lmer(meanb~1+(1|tipsgenre),
                  data=data_lme,
                  REML= F)

AICcmodavg::AICc(null_model, return.K = F, second.ord=F) - AICcmodavg::AICc(model_mixed, return.K = F, second.ord=F)

confint(model_mixed)

#--------- phylogenetic ANCOVA to test for differences in brightness while accounting for specimen age --------------

data_lme$group = as.factor(paste0(data_lme$sex,data_lme$view))

corBM<-corBrownian(phy=subtree,form=~tipsgenre)

brightness.ancova<-gls(meanb~group+specimen_age,data=data_lme,
                    correlation=corBM)

summary(brightness.ancova)

post.hoc<-glht(brightness.ancova,linfct=mcp(group="Tukey"))
summary(post.hoc)
