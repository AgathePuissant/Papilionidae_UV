
#-------------------------- Import necessary libraries and functions -----------

library(patternize)
library(ape)
library(recolorize)
library(diverge)
library(raster)
library(stringr)
library(phytools)
library(foreach)
library(doParallel)

source("./code/basis_functions/match_tree.R")
source("./code/basis_functions/get_phenotype.R")
source("./code/patPCA_total.R")

colors_sis_sp <- read.delim2("./data/colors_sis_sp.txt", header=FALSE)

#---------------------- Prepare data for patternize ----------------------------

mode = "visible" #Here change to "visible" to do the analysis on the visible light images

if (mode=="UV"){
  
  list_mask = list.files(path = './UV_pictures/outlines_UV')
  
  IDlist <- list.files("./UV_pictures/UV_pictures_resized")
  IDlist <- gsub(".png","",IDlist)
  
  prepath <- './UV_pictures/UV_pictures_resized'
  extension <- '.png'
  
  imageList <- makeList(IDlist, 'image', prepath, extension)
  
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
  
  list_get_phenotype = get_phenotype(c("F"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV_pythoncalib.csv",path_coords = pc, reduce_dataset=T)
  meanphen <- list_get_phenotype[[1]]
  data_FV <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  rm(list=c("list_get_phenotype"))
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data_FV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
  subtree <- list_match[[1]]
  
  list_get_phenotype = get_phenotype(c("M"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV_pythoncalib.csv",path_coords = pc, reduce_dataset=T)
  meanphen <- list_get_phenotype[[1]]
  data_MD <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  rm(list=c("list_get_phenotype"))
  
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data_MD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
  subtree <- list_match[[1]]
  
  list_get_phenotype = get_phenotype(c("M"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_UV_pythoncalib.csv",path_coords = pc, reduce_dataset=T)
  meanphen <- list_get_phenotype[[1]]
  data_MV <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  rm(list=c("list_get_phenotype"))
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data_MV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
  subtree <- list_match[[1]]
  
  data_grayscale = rbind(data_FD,data_FV,data_MD,data_MV)
  
  
  
  data_patternize = data_grayscale
}else{
  
  list_mask = list.files(path = './visible_pictures/outlines_visible')
  
  IDlist <- list.files("./visible_pictures/visible_pictures_resized")
  IDlist <- gsub(".jpg","",IDlist)
  
  prepath <- './visible_pictures/visible_pictures_resized'
  extension <- '.jpg'
  
  imageList <- makeList(IDlist, 'image', prepath, extension)
  
  pc="./data/pca_embeddings_match_all.csv"
  lvl = "sp"
  if (lvl == "form"){
    adp=T
  }else{
    adp = F
  }
  
  list_get_phenotype = get_phenotype(c("F"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_visible.csv",path_coords = pc, reduce_dataset=T)
  meanphen <- list_get_phenotype[[1]]
  data_FD <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  rm(list=c("list_get_phenotype"))
  
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data_FD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
  subtree <- list_match[[1]]
  
  list_get_phenotype = get_phenotype(c("F"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_visible.csv",path_coords = pc, reduce_dataset=T)
  meanphen <- list_get_phenotype[[1]]
  data_FV <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  rm(list=c("list_get_phenotype"))
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data_FV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
  subtree <- list_match[[1]]
  
  list_get_phenotype = get_phenotype(c("M"),c("D"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_visible.csv",path_coords = pc, reduce_dataset=T)
  meanphen <- list_get_phenotype[[1]]
  data_MD <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  rm(list=c("list_get_phenotype"))
  
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data_MD, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
  subtree <- list_match[[1]]
  
  list_get_phenotype = get_phenotype(c("M"),c("V"), mode = 'mean', level = lvl, path_data_photos = "./data/data_photos_visible.csv",path_coords = pc, reduce_dataset=T)
  meanphen <- list_get_phenotype[[1]]
  data_MV <- list_get_phenotype[[2]]
  sp_data <- list_get_phenotype[[4]]
  rm(list=c("list_get_phenotype"))
  
  list_match <- match_tree(meanphen_match = meanphen, data_match = data_MV, add_poly=adp, tree_path = "./data/Papilionidae_MCC_clean.tre")
  subtree <- list_match[[1]]
  
  data = rbind(data_FD,data_FV,data_MD,data_MV)
  
  data_patternize = data
}

#Set patternize parameters
RGB=rep(220,3) #RGB cutoff 
vec_nk = c() #To store number of color for segmentation of visible light images


#Prepare the dataframe to receive the results
# sis=extract_sisters(subtree)
sis=read.table("./data/sis_list.csv", head=T,sep=";",row.names = 1)
colnames(sis) <- c("sp1","sp2")

listpatternize =list()
sis$ratioM = NA
sis$ratioF = NA

#------------------- Patternize analysis ---------------------------------------

for (i in c(1:length(sis$sp1))){
  
  print(i)
  
  #Get the number of colors to segment the images of sister species
  colors_sp = colors_sis_sp[(colors_sis_sp$V1==sis$sp1[i]) & (colors_sis_sp$V2==sis$sp2[i]),]$V3
  colors_sp = str_split(colors_sp, ", ")[[1]]
  nk=length(colors_sp)
  
  #Get the corresponding image id 
  idsp1 = data_patternize[data_patternize$genresp == sis$sp1[i],]$id
  idsp2 = data_patternize[data_patternize$genresp == sis$sp2[i],]$id
  idsp1_dv = paste0(data_patternize[data_patternize$genresp == sis$sp1[i],]$id,data_patternize[data_patternize$genresp == sis$sp1[i],]$view)
  idsp2_dv = paste0(data_patternize[data_patternize$genresp == sis$sp2[i],]$id,data_patternize[data_patternize$genresp == sis$sp2[i],]$view)
  
  
  nsp1 = which(gsub("[D|V]","",IDlist) %in% idsp1)
  nsp2 = which(gsub("[D|V]","",IDlist) %in% idsp2)
  
  
  
  sublist = imageList[c(nsp1,nsp2)] #Get images from the list
  target <- sublist[[1]] #Define the first image as the target image to align the other wings
  
  
  
  if (mode=="UV"){
    
    foutline <-read.table(paste0('./UV_pictures/outlines_UV/',list_mask[grep(gsub("_[1-9]","",names(target)[1]), list_mask)]), h= F) #Read the corresponding outline to remove background
    
    #Aligne the wings using patternize
    foutline[,1] = xmax(sublist[[1]]) - foutline[,1]
    rasterList_regRGB <- patRegRGB(sublist, target, RGB, plot = "stack", removebgR = 254, maskOutline = foutline, useBlockPercentage = 100)
    patternize_list <- rasterList_regRGB
    
    #sample ids of dorsal and ventral side for each species
    Dsp1 = idsp1_dv[grep("D", idsp1_dv)]
    Vsp1 = idsp1_dv[grep("V", idsp1_dv)]
    Dsp2 = idsp2_dv[grep("D", idsp2_dv)]
    Vsp2 = idsp2_dv[grep("V", idsp2_dv)]
    Dsp1 = Dsp1[!is.na(Dsp1)]
    Vsp1 = Vsp1[!is.na(Vsp1)]
    Dsp2 = Dsp2[!is.na(Dsp2)]
    Vsp2 = Vsp2[!is.na(Vsp2)]
    
    #Compute PCA on binary dataframes and plot the output "morpho-space"
    pcaOut <- patPCA(rasterList_regRGB,
                     list(Dsp1,Vsp1,Dsp2,Vsp2),
                     c("darkblue","darkred","blue","red"),
                     plot=F,#Change to T to plot
                     plotChanges = F, #Change to T to plot the area of variation corresponding to the axis
                     normalized = TRUE)
    
    
    
  }
  if (mode=="visible"){
    
    foutline <-read.table(paste0('./visible_pictures/outlines_visible/',list_mask[grep(gsub("_[1-9]","",names(target)[1]), list_mask)]), h= F) #Read the corresponding outline to remove background
    
    #Align the wings using patternize
    rasterList_reg <- alignReg(sublist, target, plotTransformed = F, removebgR = 254, maskOutline = foutline, useBlockPercentage = 100)
    
    
    #This workflow to use recolorize to impose a color palette for further patternize analysis is taken from :
    # https://hiweller.rbind.io/post/recolorize-patternize-workflow/
    
    # convert from RasterBricks to image arrays using the brick_to_array function:
    imgs <- lapply(rasterList_reg, brick_to_array)
    names(imgs) <- names(rasterList_reg)
    
    # save raster extents for later conversion:
    extent_list <- lapply(rasterList_reg, raster::extent)
    
    # make an empty list for storing the recolorize objects
    rc_list <- vector("list", length(imgs))
    names(rc_list) <- names(imgs)
    
    # for every image, run the same recolorize2 function to fit a recolorize object:
    for (q in 1:length(imgs)) {
      rc_list[[q]] <- recolorize2(imgs[[q]], bins = 3,
                                  cutoff = 35, plotting = F)
    }
    
    # get a dataframe of all colors:
    all_palettes <- do.call(rbind, lapply(rc_list, function(r) r$centers))
    
    # and for cluster sizes (as a proportion of their original image):
    all_sizes <- do.call(c, lapply(rc_list, function(s) s$sizes))
    
    # plot colors using hclust and return grouping list:
    par(mar = rep(2, 4))
    cluster_list <- hclust_color(all_palettes, n_final = nk)
    
    
    # make an empty matrix for storing the new palette
    butterfly_palette <- matrix(NA, ncol = 3, nrow = length(cluster_list))
    
    # for every color in cluster_list...
    for (o in 1:length(cluster_list)) {
      
      # get the center indices
      idx <- cluster_list[[o]]
      
      # get the average value for each channel, using cluster size to get a weighted average
      ctr <- apply(all_palettes, 2, 
                   function(p) weighted.mean(p[idx], 
                                             w = all_sizes[idx]))
      
      # store in the palette matrix
      butterfly_palette[o, ] <- ctr
    }
    
    impose_list <- lapply(imgs, function(k) imposeColors(k, butterfly_palette, 
                                                         adjust_centers = F, 
                                                         plotting = F))
    
    # convert to patternize:
    patternize_list <- lapply(impose_list, recolorize_to_patternize)
    
    # and set extents again:
    for (m in 1:length(patternize_list)) {
      for (n in 1:length(patternize_list[[1]])) {
        raster::extent(patternize_list[[m]][[n]]) <- extent_list[[m]]
      }
    }
    
    rasterList_regK <- patternize_list
    
    #The paragraph below can be used to visualize the area of variation on the images for each color
    
    # summedRaster_regK <- sumRaster(rasterList_regK, c(idsp1_dv,idsp2_dv), type = 'k')
    # summedRaster_masked <- lapply(summedRaster_regK, function(x) maskOutline(x, foutline, refShape = 'target', flipOutline = 'y', imageList = imageList))
    # 
    # library(viridis)
    # colfunc <- inferno(100)
    # plotHeat(summedRaster_masked, c(idsp1_dv,idsp2_dv), refShape = 'target', imageList = imageList, colpalette = colfunc)
    
    pcaOut <- patPCA_total(patternize_list, quietly = FALSE)
  }
  
  
  # first, make a blank plot
  PCx <- 1; PCy <- 2
  pca_summary <- summary(pcaOut)
  limits <- apply(pcaOut$x[ , c(PCx, PCy)], 2, range)
  par(mar = c(4, 4, 2, 1))
  
  #Get the grouping factors for the plot
  idview = rownames(pcaOut$x)
  mask = match(idview,paste0(data_patternize$id,data_patternize$view))
  idview = paste0(data_patternize[mask,]$view,data_patternize[mask,]$sex)
  idsp = as.factor(paste0(data_patternize[mask,]$tipsgenre))
  idview=as.factor(idview)
  colorspal = c("darkblue", "darkred","blue","red") #Define the color for each grouping factor
  names(colorspal) <- c("DM","DF","VM","VF") #D: dorsal, V: ventral, M: males, F:females
  idview = colorspal[match(idview,names(colorspal))]
  
  shapepal = c("16","17")
  names(shapepal) = levels(idsp)
  idsp = shapepal[match(idsp,names(shapepal))]
  
  plot(pcaOut$x[ , c(PCx, PCy)], type = "n",
       asp = 1,
       xlim = limits[ , 1] + c(-5, 5), 
       ylim = limits[ , 2] + c(-10, 10),
       xlab=paste0('PC1 (', round(pca_summary$importance[2, PCx]*100, 1), ' %)'),
       ylab=paste0('PC2 (', round(pca_summary$importance[2, PCy]*100, 1), ' %)'))
  
  # then add images:
  for (p in c(1:dim(pcaOut$x)[1])) {
    
    if (mode=="visible"){
    add_image(impose_list[[p]]$original_img,
              x = pcaOut$x[p, PCx],
              y = pcaOut$x[p, PCy],
              width = 20)
    }
    
    points(x = pcaOut$x[p, PCx]-5,
           y = pcaOut$x[p, PCy]-5,
           pch = as.numeric(idsp[p]),
           col = idview[p],
           cex=1.5)
  }
  
  
  #Get mean coordinates per side of the wings, sex and species
  data_patternize$idview = paste0(data_patternize$id,data_patternize$view)
  coord_pca = as.data.frame(pcaOut$x)
  coord_pca$idview = names(sublist)
  coord_pca = merge(coord_pca, data_patternize, by="idview")
  mean_coord_pca = coord_pca %>% group_by(sex,view,tipsgenre) %>% summarise_if(is.numeric, mean, na.rm = TRUE) 
  
  #Compute distances between each group of sex/species/wing side
  dist_coord_pca = as.matrix(dist(mean_coord_pca[,grep("PC",colnames(mean_coord_pca))]))
  
  #Label the two species as species 1 or 2, sex and wing side
  mean_coord_pca$nsp =  cut(as.numeric(mean_coord_pca$tipsgenre), 2, labels=c("1","2"))
  colnames(dist_coord_pca) <- paste0(mean_coord_pca$nsp,mean_coord_pca$view,mean_coord_pca$sex)
  rownames(dist_coord_pca) <- paste0(mean_coord_pca$nsp,mean_coord_pca$view,mean_coord_pca$sex)
  
  #Compute ratio of divergence of ventral sides over divergence of dorsal sides between species for each sex
  sis[i,]$ratioM <- dist_coord_pca["2VM","1VM"]/dist_coord_pca["2DM","1DM"]
  sis[i,]$ratioF <- dist_coord_pca["2VF","1VF"]/dist_coord_pca["2DF","1DF"]
  
  #Add the results of the PCA to a list for further analyses
  listpatternize = append(listpatternize,list(patternize_list, pcaOut))

  
}


#------------------- Permutation tests -----------------------------------------

# Preprocessing of the list of results
coordpca_list <- lapply(seq(2, length(listpatternize), 2), function(pca_n) {
  coordpca_perm <- as.data.frame(listpatternize[[pca_n]]$x)
  coordpca_perm$idview <- rownames(coordpca_perm)
  coordpca_perm <- merge(coordpca_perm, data_patternize, by = "idview")
  coordpca_perm
})



# Parallelized computation
# Initialize parallel backend
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

# Permutations
nsim <- 1000 #Set the numbers of permutations

#Prepare the dataframe for the output of each permutation
perm_df <- matrix(nrow = nsim, ncol = length(coordpca_list))

perm_df <- foreach(permutation = 1:nsim, .combine = "rbind") %dopar% {
  
  results <- numeric(length(coordpca_list))
  
  for (pca_n in seq_along(coordpca_list)) {
    
    coordpca_perm <- coordpca_list[[pca_n]] #For each pair of sister species get the PCA
    
    #Randomly permute dorsal and ventral coordinates of individuals
    coin_toss <- sample(c(0, 1), nrow(coordpca_perm), replace = TRUE) #Toss a coin for each PCA coordinate
    
    for (row in which(coin_toss == 1)) {
      storage <- coordpca_perm[row, "view"] #Store the wing side (dorsal or ventral) of the focal coordinate
      
      #Get the id of the opposite wing side so that the side is not flipped twice
      id_otherview <- coordpca_perm$tipsgenre == coordpca_perm[row, "tipsgenre"] &
        coordpca_perm$sex == coordpca_perm[row, "sex"] &
        coordpca_perm$view != coordpca_perm[row, "view"]
      
      #Replace the wing side of the focal coordinate by the opposite one
      coordpca_perm[row, "view"] <- unique(coordpca_perm[id_otherview, "view"])
      
      # Replace the wing side of the originally opposite wing side
      coordpca_perm[id_otherview, "view"] <- rep(storage,sum(id_otherview))
      
      #Set the coin toss of the originally opposite wing side to 0 so that it is not flipped twice
      coin_toss[id_otherview]=0
    }
    
    
    # Calculate means for columns containing "PC" in their names
    mean_coord_pca_perm <- aggregate(coordpca_perm[, grep("PC", names(coordpca_perm))],
                                     by = list(sex = coordpca_perm$sex, view = coordpca_perm$view, tipsgenre = coordpca_perm$tipsgenre),
                                     FUN = mean, na.rm = TRUE)
    mean_coord_pca_perm <- mean_coord_pca_perm[order(mean_coord_pca_perm$view), ]
    
    #Calculate distances and compute ratios for females and males
    dist_coord_pca_perm <- as.matrix(dist(mean_coord_pca_perm[, grep("PC", colnames(mean_coord_pca_perm))]))
    
    mean_coord_pca_perm$nsp <- cut(as.numeric(mean_coord_pca_perm$tipsgenre), 2, labels = c("1", "2"))
    colnames(dist_coord_pca_perm) <- paste0(mean_coord_pca_perm$nsp, mean_coord_pca_perm$view, mean_coord_pca_perm$sex)
    rownames(dist_coord_pca_perm) <- paste0(mean_coord_pca_perm$nsp, mean_coord_pca_perm$view, mean_coord_pca_perm$sex)
    
    
    ratioM_perm <- dist_coord_pca_perm["2VM", "1VM"] / dist_coord_pca_perm["2DM", "1DM"]
    ratioF_perm <- dist_coord_pca_perm["2VF", "1VF"] / dist_coord_pca_perm["2DF", "1DF"]
    
    #Uncomment the line for the sex for which you want to plot the permutation results
    results[pca_n] <- ratioF_perm
    # results[pca_n] <- ratioF_perm
  }
  results
}

# Clean up parallel backend
stopCluster(cl)

#Compute median ratios for each permutation
perm_df = as.data.frame(perm_df)
median_perm = perm_df %>%
  rowwise() %>%
  mutate(row_median = median(c_across(where(is.numeric)), na.rm=TRUE))


#-------------------- Plot permutations and compute p-values -------------------

#Uncomment the line for the sex of which you want to plot the permutations results
stattotest = median(sis$ratioF, na.rm=T)
# stattotest = median(sis$ratioF, na.rm=T)


pval = sum(median_perm$row_median>stattotest)/nsim 
pval

ggplot(median_perm, aes(x=row_median))+
  geom_density(fill="darkblue", alpha=0.5, bw=0.01)+
  theme_classic()+
  geom_vline(xintercept=stattotest, size=1.3, linetype="solid", color = "red")

#----------------------------------