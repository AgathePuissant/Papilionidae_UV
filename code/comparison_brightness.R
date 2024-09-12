library(ggplot2)
library(ggpmisc)

mean_brightness <- read.csv("./data/mean_brightness.csv", sep=";")
meanb <- read.csv("./data/old/meanb.csv", sep=";")

df_plot = merge(mean_brightness,meanb,by=c("id","view"))


ggplot(df_plot,aes(meanb.x/255,meanb.y))+
  geom_point()+
  xlab("Mean brightness for each image obtained using the new setup")+
  ylab("Mean brightness for each image obtained using the previous setup")+
  theme_classic()+
  stat_poly_line() +
  stat_poly_eq()

