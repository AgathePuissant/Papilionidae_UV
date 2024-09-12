library(ggplot2)
library(ggpmisc)

mean_brightness <- read.csv("./data/mean_brightness.csv", sep=";")
meanb <- read.csv("./data/mean_brightness_before.csv", sep=";")

df_plot = merge(mean_brightness,meanb,by=c("id","view"))


ggplot(df_plot,aes(meanb.x,meanb.y))+
  geom_point()+
  xlab("Mean brightness for each image obtained using the new setup")+
  ylab("Mean brightness for each image obtained using the previous setup")+
  theme_classic()+
  stat_poly_line() +
  stat_poly_eq()

res = lm(df_plot$meanb.y~df_plot$meanb.x)

RSS <- c(crossprod(res$residuals))
MSE <- RSS / length(res$residuals)
RMSE <- sqrt(MSE)
RMSE
