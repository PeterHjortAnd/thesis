## Samlet KI
# Run after all other have been run
estall <- matrix(rep(0,120),nrow= 30, ncol =5)
for (i in 1:5){
  estall[,i] <- c(est15[,i],est16[,i],est17[,i],est19[,i],est20[,i]) 
}
estall_data <- as.data.frame(estall)
colnames(estall_data) <- c(c('model', 'est','lower','upper','year'))

sp <- ggplot(estall_data, aes(x = model, y = est, color = model)) +        # ggplot2 plot with confidence intervals
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  labs(x = "Models", y = "Estimates of non-observed bears") +
  scale_x_continuous(breaks = 1:6)
sp
sp + facet_grid(. ~ year) + theme_bw(base_size = 18)
