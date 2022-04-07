
library(viridis)

k <- 100
x<-df1$x
y<-df1$y
contour_cols <- viridis(k, alpha = 0.5)
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

dens <- get_density(x, y, k)
## 2. Generating plots -----
def.par<-par()
#png(paste0(cp ,data1," vs ",data2," - ",i1," - ",i2, "Protein ratios.png"), width=4,height=4,units="in",res=300)
par(mar=c(5,6,2,2))
plot(x, y,
     col = contour_cols[findInterval(dens, seq(0, max(dens), length.out = k))], pch = 16, xlab="", ylab="",
     #xlim=c(-2.3,3.5), ylim=c(-1.3,2.5),
     cex.axis=2)


text(-2.5, 2.30, paste0(i1,"/",i2, " (",sum(!is.na(data.RA.l[,2])),")"),pos = 4, cex=1.75)
text(-2.5, 1.80, paste0("?? = ", round(cor(x,y),2)), pos = 4, cex=2)
mtext(paste0(data1,", log2"), side=2, padj = -3, cex=1.75)
mtext(paste0(data2,", log2"), side=1, padj = 3, cex=1.75)
#abline(a=0, b=1)
abline(lm(y~x), col="black", lty=2)
#print(lm(y~x))
#dev.off()

contour_cols

options(repr.plot.width = 1, repr.plot.height = 8)   
p_xq<- ggplot(df1) +
  geom_point(aes(x=x, y=y), 
             col = contour_cols[findInterval(dens, seq(0, max(dens), length.out = k))], 
             pch = 16) +
  geom_line(data=predictions_quantiles_m4_long, 
            aes(x=x, y=value, group=variable),
            color='#082133') +
  labs(title='GAMLSS',
       # subtitle='Location, scale, and shape\nare modeled as functions of x',
       x = 'Age', 
       y='Secret brain region')  + 
  annotate("text", x = 35, y = 30, size=5,
           label = expression(paste(
             'y ~ D(',mu,',',sigma,',',nu,',',tau,')',
             )))+
  ylim(c(20,60))+
  theme_minimal(); 
p_xq
pngfile <- fs::path(knitr::fig_path(),  "theming2.png")
agg_png(pngfile, width = 20, height = 36, units = "cm", res = 300)
plot(p_xq)
invisible(dev.off())
knitr::include_graphics(pngfile)





# testing ---------------------------------------------------------------------------


p_xq<- ggplot(df1) +
  geom_point(aes(x=x, y=y), 
             col = contour_cols[findInterval(dens, seq(0, max(dens), length.out = k))], 
             pch = 16) +
  geom_line(data=predictions_quantiles_m4_long, 
            aes(x=x, y=value, group=variable),
            color='#082133') +
  labs(title=expression(paste('GAMLSS | y ~ D(',mu,',',sigma,',',nu,',',tau,')')),
       # subtitle='Location, scale, and shape\nare modeled as functions of x',
       x = 'Age',
       y='Secret brain region')  + 
  annotate("text", x = 80, y = 55, size=5,
           label = paste(
             '\nSkewness =', '-1.5',
             '\nKurtosis =', '1'
           ))+
  ylim(c(min(y)*0.9,max(y)*1.1))+
  theme_minimal(); 
p_xq
