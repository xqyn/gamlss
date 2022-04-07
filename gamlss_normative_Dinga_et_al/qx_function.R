library(gamlss.dist)
library(mgcv)
library(reshape2)
library(ggplot2)
library(cowplot)
library(patchwork)
library(viridis)



params_to_quantiles_norm <- function(quantiles, params){
  as.data.frame(sapply(quantiles, 
                       function(q){
                         qnorm(p=q, mean=params[,1], sd = exp(params[,2]))
                       }))
}


params_to_quantiles_shash <- function(quantiles, params, qshash){
  as.data.frame(sapply(quantiles, 
                       function(q){
                         qshash(p=q, 
                                # param is called mu, but it expects 
                                # a vector of all 4 shash parameters
                                mu=params)
                       }))
}
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


xq_gamlss_stupid_func <- function(number_of_data, 
                                  intercept,
                                  B1,
                                  B2,
                                  skewness, 
                                  kurtosis, 
                                  density){
  n <- number_of_data
  x <- seq(0, 100, length.out = n)
  sigmas <- 0.5 + 1.5*seq(from = -1, to = 1, length.out = n)**2
  y  <- intercept + B1*x - B2*(x**2) + rSHASHo2(n = n, sigma = sigmas, nu=skewness, tau=kurtosis)
  df <- data.frame('x' = x, 'y' = y, 'dataset'='Training site')
  model <- gam(list(y ~ s(x), # fit mu as a smooth function of x
                    ~ s(x), # fit sigma as a smooth function of x
                    ~ 1, # fit nu (skewness) as an intercept
                    ~ 1), # fit tau (kurtosis) as an intercept
               family=shash(), # shash distribution instead of gaussian 
               data=df)
  predictions_params_model <- predict(model, newdata = df)
  quantiles <- pnorm(c(-2:2))
  qshash <- model$family$qf
  predictions_quantiles_model <- params_to_quantiles_shash(quantiles, 
                                                           predictions_params_model,
                                                           qshash)
  reshape_quantiles_to_long <- function(quantiles_df, x_var){
    quantiles_df$x <- x_var
    return(reshape2::melt(quantiles_df, id.vars = c('x')))}
  predictions_quantiles_model_long <- reshape_quantiles_to_long(predictions_quantiles_model, df$x)
  k <- density
  x<-df$x
  y<-df$y
  contour_cols <- viridis(k, alpha = 0.5)
  dens <- get_density(x, y, k)
  ## 2. Generating plots -----
  def.par<-par()
  p_xq<- ggplot(df) +
    geom_point(aes(x=x, y=y), 
               col = contour_cols[findInterval(dens, seq(0, max(dens), length.out = k))], 
               pch = 16, 
               alpha = .35,
               size=3) +
    geom_line(data=predictions_quantiles_model_long, 
              aes(x=x, y=value, group=variable),
              color='#082133') +
    labs(title=expression(paste('GAMLSS | y ~ D(',mu,',',sigma,',',nu,',',tau,')')),
         # subtitle='Location, scale, and shape\nare modeled as functions of x',
         x = 'Age',
         y='Secret brain region')  + 
    annotate("text", x = 80, y = max(y)+5, size=4,
             label = paste(
               '\nSkewness =', skewness,
               '\nKurtosis =', kurtosis
             ))+
    ylim(c(min(y)*0.9,max(y)*1.1))+
    theme_minimal(); 
}


#plot <- xq_gamlss_stupid_func(number_of_data = ,
#                                 intercept = ,
#                                 B1 = ,
#                                 B2 = ,
#                                 skewness = ,
#                                 kurtosis = ,
#                                 density = )

plot_xq <- xq_gamlss_stupid_func(number_of_data = 1000,
                                 intercept = 50,
                                 B1 = 0.15,
                                 B2 = 0.003,
                                 skewness = -1.5,
                                 kurtosis = 1,
                                 density = 100)
plot_xq


# Skewness --------------------------------------------------------------------------


p1 <- xq_gamlss_stupid_func(2000,50,0.15,0.003,-1.5,1,100)
p2 <- xq_gamlss_stupid_func(2000,50,0.15,0.003,-0.5,1,100)
p3 <- xq_gamlss_stupid_func(2000,50,0.15,0.003,1.5,1,100)
p4 <- xq_gamlss_stupid_func(2000,50,0.15,0.003,3,1,100)

p_all <- p1 + p2 + p3 + p4 + 
  plot_annotation(tag_levels = 'A') & 
  theme_cowplot() &
  theme(text=element_text(size=8),
        axis.text = element_text(size=9)); p_all


# Kurtosis --------------------------------------------------------------------------

p1 <- xq_gamlss_stupid_func(2000,50,0.15,0.003,-1,-1,100)
p2 <- xq_gamlss_stupid_func(2000,50,0.15,0.003,-1,0.5,100)
p3 <- xq_gamlss_stupid_func(2000,50,0.15,0.003,-1,1.5,100)
p4 <- xq_gamlss_stupid_func(2000,50,0.15,0.003,-1,2.5,100)

p_all <- p1 + p2 + p3 + p4 + 
  plot_annotation(tag_levels = 'A') & 
  theme_cowplot() &
  theme(text=element_text(size=8),
        axis.text = element_text(size=9)); p_all
