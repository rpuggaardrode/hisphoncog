library(tidyverse)
library(viridis)

df <- read.csv('vot_voi.csv', sep=';')
dia_polygon <- read.csv('polygons_dialectareas.csv', sep=';')
glimpse(df)
df %>% group_by(stop) %>% summarize(count = n())

lats <- unique(df$lat)
longs <- unique(df$long)
locs <- data.frame(lat = lats, lon = longs)

wave_tg <- function(sound, textgrid, tier,
                    start=0, end=Inf) {

  tg <- rPraat::tg.read(textgrid)
  snd <- rPraat::snd.read(sound, from=start, to=end, units='seconds')

  sig <- snd$sig
  fs <- snd$fs
  t <- snd$t - start
  lab <- tg[[tier]]$label

  t1 <- tg[[tier]]$t1 - start
  t2 <- tg[[tier]]$t2 - start
  line_vec <- t1[-1]
  rn <- diff(range(sig))

  h <- floor(max(sig)*10)/10
  l <- ceiling(min(sig)*10)/10
  if ((h*10)-(l*10) != 1) {
    r <- (h*10)-(l*10)
  } else {
    l <- l-0.05
    r <- 3
  }

  plot(t, sig, type='l',
       ylim=c(min(sig)-rn*0.15, max(sig)),
       ylab='Amplitude', xlab='Time (s)',
       yaxp=c(h, l, r))
  abline(h=c(min(sig)-rn*0.2, min(sig)-rn*0.05))
  abline(v=line_vec, lty='dotted')
  text(x=t2-((t2-t1)/2),
       y=min(sig)-rn*0.11,
       labels=lab)

}

wave_tg('asp_annot.wav', 'asp_annot.TextGrid', tier='vot')
par(mfrow=c(2,2))
wave_tg('int_fv.wav', 'int_fv.TextGrid', tier='stop')
wave_tg('int_nfv.wav', 'int_nfv.TextGrid', tier='stop')
wave_tg('pp_fv.wav', 'pp_fv.TextGrid', tier='stop')
wave_tg('pp_nfv.wav', 'pp_nfv.TextGrid', tier='stop')
dev.off()

ggplot(locs) +
  aes(x=lon, y=lat) +
  geom_polygon(data=dia_polygon, aes(x=lat, y=long, group=dialect),
               fill=NA, color='grey', lwd=0.6) +
  coord_fixed(1.7) +
  geom_point(size=2.5, alpha=0.5, col='red') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "#FFFFFF", color = "#FFFFFF"),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks = element_blank(), legend.background = element_blank(),
        legend.title = element_blank(), legend.text = element_blank())

df %>% group_by(sg) %>% summarize(count = n())

bdg <- df %>% filter(sg == 'negsg')

bdg[which(bdg$voi==''),'voi'] <- 'vl'
bdg[which(bdg$boundary==''),'boundary'] <- 'n'

for(col in c('parish', 'gender', 'dialect', 'voi', 'pal',
             'height', 'bkn2l', 'rdness', 'stop', 'str', 'boundary')){
  bdg[[col]] <- as.factor(bdg[[col]])
}

contrast <- cbind(c(-0.5,+0.5))
colnames(contrast) <- '-nst+st'
contrasts(bdg$str) <- contrast

colnames(contrast) <- '-f+m'
contrasts(bdg$gender) <- contrast

colnames(contrast) <- '-nrd+rd'
contrasts(bdg$rdness) <- contrast

colnames(contrast) <- '-bk+nbk'
contrasts(bdg$bkn2l) <- contrast

colnames(contrast) <- '-nbd+bd'
contrasts(bdg$boundary) <- contrast

colnames(contrast) <- '-npal+pal'
contrasts(bdg$pal) <- contrast

contrast <- cbind(c(+1/3, +1/3, -2/3), c(+0.5, -0.5, 0))
colnames(contrast) <- c('+AB-V', '-A+B')
contrasts(bdg$stop) <- contrast

contrast = cbind(c(+2/3, -1/3, -1/3), c(0, -0.5, +0.5))
colnames(contrast) = c("-ml+h", "-l+m")
contrasts(bdg$height) = contrast

save(bdg, file='bdg.Rda')

gam_voi <- mgcv::bam(voi ~ stop + str + gender + rdness +
                       bkn2l + boundary + height + pal +
                       s(long,lat) +
                       s(parish, bs='re') +
                       s(parish, by=stop, bs='re') +
                       s(parish, by=str, bs='re') +
                       s(parish, by=rdness, bs='re') +
                       s(parish, by=bkn2l, bs='re') +
                       s(parish, by=boundary, bs='re') +
                       s(parish, by=height, bs='re') +
                       s(parish, by=pal, bs='re'),
                     data=bdg,
                     discrete=T,
                     family=binomial(link='logit'),
                     nthreads=10,
                     control=list(trace=T))

gam_voi_ml <- mgcv::bam(voi ~ stop + str + gender + rdness +
                       bkn2l + boundary + height +
                       s(parish, bs='re') +
                       s(parish, by=stop, bs='re') +
                       s(parish, by=str, bs='re') +
                       s(parish, by=rdness, bs='re') +
                       s(parish, by=bkn2l, bs='re') +
                       s(parish, by=boundary, bs='re') +
                       s(parish, by=height, bs='re') +
                       s(long,lat),
               data=bdg,
               method='ML',
               family=binomial(link='logit'),
               nthreads=10,
               control=list(trace=T))

gam_voi_null <- mgcv::bam(voi ~ stop + str + gender + rdness +
                       bkn2l + boundary + height +
                       s(parish, bs='re') +
                       s(parish, by=stop, bs='re') +
                       s(parish, by=str, bs='re') +
                       s(parish, by=rdness, bs='re') +
                       s(parish, by=bkn2l, bs='re') +
                       s(parish, by=boundary, bs='re') +
                       s(parish, by=height, bs='re'),
                     data=bdg,
                     # discrete=T,
                     method='ML',
                     family=binomial(link='logit'),
                     nthreads=10,
                     control=list(trace=T))

load('bdg_models.Rda')

#Read up on interpreting residuals for logistic GAMMs
itsadug::compareML(gam_voi, gam_voi_null, suggest.report=T)

check <- mgcv::gam.check(gam_vot)

gam_map_sig <- function(mod, sig_level=1) {
  viz <- mgcViz::getViz(mod)
  plot(mgcViz::sm(viz, 1)) + coord_fixed(1.7) +
    mgcViz::l_fitRaster(pTrans = function(.p) .p<sig_level) +
    scale_fill_gradientn(colours=c('blue','white','red'),
                         na.value='grey', name=NULL) +
    geom_polygon(data=dia_polygon, aes(x=lat, y=long, group=dialect),
                          fill=NA, color='black', inherit.aes=F) +
    ggtitle('') +
    theme(panel.border=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank(),
          text = element_text(size=18))
}




viz <- mgcViz::getViz(gam_vot)

viz <- mgcViz::getViz(gam_voi)
plot(mgcViz::sm(viz, 1)) + ggplot2::coord_fixed(1.7) +
  mgcViz::l_fitRaster(pTrans = function(.p) .p<1) +
  #mgcViz::l_fitRaster() +
  scale_fill_gradientn(colours=c("blue","white","red"), na.value='grey', name=NULL) +
  ggplot2::geom_polygon(data=dia_polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                        fill=NA, color='black', inherit.aes=F) +
  #ggtitle('Fitted log likelihood of voicing\n(white when p>0.05)') +
  ggtitle('') +
  theme(panel.border=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        text = element_text(size=18))


summary(gam_voi, re.test=F)

#####

ptk <- df %>% filter(sg == 'possg')

for(col in c('parish', 'gender', 'dialect', 'voi', 'pal',
             'height', 'bkn2l', 'rdness', 'stop', 'str', 'boundary')){
  ptk[[col]] <- as.factor(ptk[[col]])
}

contrast <- cbind(c(-0.5,+0.5))
colnames(contrast) <- '-nst+st'
contrasts(ptk$str) <- contrast

colnames(contrast) <- '-f+m'
contrasts(ptk$gender) <- contrast

colnames(contrast) <- '-npal+pal'
contrasts(ptk$pal) <- contrast

colnames(contrast) <- '-nrd+rd'
contrasts(ptk$rdness) <- contrast

colnames(contrast) <- '-bk+nbk'
contrasts(ptk$bkn2l) <- contrast

colnames(contrast) <- '-nbd+bd'
contrasts(ptk$boundary) <- contrast

contrast = cbind(c(+2/3, -1/3, -1/3), c(0, -0.5, +0.5))
colnames(contrast) = c("-ml+h", "-l+m")
contrasts(ptk$height) = contrast

contrast <- cbind(c(+2/3, -1/3, -1/3), c(0, -0.5, +0.5))
colnames(contrast) <- c('-AB+V', '-B+A')
contrasts(ptk$stop) <- contrast

save(ptk, file='ptk.Rda')

gam_vot <- mgcv::bam(dur ~ stop + str + gender +
                       pal + bkn2l + rdness + height +
                       boundary +
                       s(long,lat) +
                       s(parish, bs='re') +
                       s(parish, by=stop, bs='re') +
                       s(parish, by=str, bs='re') +
                       s(parish, by=pal, bs='re') +
                       s(parish, by=bkn2l, bs='re') +
                       s(parish, by=rdness, bs='re') +
                       s(parish, by=height, bs='re') +
                       s(parish, by=boundary, bs='re'),
                     data=ptk,
                     discrete=T,
                     family='scat',
                     nthreads=10,
                     control=list(trace=T))

gam_vot_null <- mgcv::bam(dur ~ stop + str + gender +
                       pal + bkn2l + rdness + height +
                       s(parish, bs='re') +
                       s(parish, by=stop, bs='re') +
                       s(parish, by=str, bs='re') +
                       s(parish, by=pal, bs='re') +
                       s(parish, by=bkn2l, bs='re') +
                       s(parish, by=rdness, bs='re') +
                       s(parish, by=height, bs='re') +
                       s(parish, by=boundary, bs='re'),
                     data=ptk_dat,
                     discrete=T,
                     family='scat',
                     control=list(trace=T))

load('ptk_models.Rda')
load('Z:/jys/ptk_mod.Rda')

itsadug::compareML(gam_vot, gam_vot_null, suggest.report=T)

viz <- mgcViz::getViz(gam_vot)
plot(mgcViz::sm(viz, 1)) + ggplot2::coord_fixed(1.7) +
  mgcViz::l_fitRaster(pTrans = function(.p) .p<0.05) +
  #mgcViz::l_fitRaster() +
  scale_fill_gradientn(colours=c("red","white","blue"), na.value='grey', name=NULL,
                       breaks=c(6, 3, 0, -3, -6)) +
  ggplot2::geom_polygon(data=dia_polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                        fill=NA, color='black', inherit.aes=F) +
  # ggtitle('Fitted VOT\n(white when p>0.05)') +
  ggtitle('') +
  theme(panel.border=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        text = element_text(size=18))

summary(gam_voi, re.test=F)
summary(gam_vot, re.test=F)

###

plot.df <- mgcv::plot.gam(gam_vot, select=1, rug=F)
plot.df <- lapply(plot.df, function(x) {
  x$raw <- NULL
  x
})
plot.df <- plot.df[[1]]
plot_data <- expand.grid(plot.df$x, plot.df$y)
plot_data$fit <- as.vector(plot.df$fit)
plot_data$ci_lower <- plot_data$fit - qnorm(1-0.05/2)*as.vector(plot.df$se)
plot_data$ci_upper <- plot_data$fit + qnorm(1-0.05/2)*as.vector(plot.df$se)
legend_limits <- c(min(plot_data$ci_lower, na.rm = T), max(plot_data$ci_upper, na.rm = T))

gg_fit <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "fit")) +
  geom_tile() +
  ggtitle("Fitted log likelihood\n of voicing") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(limits = legend_limits,
                       low = 'red', mid = "white", high = 'blue',
                       name = "estimate") +
  coord_fixed(1.7) +
  geom_polygon(data=dia_polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                                  fill=NA, color='black', inherit.aes=F) +
  theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank())

gg_lower <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "ci_lower")) + geom_tile() +
  ggtitle("Lower CI boundary") + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(limits = legend_limits, low = 'red', mid = "white", high = 'blue', name = "estimate") +
  coord_fixed(1.7) + geom_polygon(data=dia_polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                                  fill=NA, color='black', inherit.aes=F) + theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank())

gg_upper <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "ci_upper")) + geom_tile() +
  ggtitle("Upper CI boundary") + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(limits = legend_limits, low = 'red', mid = "white", high = 'blue', name = "estimate") +
  coord_fixed(1.7) + geom_polygon(data=dia_polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                                  fill=NA, color='black', inherit.aes=F) + theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank())

library(patchwork)
combined <- gg_lower + gg_fit + gg_upper & theme(legend.position='right')
combined + plot_layout(guides='collect')

###

gam_map_extract_fit <- function(mod, select=1, ci=0, sig_level=1) {
  png('tmp')
  plot.df <- mgcv::plot.gam(mod, select=1, rug=F)
  dev.off()
  unlink('tmp', recursive=T)
  plot.df <- lapply(plot.df, function(x) {
    x$raw <- NULL
    x
  })
  plot.df <- plot.df[[select]]
  plot_data <- expand.grid(plot.df$x, plot.df$y)
  plot_data$fit <- as.vector(plot.df$fit)
  plot_data$se <- as.vector(plot.df$se)
  plot_data$p <- 1 - pnorm(abs(plot_data$fit)/plot_data$se)
  plot_data$alpha <- plot_data$p < sig_level

  cip <- 1-ci/100
  plot_data$ci_lower <- plot_data$fit - qnorm(1-cip/2)*as.vector(plot.df$se)
  plot_data$ci_upper <- plot_data$fit + qnorm(1-cip/2)*as.vector(plot.df$se)

  return(plot_data)
}

gam_map <- function(mod, poly_obj, select=1, sig_level=1, fill='fit',
                    legend_name='estimate', ...) {
  p <- gam_map_extract_fit(mod, select, sig_level=sig_level)
  ggplot(p, aes_string(x='Var1', y='Var2', fill = fill)) +
    geom_tile(alpha = p$alpha) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue',
                         na.value='grey',
                         name = legend_name, ...) +
    coord_fixed(1.7) +
    geom_polygon(data=poly_obj, aes(x=lat, y=long, group=dialect),
                 fill=NA, color='black', inherit.aes=F) + theme_bw() +
    theme(panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank())
}

gam_map_ci <- function(mod, poly_obj, select=1, ci=95) {
  p <- gam_map_extract_fit(mod, select, ci=ci)
  legend_lim <- c(min(p$ci_lower, na.rm = T), max(p$ci_upper, na.rm = T))

  fit <- gam_map(mod, poly_obj, select, limits=legend_lim) +
    ggtitle('Fit') +
    theme(plot.title = element_text(hjust = 0.5))
  lower <- gam_map(mod, poly_obj, select, limits=legend_lim, fill='ci_lower') +
    ggtitle(paste0('Lower ', ci, '%\nCI boundary')) +
    theme(plot.title = element_text(hjust = 0.5))
  upper <- gam_map(mod, poly_obj, select, limits=legend_lim, fill='ci_upper') +
    ggtitle(paste0('Upper ', ci, '%\nCI boundary')) +
    theme(plot.title = element_text(hjust = 0.5))

  combined <- lower + fit + upper & theme(legend.position='right')
  return(combined + patchwork::plot_layout(guides='collect'))
}

pTrans <- function(.p) {.p<0.05}
pTrans(1 - pnorm(abs(plot_data$fit)/as.vector(plot.df$se)))

gam_map_with_cis

plot.df <- mgcv::plot.gam(gam_vot, select=1, rug=F)
plot.df <- lapply(plot.df, function(x) {
  x$raw <- NULL
  x
})
plot.df <- plot.df[[1]]
plot_data <- expand.grid(plot.df$x, plot.df$y)
plot_data$fit <- as.vector(plot.df$fit)
plot_data$ci_lower <- plot_data$fit - qnorm(1-0.05/2)*as.vector(plot.df$se)
plot_data$ci_upper <- plot_data$fit + qnorm(1-0.05/2)*as.vector(plot.df$se)
legend_limits <- c(min(plot_data$ci_lower, na.rm = T), max(plot_data$ci_upper, na.rm = T))

gg_fit <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "fit")) + geom_tile() +
  ggtitle("Fitted VOT") + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(limits = legend_limits, low = 'red', mid = "white", high = 'blue', name = "estimate") +
  coord_fixed(1.7) + geom_polygon(data=dia_polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                                  fill=NA, color='black', inherit.aes=F) + theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank())

gg_lower <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "ci_lower")) + geom_tile() +
  ggtitle("Lower CI boundary") + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(limits = legend_limits, low = 'red', mid = "white", high = 'blue', name = "estimate") +
  coord_fixed(1.7) + geom_polygon(data=dia_polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                                  fill=NA, color='black', inherit.aes=F) + theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank())

gg_upper <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "ci_upper")) + geom_tile() +
  ggtitle("Upper CI boundary") + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  scale_fill_gradient2(limits = legend_limits, low = 'red', mid = "white", high = 'blue', name = "estimate") +
  coord_fixed(1.7) + geom_polygon(data=dia_polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                                  fill=NA, color='black', inherit.aes=F) + theme_bw() +
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank())

library(patchwork)
combined <- gg_lower + gg_fit + gg_upper & theme(legend.position='right')
combined + plot_layout(guides='collect')

#####


gam_map_with_cis <- function(gam_obj, select, polygon, fit_title) {
  plot.df <- mgcv::plot.gam(gam_obj, select=1, rug=F)
  plot.df <- lapply(plot.df, function(x) {
    x$raw <- NULL
    x
  })
  plot.df <- plot.df[[select]]
  plot_data <- expand.grid(plot.df$x, plot.df$y)
  plot_data$fit <- as.vector(plot.df$fit)
  plot_data$ci_lower <- plot_data$fit - qnorm(1-0.05/2)*as.vector(plot.df$se)
  plot_data$ci_upper <- plot_data$fit + qnorm(1-0.05/2)*as.vector(plot.df$se)
  legend_limits <- c(min(plot_data$ci_lower, na.rm = T),
                     max(plot_data$ci_upper, na.rm = T))

  gg_fit <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "fit")) +
    geom_tile() +
    ggtitle(fit_title) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_gradient2(limits = legend_limits,
                         low = 'red', mid = "white", high = 'blue',
                         name = "estimate") +
    coord_fixed(1.7) +
    geom_polygon(data=polygon,
                 ggplot2::aes(x=lat, y=long, group=dialect),
                                    fill=NA, color='black', inherit.aes=F) +
    theme_bw() +
    theme(panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank())

  gg_lower <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "ci_lower")) +
    geom_tile() +
    ggtitle("Lower CI boundary") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_gradient2(limits = legend_limits,
                         low = 'red', mid = 'white', high = 'blue',
                         name = 'estimate') +
    coord_fixed(1.7) +
    geom_polygon(data=polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                                    fill=NA, color='black', inherit.aes=F) +
    theme_bw() +
    theme(panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank())

  gg_upper <- ggplot(plot_data, aes_string(x="Var1", y="Var2", fill = "ci_upper")) +
    geom_tile() +
    ggtitle("Upper CI boundary") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_gradient2(limits = legend_limits,
                         low = 'red', mid = "white", high = 'blue',
                         name = "estimate") +
    coord_fixed(1.7) +
    geom_polygon(data=polygon, ggplot2::aes(x=lat, y=long, group=dialect),
                                    fill=NA, color='black', inherit.aes=F) +
    theme_bw() +
    theme(panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank())

  combined <- gg_lower + gg_fit + gg_upper & theme(legend.position='right')
  return(combined + patchwork::plot_layout(guides='collect'))
}

###

gam_map_with_cis(gam_voi, 16, dia_polygon, 'Fitted log likelihood\n of voicing')
gam_map_with_cis(gam_vot, 1, dia_polygon, 'Fitted VOT')

###

sum <- gam_voi %>%
  summary(re.test=F) %>%
  `[[`('p.table') %>%
  as.data.frame %>%
  mutate(OR = ifelse(Estimate < 0,
                     paste0('1:', round(exp(abs(Estimate)), 2)),
                     paste0(round(exp(abs(Estimate)), 2), ':1')))
  #%>%
  #extract(p.table)#%>%
  #extract2(p.table) #%>%
  #as.data.frame() %>%
  #mutate(OR = exp(abs(Estimate)))
