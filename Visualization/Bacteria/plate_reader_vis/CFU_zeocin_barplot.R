library(ggplot2)

## Page 125's data (Zeocin CFU)
df <- data.frame(
	cond=c('None', 'LacO+RBS', 'LacO'),
	colony=c(447, 283, 481),
	CFU=c(4.47, 2.83, 4.81))

df$cond <- factor(df$cond, levels=c('None', 'LacO+RBS', 'LacO'))
g1 <- ggplot(df, aes(x=cond, y=CFU)) +
    geom_bar(fill='gold2', stat='identity', position='dodge', width=0.7) +
    theme_classic() +
    theme(
    	axis.ticks.x=element_blank(),
    	axis.ticks=element_line(colour='black')) +
    scale_y_continuous(expand=c(0, 0), limit=c(0, 6)) +
    labs(x='', y='')

ggsave(plot=g1, file='Zeo_CUF.pdf', w=2, h=3)
