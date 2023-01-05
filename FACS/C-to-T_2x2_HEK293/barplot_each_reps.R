library(ggplot2)
library(scales)
library(dplyr)

FONT.SIZE.TITLE <- 9
FONT.SIZE.LEGEND <- 7
LINE.W <- 0.3481 # equivalent to 0.75pt in Keynote

df <- read.csv('11072021_CA053_summary.tsv', sep='\t', header=TRUE)
df$System <- factor(df$System, levels=c('CtoT_BC','CRISPRa_BC','gRNA_BC'))
df$gRNA_type <- factor(df$gRNA_type, levels=c('OT_gRNA', 'NT_gRNA'))

plt.barplot <- function(d, y.lim.max, ...) {
    g <- ggplot(d, aes(y=Percent_FITC, x=Replicate)) + 
        facet_wrap(~gRNA_type, scales='fixed') +
        geom_bar(fill='green3', stat='identity', width=0.5) +
        theme_classic() +
        scale_y_continuous(expand=c(0, 0), limits=c(0, y.lim.max)) + 
        labs(x='', y='%EGFP\npositive cells') + 
        theme(
           axis.ticks.y=element_line(colour='black', size=LINE.W),
            axis.ticks.x=element_line(colour=NA),
            axis.line=element_line(size=LINE.W, colour='black'),
            axis.title=element_text(size=FONT.SIZE.TITLE), 
    	   strip.background=element_blank(),
    	   axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, colour='black', size=FONT.SIZE.LEGEND),
    	   axis.text.y=element_text(colour='black', size=FONT.SIZE.LEGEND),
    	   legend.position = "none")
    return(g)
}

for (qbc in c('V4-BC2', 'V4-BC4')) {
    g1 <- plt.barplot(df%>% filter(Barcode==qbc, System=='CtoT_BC'), 20)
    ggsave(plot=g1, filename=sprintf('CA035_CtoT_BC_%s.pdf', qbc), w=3.5, h=1.5)
    g2 <- plt.barplot(df%>% filter(Barcode==qbc, System=='CRISPRa_BC'), 100)
    ggsave(plot=g2, filename=sprintf('CA035_CRISPRa_BC_%s.pdf', qbc), w=3.5, h=1.5)
    g3 <- plt.barplot(df%>% filter(Barcode==qbc, System=='gRNA_BC'), 8)
    ggsave(plot=g3, filename=sprintf('CA035_gRNA_BC_%s.pdf', qbc), w=3.5, h=1.5)
}




#ggsave(filename=sprintf('CA053_barplot_%s_%s.pdf', target.BC, gRNA.type), w=3.5, h=1.5)




