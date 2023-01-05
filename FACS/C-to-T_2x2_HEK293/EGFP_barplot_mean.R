library(ggplot2)
library(dplyr)
library(scales)

df <- read.table("11072021_CA053_summary.tsv", sep = "\t", header = TRUE)

FONT.SIZE <- 8
LABEL.FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

plt.mean.bar <- function(x, system, xlim) {
    g1 <- ggplot(df %>% filter(System == system), aes(y = Percent_FITC, x = Query_gRNA)) +
        facet_grid(~Barcode) +
        geom_bar(aes(fill = gRNA_type), position = "dodge", stat = "summary", fun = "mean") +
        geom_point(colour = "black", size = 2.0, shape = 1) +
        theme_classic() +
        scale_x_discrete(expand = c(0, 0), labels = c("V4-BC2" = "BC-C1", "V4-BC4" = "BC-C2")) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, xlim)) +
        labs(x = "", y = "") +
        scale_fill_manual(values = c("gray", "green3")) +
        theme(
            panel.background = element_rect(colour = "transparent", fill = "transparent"),
            plot.background = element_rect(colour = "transparent", fill = "transparent"),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.ticks = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(size = LINE.W),
            axis.title = element_text(size = FONT.SIZE),
            axis.text.y = element_text(colour = "black", size = LABEL.FONT.SIZE),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = LABEL.FONT.SIZE, colour = "black")
        ) +
        guides(fill = "none")
    return(g1)
}

g.CtoT <- plt.mean.bar(df, "CtoT_BC", 100)
g.CRISPRa <- plt.mean.bar(df, "CRISPRa_BC", 100)
g.gRNA <- plt.mean.bar(df, "gRNA_BC", 100)

fig.w <- 1
fig.h <- 2
ggsave(plot = g.CtoT, filename = "./plt_mean_EGFP/CtoT_mean_EGFP_11032022.pdf", w = fig.w, h = fig.h, device = cairo_pdf)
ggsave(plot = g.CRISPRa, filename = "./plt_mean_EGFP/CRISPRa_mean_EGFP_11032022.pdf", w = fig.w, h = fig.h, device = cairo_pdf)
ggsave(plot = g.gRNA, filename = "./plt_mean_EGFP/gRNA_mean_EGFP_11032022.pdf", w = fig.w, h = fig.h, device = cairo_pdf)