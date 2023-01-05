library(ggplot2)
library(dplyr)
library(scales)

df <- read.table("11072021_CA053_summary.tsv", sep = "\t", header = TRUE)
df$System <- factor(df$System, levels = c("CtoT_BC", "CRISPRa_BC", "gRNA_BC"))

FONT.SIZE <- 7
# LINE.W <- 0.3481 # equivalent to 0.75pt in Keynote
LINE.W <- 0.232 # 0.5pt

breaks_log10 <- function(x) {
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    10^(seq.int(low, high))
}

g1 <- ggplot(df %>% filter(gRNA_type == "OT_gRNA"), aes(x = Query_gRNA, y = Median_GFP_exp)) +
    facet_grid(~System) +
    geom_point(aes(fill = Barcode), size = 2, shape = 21, colour = "black") +
    scale_fill_manual(values = c("V4-BC2" = "#D53E4F", "V4-BC4" = "#FC8E59")) +
    scale_y_log10(
        breaks = c(10^3, 10^4, 10^5, 1 * 10^6),
        labels = scales::trans_format("log10", math_format(10^.x)),
        expand = c(0, 0),
        limits = c(1000, 1000000)
    ) +
    labs(x = "", y = "Median EGFP fluorescence intensity (a.u.)") +
    theme_classic() +
    theme(
        panel.background = element_rect(colour = "transparent", fill = "transparent"),
        plot.background = element_rect(colour = "transparent", fill = "transparent"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = LINE.W),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.text.x = element_blank()
    ) +
    guides(fill = "none")
ggsave(plot = g1, filename = "EGFP_median_MFI_09262022.pdf", w = 1.2, h = 1.95)