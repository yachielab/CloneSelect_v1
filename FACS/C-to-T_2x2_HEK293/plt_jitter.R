library(ggplot2)
library(scales)
library(dplyr)
library(stringr)

df <- read.table("11072021_CA053_summary.tsv", sep = "\t", header = TRUE)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

plt.jitter <- function(target.system) {
    g1 <- ggplot(df %>% filter(System == target.system), aes(x = gRNA_type, y = Percent_FITC)) +
        geom_jitter(aes(colour = gRNA_type), size = 0.5) +
        scale_colour_manual(values = c("OT_gRNA" = "green3", "NT_gRNA" = "gray")) +
        theme_classic() +
        theme(
            panel.background = element_rect(colour = "transparent", fill = "transparent"),
            plot.background = element_rect(colour = "transparent", fill = "transparent"),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust = 1, size = FONT.SIZE),
            axis.text.y = element_text(colour = "black", size = FONT.SIZE),
            axis.title = element_text(size = FONT.SIZE, colour = "black"),
            axis.line = element_line(size = LINE.W, colour = "black"),
            axis.ticks = element_line(size = LINE.W, colour = "black"),
        ) +
        # scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(limits = c(-1, 100), breaks = c(0, 50, 100)) +
        labs(y = "", x = "") +
        guides(colour = "none")
    return(g1)
}
ggsave(plot = plt.jitter("CtoT_BC"), filename = "./plots/C-to-T_jitter_01102023.pdf", h = 1.2, w = 0.99, device = cairo_pdf)

ggsave(plot = plt.jitter("CRISPRa_BC"), filename = "./plots/CRISPRa_jitter_01102023.pdf", h = 1.2, w = 0.99, device = cairo_pdf)

ggsave(plot = plt.jitter("gRNA_BC"), filename = "./plots/gRNA_jitter_01102023.pdf", h = 1.2, w = 0.99, device = cairo_pdf)