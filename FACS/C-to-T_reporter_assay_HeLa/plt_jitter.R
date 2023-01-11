library(dplyr)
library(ggplot2)
library(scales)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

df <- read.table("data_summary_01302020.tsv", header = TRUE, sep = "\t")
df <- df %>% mutate(gRNA_type = case_when(Cell_line == gRNA ~ "OT", Cell_line != gRNA ~ "NT"))

plt.jitter <- function(df) {
    g1 <- ggplot(df, aes(x = gRNA_type, y = EGFP_int)) +
        geom_jitter(aes(colour = gRNA_type), size = 0.5) +
        scale_colour_manual(values = c("OT" = "green3", "NT" = "gray")) +
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
        scale_y_continuous(limits = c(-1, 50), breaks = c(0, 25, 50)) +
        labs(y = "", x = "") +
        guides(colour = "none")
    return(g1)
}

ggsave(plot = plt.jitter(df), filename = "jitter_01102023.pdf", h = 0.9, w = 0.93, device = cairo_pdf)