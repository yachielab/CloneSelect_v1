library(dplyr)
library(ggplot2)
library(scales)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

df <- read.csv("GFP_stats.csv")

plt.jitter <- function(target.system) {
    g1 <- ggplot(df %>% filter(AID_type == target.system), aes(x = Target_type, y = percent)) +
        geom_jitter(aes(colour = Target_type), size = 0.5) +
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
        scale_y_continuous(limits = c(-1, 30), breaks = c(0, 30)) +
        labs(y = "", x = "") +
        guides(colour = "none")
    return(g1)
}

ggsave(plot = plt.jitter("nCas9-AID"), filename = "./plots/nCas9_AID_jitter_01102023.pdf", h = 0.9, w = 0.70, device = cairo_pdf)
ggsave(plot = plt.jitter("dCas9-AID"), filename = "./plots/dCas9_AID_jitter_01102023.pdf", h = 0.9, w = 0.70, device = cairo_pdf)
ggsave(plot = plt.jitter("nCas9"), filename = "./plots/nCas9_jitter_01102023.pdf", h = 0.9, w = 0.70, device = cairo_pdf)