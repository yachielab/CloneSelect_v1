library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

stat.df <- read.csv("./processed_csv/GFP_stat_ABE3x3_04112022.csv", header = TRUE)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

# calculate fold-change between OT and NT
fc.df <- stat.df %>%
    mutate(percent = replace(percent, percent == 0, 0.001)) %>%
    group_by(system, tx_reagent, rep, is_matched) %>%
    summarise(mean.GFP = mean(percent)) %>%
    pivot_wider(values_from = "mean.GFP", names_from = "is_matched") %>%
    mutate(fc = OT / NT, log2fc = log2(fc))

# Summarized jitter plot
fc.summary <- fc.df %>%
    group_by(system) %>%
    summarise(NT_mean = median(NT), OT_mean = median(OT), fc_mean = median(fc), log2fc = median(log2fc))

fc.df$System <- factor(fc.df$system, levels = c("A-to-G_BC", "CRISPRa_BC", "gRNA_BC"))
g.jitter <- ggplot(fc.df, aes(x = system, y = log2fc)) +
    geom_jitter(aes(colour = system), alpha = 0.7) +
    geom_crossbar(
        data = fc.summary, aes(ymin = log2fc, ymax = log2fc),
        size = 0.4, col = "black", width = 0.7
    ) +
    theme_classic() +
    theme(
        panel.spacing.x = unit(0.6, "mm"),
        panel.background = element_rect(colour = "transparent", fill = "transparent"),
        plot.background = element_rect(colour = "transparent", fill = "transparent"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust = 1, size = FONT.SIZE),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE, margin = margin(0, 1, 0, 0)),
        axis.title = element_text(size = FONT.SIZE),
        axis.line = element_line(size = LINE.W, colour = "black"),
        axis.ticks = element_line(size = LINE.W, colour = "black"),
        axis.ticks.length = unit(1, "mm"),
    ) +
    scale_color_manual(values = c("A-to-G_BC" = "#c5921b", "gRNA_BC" = "#7EBC41", "CRISPRa_BC" = "#4293C3")) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 20)) +
    scale_x_discrete(labels = c("A-to-G_BC" = "A-to-G\nreporter", "gRNA_BC" = "High-copy\nCRISPRa", "CRISPRa_BC" = "Low-copy\nCRISPRa")) +
    guides(colour = "none") +
    labs(y = "Log2 fold-change of EGFP positive cells", x = "")

ggsave(plot = g.jitter, filename = "plots/fc_jitter_3x2_04142022.pdf", device = cairo_pdf, w = 1.2, h = 2.4)