library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(readxl)

FONT.SIZE <- 7
LINE.W <- 0.232

df <- read_excel("merged_96lib_pDNA_gDNA_01082023.xlsx")

df.rank <- df %>%
    arrange(-gDNA_RPM) %>%
    mutate(rank.g = 1:length(barcode))

g.bc.rank <- ggplot(df.rank, aes(y = gDNA_freq, x = rank.g)) +
    geom_point(data = df.rank %>% filter(Is_used_for_isolation == "No"), aes(colour = Is_used_for_isolation), alpha = 0.6) +
    geom_point(data = df.rank %>% filter(Is_used_for_isolation == "Yes"), aes(colour = Is_used_for_isolation)) +
    scale_colour_manual(values = c("Yes" = "deeppink2", "No" = "gray90")) +
    theme_classic() +
    scale_y_continuous(
        breaks = c(0, 0.01, 0.02, 0.03),
        limits = c(0, 0.03),
        expand = c(0, 0.001),
        label = scales::percent_format(suffix = "")
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 120)) +
    labs(x = "Barcode rank", y = "Clone frequency (%)") +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0.5, 3, 0.5, 0.5, "mm"),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE)
    ) +
    guides(color = "none")
ggsave(plot = g.bc.rank, filename = "plots/barcode_rank_03142022.pdf", w = 3.1, h = 1.6)