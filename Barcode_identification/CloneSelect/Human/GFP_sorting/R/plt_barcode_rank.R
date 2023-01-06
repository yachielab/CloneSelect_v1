library(ggplot2)
library(scales)
library(dplyr)
library(readxl)

df.rank <- read.csv("./barcode_rank_plot_07182022.csv", header = TRUE, sep = ",")
# Add human readable clone name
df.rank <- df.rank %>% mutate(bc_rank_name = sprintf("clone %03d", rank))

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

g.bc.rank <- ggplot(df.rank, aes(y = freq, x = rank)) +
    geom_point(data = df.rank %>% filter(is_isolated == "No"), aes(label = new_name, colour = is_isolated), alpha = 0.6) +
    geom_point(data = df.rank %>% filter(is_isolated == "Yes"), aes(label = new_name, colour = is_isolated)) +
    scale_colour_manual(values = c("Yes" = "deeppink2", "No" = "gray90")) +
    theme_classic() +
    scale_y_continuous(
        breaks = c(0, 0.01, 0.02, 0.03),
        limits = c(0, 0.03),
        expand = c(0, 0.001),
        label = scales::percent_format(suffix = "")
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 80)) +
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
ggsave(plot = g.bc.rank, filename = "barcode_rank_10072022.pdf", w = 3.0, h = 1.6)