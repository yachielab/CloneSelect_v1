library(ggplot2)
library(dplyr)

FONT.SIZE <- 7
LINE.W <- 0.232

df <- read.table("data/CAN1_edit_ratio.tsv", header = TRUE)
g <- ggplot(df %>% filter(Condition == "Can(-)"), aes(y = Edit_ratio, x = Enzyme)) +
    geom_bar(aes(group = Enzyme, fill = Condition), position = "dodge", stat = "summary", fun = "mean") +
    geom_point(colour = "black", size = 1) +
    theme_classic() +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W, colour = "black"),
        axis.title = element_text(size = FONT.SIZE),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black", size = FONT.SIZE),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        legend.position = "none"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
    scale_fill_manual(values = c("dodgerblue")) +
    labs(x = "", y = "Editing frequency (%)") +
    guides(fill = "black")

ggsave(plot = g, filename = "./plots/CAN1_assay_edit_freq_05042022.pdf", w = 1.25, h = 2.9, device = cairo_pdf)