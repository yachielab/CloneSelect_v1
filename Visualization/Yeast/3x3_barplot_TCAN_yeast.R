library(ggplot2)
library(scales)
library(dplyr)

FONT.SIZE <- 7
FONT.SIZE.LABEL <- 8
LINE.W <- 0.232

df <- read.table("data/3x3_TCAN.tsv", header = TRUE, sep = "\t")
df <- df %>% mutate(is_match = ifelse(gRNA == Cell, "OT", "NT"))

g <- ggplot(df, aes(x = gRNA, y = Int / 1000)) +
    facet_wrap(~Cell) +
    geom_bar(aes(group = as.factor(No), fill = is_match), position = "dodge", stat = "summary", fun = "mean") +
    scale_fill_manual(values = c("OT" = "#ff5c5c", "NT" = "gray")) +
    geom_point(colour = "black", size = 1.5, shape = 1) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1200), label = scales::comma_format(suffix = "K")) +
    scale_x_discrete(expand = c(0, 0), labels = c("BC2" = "yBC-C1", "BC4" = "yBC-C2", "BC8" = "yBC-C3")) +
    theme(
        panel.spacing.x = unit(0.8, "mm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE.LABEL),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.text.x = element_text(colour = "black", size = FONT.SIZE, angle = 90, vjust = 0.5, hjust = 1),
    ) +
    labs(x = "", y = "") +
    guides(fill = "none")

ggsave(plot = g, filename = "./plots/3x3_RFP_TCAN_08172022.pdf", w = 1.71, h = 1.9, device = cairo_pdf)