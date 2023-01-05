library(ggplot2)
library(dplyr)

df <- read.table("data_summary_01302020.tsv", header = TRUE, sep = "\t")

df <- df %>% mutate(gRNA_type = case_when(Cell_line == gRNA ~ "OT", Cell_line != gRNA ~ "NT"))

FONT.SIZE <- 7
LINE.W <- 0.232 # equivalent to 0.75pt in Keynote

g1 <- ggplot(df, aes(y = EGFP_int, x = gRNA)) +
    facet_grid(~Cell_line) +
    geom_bar(aes(fill = gRNA_type), position = "dodge", stat = "summary", fun = "mean") +
    geom_point(colour = "black", size = 2.0, shape = 1) +
    theme_classic() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +
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
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = FONT.SIZE, colour = "black")
    ) +
    guides(fill = "none")

ggsave(plot = g1, filename = "GFP_hela_barplt_11032022.pdf", h = 2, w = 0.93, device = cairo_pdf)