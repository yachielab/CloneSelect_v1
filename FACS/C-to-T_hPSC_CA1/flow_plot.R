library(ggplot2)

FONT.SIZE <- 7
FONT.SIZE.LABEL <- 8
LINE.W <- 0.232

df <- read.csv("perGFP+v2.csv", header = TRUE)
df$condition <- factor(df$condition, c("OT gRNA", "NT gRNA"))
g <- ggplot(df, aes(x = condition, y = GFP_percent)) +
    geom_bar(aes(group = as.factor(condition), fill = condition), position = "dodge", stat = "summary", fun = "mean") +
    geom_point(
        colour = "black", size = 2,
        shape = 1
    ) +
    scale_fill_manual(values = c("NT gRNA" = "gray", "OT gRNA" = "green3")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +
    theme(
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
    labs(x = "", y = "%EGFP positive cells") +
    guides(fill = "none")

ggsave(plot = g, filename = "10052022_2FACS_Omar_barplot.pdf", w = 0.8, h = 2.5, device = cairo_pdf)