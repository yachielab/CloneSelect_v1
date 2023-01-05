library(ggplot2)

df <- read.csv("./CA064_stat_10212022.csv", header = TRUE)

FONT.SIZE <- 7
FONT.SIZE.LABEL <- 8
LINE.W <- 0.232


df$gRNA_type <- factor(df$gRNA_type, levels = c("OT", "NT"))

g <- ggplot(df, aes(x = gRNA_type, y = percent)) +
    geom_bar(aes(group = gRNA_type, fill = gRNA_type), position = "dodge", stat = "summary", fun = "mean") +
    geom_point(
        colour = "black", size = 2,
        shape = 1
    ) +
    scale_fill_manual(values = c("OT" = "#ff4f4f", "NT" = "gray")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +
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
    labs(x = "", y = "%mCherry positive cells") +
    guides(fill = "none")
ggsave(plot = g, filename = "10212022_CA064_barplot.pdf", w = 0.8, h = 2.1, device = cairo_pdf)