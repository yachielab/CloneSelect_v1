library(ggplot2)
library(scales)
library(dplyr)

LABEL.FONT.SIZE <- 7
LEGENG.FONT.SIZE <- 7
LINE.W <- 0.232 # equivalent to 0.75pt in Keynote

df <- read.csv("700K_colony_count.csv", header = TRUE, sep = ",")

breaks_log10 <- function(x) {
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    10^(seq.int(low, high))
}

g <- ggplot(df, aes(x = Rep, y = Count)) +
    geom_bar(stat = "identity", position = "dodge", fill = "dodgerblue3") +
    theme_classic() +
    labs(x = "", y = "Estimated library complexity") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_log10(
        breaks = breaks_log10,
        expand = c(0, 0),
        limits = c(1, 1 * 10^7),
        labels = trans_format("log10", math_format(10^.x))
    ) +
    theme_classic() +
    theme(
        # plot.margin = margin(1, 1, 1, 1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = LEGENG.FONT.SIZE),
        axis.text = element_text(colour = "black", size = LABEL.FONT.SIZE),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = LABEL.FONT.SIZE)
    )
ggsave(plot = g, filename = "700K_colony_count_09162021.pdf", w = 0.9, h = 2.7)