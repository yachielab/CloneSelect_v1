library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scales)
library(jcolors)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt in Keynote
DATE <- "04192022"

# loading data
df <- read.csv("./dose_error_activation_mean_04192022.csv", sep = ",", header = TRUE)

# manual color map
Reds <- brewer.pal(n = 11, name = "PiYG")
Blues <- brewer.pal(n = 7, name = "Blues")
Greens <- brewer.pal(n = 11, name = "RdBu")
Spectral <- brewer.pal(n = 7, name = "Spectral")

# plot ROC for each barcode
plt.mean.ROC <- function(df, xlim = c(0, 1.0), ylim = c(0, 1.0), colmap) {
    g <- ggplot(df, aes(y = activation, x = error)) +
        facet_wrap(~barcode) +
        geom_abline(intercept = 0, linetype = "dashed", colour = "black", size = LINE.W) +
        geom_line(aes(colour = dose_str), size = 0.7) +
        scale_colour_manual(values = colmap) +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        scale_y_continuous(
            labels = percent_format(accuracy = 1, suffix = ""),
            # breaks = c(0, 0.1, 0.2, 0.3, 0.4),
            expand = c(0, 0)
        ) +
        scale_x_continuous(
            labels = percent_format(accuracy = 0.1, suffix = ""),
            breaks = function(x) pretty(x, n = 3),
            expand = c(0.005, 0)
        ) +
        labs(x = "%EGFP error", y = "%EGFP activation") +
        theme_classic() +
        theme(
            # aspect.ratio = 1,
            panel.background = element_rect(fill = "transparent", colour = "transparent"),
            plot.background = element_rect(fill = "transparent", , colour = "transparent"),
            plot.title = element_text(colour = "black", size = FONT.SIZE),
            strip.background = element_rect(colour = "transparent", fill = "transparent"),
            panel.border = element_blank(),
            axis.title = element_text(colour = "black", size = FONT.SIZE),
            axis.ticks = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(colour = "black", size = LINE.W),
            axis.text = element_text(colour = "black", size = FONT.SIZE)
        ) +
        guides(colour = "none")
    return(g)
}

# For main figure - cropped
max.x <- 0.2
max.y <- 0.4
g.roc <- plt.mean.ROC(df, xlim = c(0, max.x), ylim = c(0, max.y), Blues) + labs(x = "", y = "")

w <- 4.0
h <- 1.9
# ggsave(plot = g.roc, w = w, h = h, filename = "plots/dose_error_activation_plot_04192022_v2.pdf")