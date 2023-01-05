library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scales)
# library(patchwork)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt in Keynote
DATE <- "04052022"

# loading data
df.mean <- read.csv("processed_csv/ABE_error_activation_mean_03072022.csv", sep = ",", header = TRUE)
df.mean.ABE <- df.mean %>%
    dplyr::filter(system == "A-to-G_BC")

df.mean.gRNA <- df.mean %>%
    dplyr::filter(system == "gRNA_BC")

df.mean.CRISPRa <- df.mean %>%
    dplyr::filter(system == "CRISPRa_BC")

plt.mean.ROC <- function(df, xlim = c(0, 1.0), ylim = c(0, 1.0), cols) {
    g <- ggplot(df, aes(y = mean.activation, x = mean.error)) +
        geom_abline(intercept = 0, linetype = "dashed", colour = "black", size = LINE.W) +
        geom_line(aes(colour = cell_line), size = 0.7) +
        scale_colour_manual(values = cols) +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        scale_y_continuous(
            labels = percent_format(accuracy = 1, suffix = ""),
            # breaks = c(0, 0.1, 0.2, 0.3, 0.4),
            expand = c(0, 0)
        ) +
        scale_x_continuous(
            labels = percent_format(accuracy = 1, suffix = ""),
            # breaks = function(x) pretty(x, n = 3),
            expand = c(0.005, 0)
        ) +
        labs(x = "%EGFP error", y = "%EGFP activation") +
        theme_classic() +
        theme(
            # aspect.ratio = 1,
            panel.background = element_rect(fill = "transparent", colour = "transparent"),
            plot.background = element_rect(fill = "transparent", , colour = "transparent"),
            plot.title = element_text(colour = "black", size = FONT.SIZE),
            strip.background = element_rect(colour = NA, fill = "gray90"),
            panel.border = element_blank(),
            axis.title = element_text(colour = "black", size = FONT.SIZE),
            axis.ticks = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(colour = "black", size = LINE.W),
            axis.text = element_text(colour = "black", size = FONT.SIZE),
        ) +
        guides(colour = "none")
    return(g)
}

plt.zoom <- function(gg) {
    # xlim: (0, 0.05)
    # ylim: (0, 0.5)
    g1 <- gg +
        scale_y_continuous(
            expand = c(0, 0), breaks = c(0, 0.5),
            labels = percent_format(accuracy = 1, suffix = "")
        ) +
        scale_x_continuous(
            expand = c(0.015, 0), breaks = c(0, 0.05),
            labels = percent_format(accuracy = 1, suffix = "")
        ) +
        theme(panel.background = element_rect(fill = "transparent")) +
        labs(x = "", y = "") + guides(colour = "none")
    return(g1)
}

# manual color map
Pals <- brewer.pal(n = 3, name = "Dark2")

xmax <- 0.01
ymax <- 0.4
g.ABE <- plt.mean.ROC(df.mean.ABE, xlim = c(0, xmax), ylim = c(0, ymax), Pals) + labs(x = "", y = "")
g.gRNA <- plt.mean.ROC(df.mean.gRNA, xlim = c(0, xmax), ylim = c(0, ymax), Pals) + labs(x = "", y = "")
g.CRISPRa <- plt.mean.ROC(df.mean.CRISPRa, xlim = c(0, xmax), ylim = c(0, ymax), Pals) + labs(x = "", y = "")

# Main figure with zoom-in figure
w <- 1.5
h <- 1.9
ggsave(plot = g.ABE, w = w, h = h, filename = "plots/ABE_error_activation_inset_plot_03142022_v2.pdf")
ggsave(plot = g.gRNA, w = w, h = h, filename = "plots/gRNA_error_activation_inset_plot_03142022.pdf")
ggsave(plot = g.CRISPRa, w = w, h = h, filename = "plots/CRIPSRa_error_activation_inset_plot_03142022.pdf")


# For 100%-100% scale plot
# Suppelementary
g.ABE <- plt.mean.ROC(df.mean.ABE, xlim = c(0, 1.0), ylim = c(0, 1.0), Pals) + labs(x = "", y = "")
g.gRNA <- plt.mean.ROC(df.mean.gRNA, xlim = c(0, 1.0), ylim = c(0, 1.0), Pals) + labs(x = "", y = "")
g.CRISPRa <- plt.mean.ROC(df.mean.CRISPRa, xlim = c(0, 1.0), ylim = c(0, 1.0), Pals) + labs(x = "", y = "")

w <- 1.95
h <- 1.95
ggsave(plot = g.ABE, w = w, h = h, filename = "plots/ABE_error_activation_plot_07182022_v2.pdf")
ggsave(plot = g.gRNA, w = w, h = h, filename = "plots/gRNA_error_activation_plot_07182022.pdf")
ggsave(plot = g.CRISPRa, w = w, h = h, filename = "plots/CRIPSRa_error_activation_plot_07182022.pdf")