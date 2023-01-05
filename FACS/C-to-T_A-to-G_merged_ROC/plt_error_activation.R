library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scales)
library(jcolors)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt in Keynote
DATE <- "04052022"

# loading data
df.mean <- read.csv("6x6_error_activation_mean_04132022.tsv", sep = ",", header = TRUE)

df.mean.CBE <- df.mean %>%
    dplyr::filter(System == "C-to-T")

df.mean.gRNA <- df.mean %>%
    dplyr::filter(System == "gRNA")

df.mean.CRISPRa <- df.mean %>%
    dplyr::filter(System == "CRISPRa")

# manual color map
Reds <- brewer.pal(n = 11, name = "PiYG")
Blues <- brewer.pal(n = 11, name = "PiYG")
Greens <- brewer.pal(n = 11, name = "RdBu")
Spectral <- brewer.pal(n = 6, name = "Spectral")

# plot ROC for each barcode
plt.mean.ROC <- function(df, xlim = c(0, 1.0), ylim = c(0, 1.0), colmap) {
    g <- ggplot(df, aes(y = mean.activation, x = mean.error)) +
        geom_abline(intercept = 0, linetype = "dashed", colour = "black", size = LINE.W) +
        geom_line(aes(colour = col_str), size = 0.7) +
        # scale_color_jcolors(palette = "pal5") +
        scale_colour_manual(values = colmap) +
        # scale_colour_brewer(palette = "Set1") +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        scale_y_continuous(
            labels = percent_format(accuracy = 1, suffix = ""),
            breaks = c(0, 0.1, 0.2, 0.3, 0.4),
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
            strip.background = element_rect(colour = NA, fill = "gray90"),
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
max.x <- 0.01
max.y <- 0.4
g.CBE <- plt.mean.ROC(df.mean.CBE, xlim = c(0, max.x), ylim = c(0, max.y), Spectral) + labs(x = "", y = "")
g.gRNA <- plt.mean.ROC(df.mean.gRNA, xlim = c(0, max.x), ylim = c(0, max.y), Spectral) + labs(x = "", y = "")
g.CRISPRa <- plt.mean.ROC(df.mean.CRISPRa, xlim = c(0, max.x), ylim = c(0, max.y), Spectral) + labs(x = "", y = "")

w <- 1.5 # 1.0 # 1.5
h <- 1.9 # 2.0 # 1.9
ggsave(plot = g.CBE, w = w, h = h, filename = "plots/CBE_error_activation_inset_plot_03142022_v2.pdf")
ggsave(plot = g.gRNA, w = w, h = h, filename = "plots/gRNA_error_activation_inset_plot_03142022.pdf")
ggsave(plot = g.CRISPRa, w = w, h = h, filename = "plots/CRIPSRa_error_activation_inset_plot_03142022.pdf")

# For supplementary figure
g.CBE <- plt.mean.ROC(df.mean.CBE, xlim = c(0, 1.0), ylim = c(0, 1.0), Spectral) + labs(x = "", y = "") + theme(aspect.ratio = 1) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), label = scales::percent_format(suffix = "")) + scale_y_continuous(expand = c(0, 0), limits = c(0, 100), label = scales::percent_format(suffix = ""))
g.gRNA <- plt.mean.ROC(df.mean.gRNA, xlim = c(0, 1.0), ylim = c(0, 1.0), Spectral) + labs(x = "", y = "") + theme(aspect.ratio = 1) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), label = scales::percent_format(suffix = "")) + scale_y_continuous(expand = c(0, 0), limits = c(0, 100), label = scales::percent_format(suffix = ""))
g.CRISPRa <- plt.mean.ROC(df.mean.CRISPRa, xlim = c(0, 1.0), ylim = c(0, 1.0), Spectral) + labs(x = "", y = "") + theme(aspect.ratio = 1) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100), label = scales::percent_format(suffix = "")) + scale_y_continuous(expand = c(0, 0), limits = c(0, 100), label = scales::percent_format(suffix = ""))

ggsave(plot = g.CBE, w = 2.0, h = 2.0, filename = "plots/CBE_error_activation_1x_plot_03142022_v2.pdf")
ggsave(plot = g.gRNA, w = 2.0, h = 2.0, filename = "plots/gRNA_error_activation_1x_plot_03142022.pdf")
ggsave(plot = g.CRISPRa, w = 2.0, h = 2.0, filename = "plots/CRIPSRa_error_activation_1x_plot_03142022.pdf")