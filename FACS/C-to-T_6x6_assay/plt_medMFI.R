library(ggplot2)
library(scales)
library(dplyr)

stat.df <- read.csv("GFP_stat_CA062_04112022.csv", header = TRUE)
stat.df <- stat.df %>% dplyr::filter(is_match == "OT")

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

# mean EGFP MFI values for crossbar
fc.summary <- stat.df %>%
    group_by(System, Cell_line2) %>%
    summarise(B525.FITC.A = mean(B525.FITC.A))

g1 <- ggplot(stat.df, aes(y = B525.FITC.A, x = Cell_line2)) +
    facet_wrap(~System) +
    geom_point(
        aes(group = System, fill = Cell_line2),
        size = 2, shape = 21,
        position = position_dodge(width = 0.9)
    ) +
    scale_y_log10(
        breaks = c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        expand = c(0, 0),
        limits = c(1, 10^6)
    ) +
    scale_fill_brewer(palette = "Spectral") +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, size = FONT.SIZE),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.title = element_text(size = 8),
        axis.line = element_line(size = LINE.W, colour = "black"),
        axis.ticks.y = element_line(size = LINE.W, colour = "black"),
        axis.ticks.x = element_blank()
    ) +
    labs(x = "", y = "Median fluorescence intensity (a.u.)") #+guides(fill = "none")

ggsave(plot = g1, filename = "6x6_MFI_plt_09242022.pdf", w = 2.3, h = 2.8, device = cairo_pdf)




stop("ok")



plt.mean.bar <- function(df) {
    g <- ggplot(df, aes(y = B525.FITC.A, x = Cell_line2)) +
        geom_jitter(aes(colour = System, group = System), alpha = 0.3, size = 2) +
        geom_crossbar(
            data = fc.summary, aes(ymin = B525.FITC.A, ymax = B525.FITC.A, colour = System),
            size = 0.4, width = 0.7
        ) +
        scale_color_manual(values = c("C-to-T" = "#C51B7D", "gRNA" = "#7EBC41", "CRISPRa" = "#4293C3")) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0.6, "mm"),
            panel.background = element_rect(colour = "transparent", fill = "transparent"),
            plot.background = element_rect(colour = "transparent", fill = "transparent"),
            strip.background = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust = 1, size = FONT.SIZE),
            axis.text.y = element_text(colour = "black", size = FONT.SIZE, margin = margin(0, 1, 0, 0)),
            axis.title = element_text(size = FONT.SIZE),
            axis.line = element_line(size = LINE.W, colour = "black"),
            axis.ticks = element_line(size = LINE.W, colour = "black"),
            axis.ticks.length = unit(1, "mm"),
        ) +
        # scale_x_discrete(expand = c(0, 0)) +
        scale_y_log10(
            breaks = c(10^4, 10^5, 10^6),
            labels = scales::trans_format("log10", math_format(10^.x)),
            limits = c(10^4, 10^6)
        ) +
        labs(y = "Median EGFP fluorescence (a.u.)", x = "") +
        guides(fill = "none", colour = "none")
    return(g)
}

fig.w <- 1.9
fig.h <- 3.0
g.MFI <- plt.mean.bar(stat.df)
ggsave(plot = g.MFI, filename = "6x6_MFI_plt_04192022.pdf", w = fig.w, h = fig.h, device = cairo_pdf)