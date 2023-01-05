library(ggplot2)
library(scales)
library(dplyr)

stat.df <- read.csv("GFP_stat_CA062_04112022.csv", header = TRUE)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

plt.mean.bar <- function(df) {
    g <- ggplot(df, aes(y = percent, x = Query_gRNA2)) +
        facet_wrap(~Cell_line2) +
        geom_bar(aes(group = Query_gRNA2, fill = is_match), position = "dodge", stat = "summary", fun = "mean") +
        geom_point(colour = "black", size = 1.5, shape = 1) +
        scale_fill_manual(values = c("NT" = "gray", "OT" = "green3")) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0.6, "mm"),
            panel.background = element_rect(colour = "transparent", fill = "transparent"),
            plot.background = element_rect(colour = "transparent", fill = "transparent"),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust = 1, size = FONT.SIZE),
            axis.text.y = element_text(colour = "black", size = FONT.SIZE, margin = margin(0, 1, 0, 0)),
            axis.title = element_text(size = FONT.SIZE),
            axis.line = element_line(size = LINE.W, colour = "black"),
            axis.ticks = element_line(size = LINE.W, colour = "black"),
            axis.ticks.length = unit(1, "mm"),
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 90), breaks = c(0, 30, 60, 90)) +
        labs(y = "", x = "") +
        guides(fill = "none")
    return(g)
}

fig.w <- 1.3
fig.h <- 2.25

# Group1 => BC2, BC4, BC6 assay space
# Group2 => BC8, BC10, BC12 assay space
for (sys in c("C-to-T", "CRISPRa", "gRNA")) {
    for (gr in c("Group1", "Group2")) {
        g <- stat.df %>%
            dplyr::filter(System == sys, Assay_group == gr) %>%
            plt.mean.bar()
        ggsave(plot = g, filename = sprintf("./plots/%s_%s_11042022.pdf", sys, gr), h = fig.h, w = fig.w, device = cairo_pdf)
    }
}