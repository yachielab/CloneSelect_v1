library(ggplot2)
library(dplyr)
library(scales)

DATE <- "04052022"
FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt in Keynote

# BC rename map
BC_renamed <- c("T3" = "BC-A1", "T4" = "BC-A2", "T5" = "BC-A3")

plt.mean.bar <- function(df) {
    #
    # For each bar plot, data are grouped by barcoded cell line
    # tx_reagenet: either query gRNA or reporter transfected
    #

    g <- ggplot(df, aes(y = percent, x = tx_reagent)) +
        facet_wrap(~cell_line) +
        geom_bar(aes(group = is_matched, fill = is_matched), position = "dodge", stat = "summary", fun = "mean") +
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
            axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, size = FONT.SIZE),
            axis.text.y = element_text(colour = "black", size = FONT.SIZE, margin = margin(0, 1, 0, 0)),
            axis.title = element_text(size = FONT.SIZE),
            axis.line = element_line(size = LINE.W, colour = "black"),
            axis.ticks = element_line(size = LINE.W, colour = "black"),
        ) +
        scale_x_discrete(expand = c(0, 0), labels = BC_renamed) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 90), breaks = c(0, 30, 60, 90)) +
        labs(y = "", x = "") +
        guides(fill = "none")
    return(g)
}

stat.df <- read.csv("processed_csv/CA059_percent_GFP_stat.csv", header = TRUE)

p.ABE <- stat.df %>%
    dplyr::filter(system == "A-to-G_BC") %>%
    # mutate(cell_line2 = paste0(cell_line, "_cell")) %>% # debug
    plt.mean.bar()

p.gRNA <- stat.df %>%
    dplyr::filter(system == "gRNA_BC") %>%
    plt.mean.bar()

p.CRISPRa <- stat.df %>%
    dplyr::filter(system == "CRISPRa_BC") %>%
    plt.mean.bar()

fig.w <- 1.3
fig.h <- 2.25
ggsave(
    plot = p.ABE, filename = sprintf("plots/ABE_GFP_barplot_%s.pdf", DATE),
    w = fig.w, h = fig.h, device = cairo_pdf
)
ggsave(
    plot = p.gRNA, filename = sprintf("plots/gRNA_GFP_barplot_%s.pdf", DATE),
    w = fig.w, h = fig.h, device = cairo_pdf
)
ggsave(
    plot = p.CRISPRa, filename = sprintf("plots/CRISPRa_GFP_barplot_%s.pdf", DATE),
    w = fig.w, h = fig.h, device = cairo_pdf
)