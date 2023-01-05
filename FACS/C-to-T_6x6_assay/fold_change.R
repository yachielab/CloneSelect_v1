library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

stat.df <- read.csv("GFP_stat_CA062_04112022.csv", header = TRUE)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

df.ca <- stat.df %>%
    filter(System == "CRISPRa") %>%
    select(c(Query_gRNA2, Cell_line2, Rep, is_match, percent))
df.c2t <- stat.df %>%
    filter(System == "C-to-T") %>%
    select(c(Query_gRNA2, Cell_line2, Rep, is_match, percent))

df.merged <- merge(df.c2t, df.ca,
    by = c("Query_gRNA2", "Cell_line2", "Rep", "is_match"),
    suffixes = c(".CtoT", ".CRISPRa")
) %>% mutate(
    percent.CtoT = replace(percent.CtoT, percent.CtoT == 0, 0.001),
    percent.CRISPRa = replace(percent.CRISPRa, percent.CRISPRa == 0, 0.001),
    fc = percent.CtoT / percent.CRISPRa
)

# calculate fold-change between OT and NT
fc.df <- stat.df %>%
    mutate(percent = replace(percent, percent == 0, 0.001)) %>%
    group_by(System, Query_gRNA2, Rep, is_match) %>%
    summarise(mean.GFP = mean(percent)) %>%
    pivot_wider(values_from = "mean.GFP", names_from = "is_match") %>%
    mutate(fc = OT / NT, log2fc = log2(fc))

# OT/NT barplot for each barcode
g.fc.barplt <- ggplot(fc.df, aes(y = log2fc, x = Query_gRNA2)) +
    facet_wrap(~System) +
    geom_bar(aes(group = Query_gRNA2, fill = Query_gRNA2), position = "dodge", stat = "summary", fun = "mean") +
    geom_point(colour = "black", size = 1.5, shape = 1) +
    theme_classic() +
    scale_fill_viridis_d(option = "F") +
    labs(y = "Log2 fold-change", x = "") +
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
    scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 15)) +
    scale_x_discrete(expand = c(0, 0)) +
    guides(fill = "none")

fig.w <- 2.0
fig.h <- 2.1
ggsave(plot = g.fc.barplt, filename = "plots/fc_barplot_04142022.pdf", device = cairo_pdf, w = fig.w, h = fig.h)

# Summarized jitter plot
fc.summary <- fc.df %>%
    group_by(System) %>%
    summarise(NT_mean = median(NT), OT_mean = median(OT), fc_mean = median(fc), log2fc = median(log2fc))

fc.df$System <- factor(fc.df$System, levels = c("C-to-T", "CRISPRa", "gRNA"))
g.jitter <- ggplot(fc.df, aes(x = System, y = log2fc)) +
    geom_jitter(aes(colour = System), alpha = 0.7) +
    geom_crossbar(
        data = fc.summary, aes(ymin = log2fc, ymax = log2fc),
        size = 0.4, col = "black", width = 0.7
    ) +
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
    scale_color_manual(values = c("C-to-T" = "#C51B7D", "gRNA" = "#7EBC41", "CRISPRa" = "#4293C3")) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 20)) +
    scale_x_discrete(labels = c("C-to-T" = "C-to-T\nreporter", "gRNA" = "High-copy\nCRISPRa", "CRISPRa" = "Low-copy\nCRISPRa")) +
    guides(colour = "none") +
    labs(y = "Log2 fold-change of EGFP positive cells", x = "")

ggsave(plot = g.jitter, filename = "plots/fc_jitter_04142022.pdf", device = cairo_pdf, w = 1.2, h = 2.4)