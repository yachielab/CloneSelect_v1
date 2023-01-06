library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales)
library(readxl)
library(reshape2)
library(stringr)

FONT.SIZE <- 7
LINE.W <- 0.3481 # equivalent to 0.75pt in Keynote

breaks_log10 <- function(x) {
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    10^(seq.int(low, high))
}

breaks_5log10 <- function(x) {
    low <- floor(log10(min(x) / 5))
    high <- ceiling(log10(max(x) / 5))
    5 * 10^(seq.int(low, high))
}

# Choose which one you plot
target <- "BC1580" # "BC1580" # or "BC100"
# target <- "BC100"

if (target == "BC1580") {
    df <- read_excel("./data/ECS_BC1550pool_heatmap_data_08052021.xlsx", skip = 4)
    df$BC_name <- factor(df$BC_name, levels = c("BC1", "BC2", "BC3", "BC4")) # no longer required
} else if (target == "BC100") {
    df <- read_excel("./data/ECS_BC100pool_heatmap_data_08052021.xlsx", skip = 4)
    df$BC_name <- factor(df$BC_name, levels = c("BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "BC7", "BC8", "BC9", "BC10")) # no longer required
}

df_zeo <- df %>% filter(plate == "+Zeo")
df_ctrl <- df %>% filter(plate == "-Zeo")

g.BC_freq.bar <- ggplot(df_zeo, aes(x = clone_id, y = BC_freq * 100)) +
    facet_grid(~ reorder(clone_id, -BC_freq), scale = "free") +
    # geom_bar(width = 0.4, position = position_dodge(width = 0.4), fill = "#4d9221", stat = "identity") +
    # We can't use geom_bar since data contains negative value when it log10 transformed
    geom_segment(aes(x = clone_id, xend = clone_id, y = 1e-4, yend = BC_freq * 100), size = 5, colour = "#4d9221") +
    theme_classic() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = NA),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 5),
        axis.line = element_line(size = 0.2),
        axis.title = element_text(size = FONT.SIZE),
        axis.ticks.length = unit(0.5, "mm"),
        panel.spacing.x = unit(0, "mm"),
        plot.margin = unit(c(t = 0.3, r = 0, b = -0.5, l = 0), "lines")
    ) +
    scale_y_log10(
        limits = c(1e-4, 3),
        expand = c(0, 0),
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(x = "", y = "") +
    guides(fill = "none")

g.zeo_freq.bar <- ggplot(df_zeo, aes(x = clone_id, BC_freq, y = Zeo_colony_percent)) +
    facet_grid(~ reorder(clone_id, -BC_freq), scale = "free_x") +
    geom_segment(aes(x = clone_id, xend = clone_id, y = 1e-4, yend = Zeo_colony_percent), size = 5, colour = "#c51b7d") +
    theme_classic() +
    scale_y_log10(
        limits = c(1e-4, 1),
        expand = c(0, 0),
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x))
    ) +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.2),
        axis.ticks.x = element_line(colour = NA),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 5),
        axis.line = element_line(size = 0.2),
        axis.title = element_text(size = FONT.SIZE),
        axis.ticks.length = unit(0.5, "mm"),
        panel.spacing.x = unit(0, "mm"),
        plot.margin = unit(c(t = 0.3, r = 0, b = -0.5, l = 0), "lines")
    ) +
    labs(x = "", y = "")

if (target == "BC100") {
    p1 <- ggarrange(g.BC_freq.bar, g.zeo_freq.bar, ncol = 1, nrow = 2, heights = c(1, 1), align = "v")
    message("Plotting BC100 barplot")
    ggsave(
        plot = p1, filename = "BC100_Zeocin_freq_bar_07182022.pdf",
        device = cairo_pdf, h = 0.8, w = 5.07
    )
} else if (target == "BC1580") {
    message("Plotting BC1580 barplot")
    p1 <- ggarrange(g.BC_freq.bar, g.zeo_freq.bar, ncol = 1, nrow = 2, heights = c(1, 1), align = "v")
    ggsave(
        plot = p1, filename = "BC1580_Zeocin_freq_bar_07182022.pdf",
        device = cairo_pdf, h = 0.8, w = 2.25
    )
}

# data for heatmap
df_zeo_melt <- df_zeo %>%
    select(id, no, clone_id, BC_name, BC_seq_detected, is_single_BC, is_CAA_conv, is_expected_BC, BC_freq) %>%
    melt(., id.vars = c("id", "clone_id", "no", "BC_name", "BC_seq_detected", "BC_freq")) %>%
    mutate(uniq_x_lab = paste0(clone_id, ":", no))

df_zeo_melt$variable <- factor(df_zeo_melt$variable, levels = c("is_CAA_conv", "is_expected_BC", "is_single_BC"))
df_zeo_melt$BC_name <- factor(df_zeo_melt$BC_name, levels = c("BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "BC7", "BC8", "BC9", "BC10"))

df_ctrl_melt <- df_ctrl %>%
    select(id, no, clone_id, BC_name, BC_seq_detected, is_single_BC, is_CAA_conv, is_expected_BC, BC_freq) %>%
    melt(., id.vars = c("id", "clone_id", "no", "BC_name", "BC_seq_detected", "BC_freq")) %>%
    mutate(uniq_x_lab = paste0(clone_id, ":", no)) %>%
    mutate(rank = str_replace(clone_id, "Clone ", "") %>% as.numeric())

df_ctrl_melt$variable <- factor(df_ctrl_melt$variable, levels = c("is_CAA_conv", "is_expected_BC", "is_single_BC"))
# df_ctrl_melt$BC_name <- factor(df_ctrl_melt$BC_name, levels = c("BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "BC7", "BC8", "BC9", "BC10"))

g.heat.zeo <- ggplot(df_zeo_melt, aes(x = reorder(uniq_x_lab, -BC_freq), y = variable)) +
    geom_tile(aes(fill = value), linejoin = "round", colour = "gray90", size = 0.5) +
    scale_fill_manual(values = c("N.D." = "gray20", "no" = "white", "yes" = "#d53e4f")) +
    theme_classic() +
    theme(
        plot.background = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        # plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "lines")
    ) +
    coord_fixed() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "", y = "") +
    guides(fill = "none")

if (target == "BC100") {
    ggsave(
        plot = g.heat.zeo, filename = "BC100_Zeocin_heatmap05102022.pdf",
        device = cairo_pdf, h = 0.8, w = 5.17
    )
} else if (target == "BC1580") {
    ggsave(
        plot = g.heat.zeo, filename = "BC1580_Zeocin_heatmap05102022.pdf",
        device = cairo_pdf, h = 0.8, w = 2.3
    )
}

g.heat.ctrl <- ggplot(df_ctrl_melt, aes(x = reorder(uniq_x_lab, rank), y = variable)) +
    geom_tile(aes(fill = value), linejoin = "round", colour = "gray90", size = 0.5) +
    scale_fill_manual(values = c("N.D." = "gray20", "no" = "white", "yes" = "#d53e4f")) +
    theme_classic() +
    theme(
        plot.background = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = FONT.SIZE, colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        # plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "lines")
    ) +
    coord_fixed() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "", y = "") +
    guides(fill = "none")

if (target == "BC100") {
    message("Plotting control BC1580 heatmap")
    ggsave(plot = g.heat.ctrl, filename = "BC100_Zeocin_ctrl_heatmap05102022.pdf", device = cairo_pdf, h = 3.17, w = 5.2)
} else if (target == "BC1580") {
    message("Plotting control BC1580 heatmap")
    ggsave(plot = g.heat.ctrl, filename = "BC1580_Zeocin_ctrl_heatmap05102022.pdf", device = cairo_pdf, h = 3.17, w = 2.3)
}