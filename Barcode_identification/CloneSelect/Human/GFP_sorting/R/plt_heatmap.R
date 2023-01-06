library(ggplot2)
library(scales)
library(dplyr)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

df.all <- read.csv("./norm_barcode_rep1_R_friendly_10072022.csv", header = TRUE, sep = ",")
df.AUG <- read.csv("./norm_AUG_barcode_rep1_R_friendly.csv", header = TRUE, sep = ",")

df.rank <- read.csv("./barcode_rank_plot_07182022.csv", header = TRUE, sep = ",")
df.rank <- df.rank %>% mutate(bc_rank_name = sprintf("clone %03d", rank))

# merge for all data (GTG and AUG)
df.all.merged <- merge(df.all, df.rank, by.x = "bc", by.y = "old_name") %>%
    arrange(rank) %>%
    mutate(new_order = 1:length(order))

# merge for AUG data
df.AUG.merged <- merge(df.AUG, df.rank, by.x = "bc", by.y = "old_name") %>%
    arrange(rank) %>%
    mutate(new_order = 1:length(order))


# plot heatmap for all barcode enrichment
g.heatmap.all <- ggplot(df.all.merged, aes(y = reorder(query_bc, -order), x = reorder(bc_rank_name, order))) +
    geom_tile(aes(fill = norm_bc), colour = NA) +
    theme_classic() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90), # element_blank(),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.line = element_blank(),
        panel.border = element_blank()
    ) +
    scale_fill_gradientn(
        colours = c("gray98", "dodgerblue", "blue", "blue2", "darkblue"),
        limits = c(0, 100), na.value = "white"
    ) +
    labs(x = "", y = "") +
    guides(fill = "none")

ggsave(plot = g.heatmap.all, filename = "BC_heatmap_10072022.pdf", device = cairo_pdf, w = 5.2, h = 2.1)

# plot heatmap for AUG edited barcode enrichment
g.heatmap.AUG <- ggplot(df.AUG.merged, aes(y = reorder(query_bc, -order), x = reorder(bc_rank_name, order))) +
    geom_tile(aes(fill = norm_bc), colour = NA) +
    theme_classic() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90), # element_blank(),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.line = element_blank(),
        panel.border = element_blank()
    ) +
    scale_fill_gradientn(
        colours = c("gray97", "green", "green2", "darkgreen"),
        limits = c(0, 100), na.value = "white"
    ) +
    labs(x = "", y = "") #+guides(fill = "none")
ggsave(plot = g.heatmap.AUG, filename = "BC_heatmap_AUG_09252022.pdf", device = cairo_pdf, w = 5.2, h = 2.1)

stop("done for plotting heatmaps")

# for all barcode types
g.jitter <- ggplot(df.all, aes(x = reorder(query_bc, order), y = norm_bc / 100)) +
    geom_jitter(data = df %>% filter(is_match == "NT"), aes(colour = is_match), width = 0.3, size = 0.7, alpha = 0.5, shape = 16) +
    geom_jitter(data = df %>% filter(is_match == "OT"), aes(colour = is_match), width = 0.1, size = 1.5, shape = 16) +
    scale_y_continuous(limits = c(-0.01, 1), expand = c(0.015, 0), position = "right") +
    scale_colour_manual(values = c("OT" = "blue2", "NT" = "gray60")) +
    theme_classic() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE, angle = 90),
        axis.line = element_line(colour = "black", size = LINE.W),
        axis.ticks = element_line(colour = "black", size = LINE.W)
    ) +
    labs(x = "", y = "") +
    guides(colour = "none")
ggsave(plot = g.jitter, filename = "BC_jitter_07192022_.pdf", device = cairo_pdf, w = 2.20, h = 1.2)

# for AUG edited barcodes
g.jitter.AUG <- ggplot(df.AUG, aes(x = reorder(query_bc, order), y = norm_bc / 100)) +
    geom_jitter(data = df %>% filter(is_match == "NT"), aes(colour = is_match), width = 0.3, size = 0.7, alpha = 0.5, shape = 16) +
    geom_jitter(data = df %>% filter(is_match == "OT"), aes(colour = is_match), width = 0.1, size = 1.5, shape = 16) +
    scale_y_continuous(limits = c(-0.01, 1), expand = c(0.015, 0), position = "right") +
    scale_colour_manual(values = c("OT" = "green3", "NT" = "gray60")) +
    theme_classic() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE, angle = 90),
        axis.line = element_line(colour = "black", size = LINE.W),
        axis.ticks = element_line(colour = "black", size = LINE.W)
    ) +
    labs(x = "", y = "") +
    guides(colour = "none")
ggsave(plot = g.jitter.AUG, filename = "BC_jitter_AUG_09252022_.pdf", device = cairo_pdf, w = 2.20, h = 1.2)