library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)
library(tidyr)

INDUCER_ORDERS <- c("No inducer", "Ara", "IPTG", "Ara+IPTG")

df <- read.table("araBADp-EGFP_T7p-nCas9_DataFrame.tsv", header = TRUE, sep = "\t")
df$inducer_type <- factor(df$inducer_type, levels = INDUCER_ORDERS)
df$gRNA <- factor(df$gRNA, levels = c("OT_gRNA", "NT_gRNA"))

df2 <- df %>%
    group_by(inducer_type, gRNA, reporter_type) %>%
    summarise(
        mean.EGFP = mean(norm_EGFP),
        sd = sd(norm_EGFP, na.rm = TRUE)
    )

# Re-order factors
df2$inducer_type <- factor(df2$inducer_type, levels = INDUCER_ORDERS)
df2$gRNA <- factor(df2$gRNA, levels = c("OT_gRNA", "NT_gRNA"))


# get fold-change stats
pivot_wider(df2 %>% filter(reporter_type == "A-to-G") %>% select(-sd), names_from = "gRNA", values_from = "mean.EGFP")




FONT.SIZE <- 7
LINE.W <- 0.232

QUERY_TYPE <- "A-to-G" # A-to-G_ctrl

g1 <- ggplot(df %>% filter(reporter_type == QUERY_TYPE), aes(y = norm_EGFP / 1000, x = gRNA)) +
    facet_grid(~inducer_type) +
    theme_classic() +
    geom_col(data = df2 %>% filter(reporter_type == QUERY_TYPE), aes(y = mean.EGFP / 1000, fill = gRNA)) +
    geom_point(colour = "black", size = 1.5, shape = 1) +
    theme(
        axis.ticks.y = element_line(colour = "black", size = LINE.W),
        axis.ticks.x = element_line(colour = NA),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, size = FONT.SIZE),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.line = element_line(size = LINE.W, colour = "black"),
        axis.title = element_text(size = FONT.SIZE),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(0.8, "mm"),
        plot.background = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(colour = "transparent", fill = "transparent")
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 125), labels = scales::comma_format(suffix = "K")) +
    scale_fill_manual(values = c("OT_gRNA" = "green3", "NT_gRNA" = "gray85")) +
    labs(x = "", y = "Normalized EGFP intensity (a.u.)") +
    guides(fill = "none")

ggsave(plot = g1, file = sprintf("plots/%s_EGFP_intensity_v2.pdf", QUERY_TYPE), w = 1.6, h = 2.8)

QUERY_TYPE <- "A-to-G_ctrl"
g2 <- ggplot(df %>% filter(reporter_type == QUERY_TYPE), aes(y = norm_EGFP / 1000, x = gRNA)) +
    facet_grid(~inducer_type) +
    theme_classic() +
    geom_col(data = df2 %>% filter(reporter_type == QUERY_TYPE), aes(y = mean.EGFP / 1000, fill = gRNA)) +
    geom_point(colour = "black", size = 1.5, shape = 1) +
    theme(
        axis.ticks.y = element_line(colour = "black", size = LINE.W),
        axis.ticks.x = element_line(colour = NA),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, size = FONT.SIZE),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.line = element_line(size = LINE.W, colour = "black"),
        axis.title = element_text(size = FONT.SIZE),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(0.8, "mm"),
        plot.background = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(colour = "transparent", fill = "transparent")
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 260), labels = scales::comma_format(suffix = "K")) +
    scale_fill_manual(values = c("OT_gRNA" = "green3", "NT_gRNA" = "gray85")) +
    labs(x = "", y = "Normalized EGFP intensity (a.u.)") +
    guides(fill = "none")

ggsave(plot = g2, file = sprintf("plots/%s_EGFP_intensity_v2.pdf", QUERY_TYPE), w = 1.6, h = 2.8)