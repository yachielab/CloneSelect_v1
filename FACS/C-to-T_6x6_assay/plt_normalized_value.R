library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

df <- read.csv("05182022_normalize_GFP.tsv", header = TRUE, sep = "\t")
df$system <- factor(df$system, levels = c("C-to-T", "A-to-G"))

FONT.SIZE <- 7
FONT.SIZE.TITLE <- 8
LINE.W <- 0.232 # 0.5pt

plt <- function(df, summary, ymax) {
    g <- ggplot(df, aes(x = system, y = norm.value)) +
        geom_point(aes(colour = system), size = 1.0, position = position_jitter(width = 0.3, seed = 42)) +
        geom_errorbar(
            data = summary, aes(ymin = norm.value, ymax = norm.value),
            size = 0.5, col = "black", width = 0.6
        ) +
        theme_classic() +
        labs(x = "", y = "Normalized\nreporter activation") +
        scale_colour_manual(values = c("C-to-T" = "#E05A99", "A-to-G" = "#87A5E0")) +
        theme(
            plot.background = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            axis.line = element_line(size = LINE.W, colour = "black"),
            axis.ticks = element_line(size = LINE.W, colour = "black"),
            axis.title = element_text(size = FONT.SIZE.TITLE, colour = "black"),
            axis.text = element_text(size = FONT.SIZE, colour = "black")
        ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
        scale_x_discrete(labels = c(
            "C-to-T" = "C-to-T\nreporter",
            "A-to-G" = "A-to-G\nreporter"
        )) +
        guides(colour = "none")
    return(g)
}

df.OT <- df %>% filter(is_match == "OT")
df.NT <- df %>% filter(is_match == "NT")

# Summarized jitter plot
summary.OT <- df.OT %>%
    group_by(system) %>%
    summarise(norm.value = median(norm.value))

summary.NT <- df.NT %>%
    group_by(system) %>%
    summarise(norm.value = median(norm.value))


# OT pair
ggsave(plot = plt(df.OT, summary.OT, 6.0), filename = "./plots/normalized_OT_08162022.pdf", w = 1.5, h = 1.5, device = cairo_pdf)

# NT pair
ggsave(plot = plt(df.NT, summary.NT, 6.0), filename = "./plots/normalized_NT_08162022.pdf", w = 1.5, h = 1.5, device = cairo_pdf)