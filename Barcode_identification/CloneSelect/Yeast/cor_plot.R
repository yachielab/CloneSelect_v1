library(ggplot2)
library(dplyr)
library(scales)
library(readxl)
library(stringr)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

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

plt.dist <- function(df, x, y, point_color) {
    g <- ggplot(df, aes_string(x = x, y = y)) +
        geom_abline(intercept = 0, colour = "black", linetype = "dashed", size = 0.2) +
        geom_point(data = df %>% filter(is.na(Used_in_isolation) == TRUE), colour = point_color, alpha = 0.7, size = 0.7, shape = 16) +
        geom_point(data = df %>% filter(Used_in_isolation == "Yes"), colour = "deeppink2", alpha = 1.0, size = 1.0, shape = 16) +
        # geom_point(aes(colour = knee_flag), alpha = 0.8, size = 0.7, shape = 16) +
        labs(x = str_replace(x, "^r", "R"), y = str_replace(y, "^r", "R")) +
        theme_classic() +
        scale_x_log10(
            breaks = breaks_log10,
            limits = c(1, 1 * 10^5),
            labels = trans_format("log10", math_format(10^.x))
        ) +
        scale_y_log10(
            breaks = breaks_log10,
            limits = c(1, 1 * 10^5),
            labels = trans_format("log10", math_format(10^.x))
        ) +
        theme(
            strip.background = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            axis.ticks = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(size = LINE.W, colour = "black"),
            axis.title = element_text(size = FONT.SIZE),
            axis.text = element_text(colour = "black", size = FONT.SIZE),
        ) +
        coord_fixed(ratio = 1)
    return(g)
}

# figure size
w <- 1.5
h <- 1.5

# BC100
# df_100 <- read.csv("data/BC100_normalized_count_06012022.csv", header = TRUE)
df_100 <- read_excel("../YeastCS_pool100_BCids_positives_09142022.xlsx")
ggsave(plot = plt.dist(df_100, "rep1", "rep2", "#65c13a"), filename = "./plots/BC100_R1R2_12212022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = plt.dist(df_100, "rep2", "rep3", "#65c13a"), filename = "./plots/BC100_R2R3_12212022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = plt.dist(df_100, "rep3", "rep1", "#65c13a"), filename = "./plots/BC100_R1R3_12212022.pdf", w = w, h = h, device = cairo_pdf)

# BC1580
# df_1580 <- read.csv("data/BC1580_normalized_count_06012022.csv", header = TRUE)
df_1580 <- read_excel("../YeastCS_pool1580_BCids_positives_09142022.xlsx")
ggsave(plot = plt.dist(df_1580, "rep1", "rep2", "#3aafc1"), filename = "./plots/BC1580_R1R2_12212022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = plt.dist(df_1580, "rep2", "rep3", "#3aafc1"), filename = "./plots/BC1580_R2R3_12212022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = plt.dist(df_1580, "rep3", "rep1", "#3aafc1"), filename = "./plots/BC1580_R1R3_12212022.pdf", w = w, h = h, device = cairo_pdf)



cor(df_100$rep1, df_100$rep2, method = "pearson")
cor(df_100$rep2, df_100$rep3, method = "pearson")
cor(df_100$rep1, df_100$rep3, method = "pearson")