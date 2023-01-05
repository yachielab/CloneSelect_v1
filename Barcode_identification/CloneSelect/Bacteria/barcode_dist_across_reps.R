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
        geom_point(data = df %>% filter(is.na(is_isolated)), colour = point_color, alpha = 0.5, size = 0.7, shape = 16) +
        geom_point(data = df %>% filter(is_isolated == "Yes"), colour = "deeppink2", alpha = 1.0, size = 1.1, shape = 16) +
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
        coord_fixed(ratio = 1) +
        guides(colour = "none")
    return(g)
}

# barcode count data
# df_100 <- read_excel("data/BC100_common_barcodes_05252022.xlsx")
df_100 <- read_excel("ECS_BC100_processed.xlsx")
df_1580 <- read_excel("ECS_BC1580_processed.xlsx")

# figure size
w <- 1.5
h <- 1.5

# BC100
ggsave(plot = plt.dist(df_100, "rep1", "rep2", "#9DB354"), filename = "./plots/BC100_R1R2_06012022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = plt.dist(df_100, "rep2", "rep3", "#9DB354"), filename = "./plots/BC100_R2R3_06012022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = plt.dist(df_100, "rep3", "rep1", "#9DB354"), filename = "./plots/BC100_R1R3_06012022.pdf", w = w, h = h, device = cairo_pdf)

# BC1580
ggsave(plot = plt.dist(df_1580, "rep1", "rep2", "#6A54B3"), filename = "./plots/BC1580_R1R2_06012022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = plt.dist(df_1580, "rep2", "rep3", "#6A54B3"), filename = "./plots/BC1580_R2R3_06012022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = plt.dist(df_1580, "rep3", "rep1", "#6A54B3"), filename = "./plots/BC1580_R1R3_06012022.pdf", w = w, h = h, device = cairo_pdf)

# cor matrix

df_100 <- df_100 %>%
    mutate(rep1 = log10(rep1 + 1), rep2 = log10(rep2 + 1), rep3 = log10(rep3 + 1))
# cor(df_100[c(4, 5, 6)], method = "pearson")

df_1580 <- df_1580 %>%
    mutate(rep1 = log10(rep1 + 1), rep2 = log10(rep2 + 1), rep3 = log10(rep3 + 1))
# cor(df_1580[c(4, 5, 6)], method = "pearson")