library(ggplot2)
library(dplyr)
library(scales)
library(readxl)

FONT.SIZE <- 7
LINE.W <- 0.232

lib <- "BC100"

if (lib == "BC100") {
  df <- read_excel("data/BC100_common_barcodes_05252022.xlsx")
  df <- df %>% mutate(rank = 1:length(barcode))
  min.barcode.count <- 10 # ignore minor barcode sequences
} else if (lib == "BC1580") {
  min.barcode.count <- 1 # this is because read depth isn't enough
  df <- read_excel("data/BC1580_common_barcodes_05252022.xlsx")
  df <- df %>% mutate(rank = 1:length(barcode))
}

g1 <- ggplot(df, aes(x = rank, y = mean_freq)) +
  geom_point(
    data = df %>% filter(knee_flag == "yes", mean_raw_count > min.barcode.count, is.na(is_isolated)),
    colour = "gray50", alpha = 0.6, shape = 16, size = 1.2
  ) +
  geom_point(
    data = df %>% filter(knee_flag == "no", mean_raw_count > min.barcode.count, is.na(is_isolated)),
    colour = "gray90", alpha = 0.6, shape = 16, size = 1.2
  ) +
  geom_point(
    data = df %>% filter(mean_raw_count > min.barcode.count, is_isolated == "Yes"),
    colour = "red2", shape = 16, size = 1.2
  ) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 0.05), labels = scales::percent_format(suffix = "")) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = LINE.W),
    axis.line = element_line(size = LINE.W, colour = "black"),
    axis.title = element_text(size = FONT.SIZE),
    axis.text = element_text(colour = "black", size = FONT.SIZE)
  ) +
  labs(x = "Barcode rank", y = "Clone frequency (%)") +
  guides(colour = "none")

if (lib == "BC100") {
  ggsave(plot = g1, file = "zeo_barcode_rank_BC100.pdf", w = 1.9, h = 1.2)
} else if (lib == "BC1580") {
  ggsave(plot = g1 + scale_y_continuous(breaks = c(0, 0.01), limits = c(0, 0.01), labels = scales::percent_format(suffix = "", accuracy = 1)), file = "zeo_barcode_rank_BC1580.pdf", w = 1.9, h = 1.2)
}