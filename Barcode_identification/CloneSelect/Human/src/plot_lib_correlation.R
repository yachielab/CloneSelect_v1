library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(stringr)
library(readxl)

FONT.SIZE <- 7
LINE.W <- 0.232
df <- read.csv("merged_96lib_pDNA_gDNA_01082023.csv", sep = ",", header = TRUE)

df.filt <- df %>%
    filter(!str_detect(clone_name, "only"))

plt.scatter <- ggplot(df.filt, aes(x = gRPM + 1, y = pRPM + 1)) +
    geom_abline(slope = 1, colour = "gray", linetype = "dashed") +
    geom_point(colour = "deeppink2", alpha = 0.5) +
    scale_x_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)),
        expand = c(0, 0), limits = c(1, 1e5)
    ) +
    scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)),
        expand = c(0, 0), limits = c(1, 1e5)
    ) +
    theme_classic() +
    theme(
        aspect.ratio = 1,
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE),
    ) +
    labs(x = "Barcode count in gDNA", y = "Barcode count in pDNA")
ggsave(plot = plt.scatter, filename = "plots/scatter_gDNA-pDNA_01082023.pdf", w = 1.7, h = 1.7)

# cor.test(log10(df$pRPM + 1), log10(df$gRPM + 1), method = "pearson")

df.pDNA <- df %>%
    select(bc, pRPM) %>%
    rename(RPM = pRPM) %>%
    arrange(-RPM) %>%
    mutate(type = "pDNA", rank = dense_rank(desc(RPM)), cum_RPM = cumsum(RPM) / sum(RPM), RPM_100 = RPM / 1000)

df.gDNA <- df %>%
    select(bc, gRPM) %>%
    rename(RPM = gRPM) %>%
    arrange(-RPM) %>%
    mutate(type = "gDNA", rank = dense_rank(desc(RPM)), cum_RPM = cumsum(RPM) / sum(RPM), RPM_100 = RPM / 1000)

df2 <- rbind(df.pDNA, df.gDNA)

# use barodes identified in both pDNA and gDNA
plt.ecdf <- ggplot(df2, aes(x = RPM_100)) +
    stat_ecdf(aes(colour = type), size = 1, geom = "line") +
    labs(x = "Barcode count (x1,000)", y = "Cumulative barcode frequency") +
    scale_colour_manual(values = c("pDNA" = "#00A4CCFF", "gDNA" = "#F95700FF")) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0.01)) +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE),
    ) +
    guides(colour = "none")
ggsave(plot = plt.ecdf, filename = "plots/barcode_ecdf_01082023.pdf", w = 1.75, h = 1.65)