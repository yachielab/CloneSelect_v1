library(ggplot2)
library(scales)
library(dplyr)

tsvFiles <- list.files(pattern = ".fcs", path = "./event_data/", full.names = TRUE)

datalist <- list()
for (tsv.file in tsvFiles) {
    sample <- unlist(strsplit(basename(tsv.file), "_"))[1]
    df <- read.csv(tsv.file, header = TRUE, skip = 1, sep = "\t")
    df$sample <- sample
    datalist[[sample]] <- df
}
agg.data <- do.call(rbind, datalist)

## each data
df.192 <- datalist[1]
df.193 <- datalist[2]
df.194 <- datalist[3]
df.195 <- datalist[4]
df.196 <- datalist[5]
df.197 <- datalist[6]

## GTG with ATG control
df1 <- rbind(df.195[[1]], df.192[[1]])
df2 <- rbind(df.196[[1]], df.193[[1]])
df3 <- rbind(df.197[[1]], df.194[[1]])

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

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

plt.hist <- function(x) {
    g <- ggplot(x, aes(x = Y610.mCHERRY.A_.FL5.A.)) +
        geom_histogram(aes(fill = sample, alpha = sample), colour = NA, bins = 200, position = "identity") +
        scale_x_log10(
            breaks = c(1e0, 1e2, 1e4, 1e6),
            expand = c(0, 0),
            limits = c(1, 1000000),
            labels = scales::trans_format("log10", math_format(10^.x))
        ) +
        theme_classic() +
        theme(
            plot.background = element_rect(fill = "transparent", colour = "transparent"),
            panel.background = element_rect(fill = "transparent", colour = "transparent"),
            plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
            axis.ticks = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(size = LINE.W),
            axis.title = element_text(size = FONT.SIZE),
            axis.text = element_text(colour = "black", size = FONT.SIZE)
        ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 600)) +
        scale_fill_manual(values = c("gray70", "firebrick2")) +
        scale_alpha_manual(values = c(1.0, 0.6), guide = FALSE) +
        labs(x = "", y = "") +
        guides(fill = "none", colour = "none")
    return(g)
}

w <- 1.4
h <- 1.8

ggsave(plot = plt.hist(df1), filename = "./plots/192_195_hist_09162021.pdf", w = w, h = h)
ggsave(plot = plt.hist(df2), filename = "./plots/193_196_hist_09162021.pdf", w = w, h = h)
ggsave(plot = plt.hist(df3), filename = "./plots/194_197_hist_09162021.pdf", w = w, h = h)
