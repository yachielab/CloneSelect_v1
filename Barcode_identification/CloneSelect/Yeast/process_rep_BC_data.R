library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)

normalize <- function(df, outname) {
    df <- df %>% filter(rep1 != "NA", rep2 != "NA", rep3 != "NA") # Use barcodes found in all replicates

    # Total count
    rep1_cnt <- sum(df$rep1)
    rep2_cnt <- sum(df$rep2)
    rep3_cnt <- sum(df$rep3)

    df <- df %>%
        group_by(barcode) %>%
        mutate(
            rep1_freq = rep1 / rep1_cnt, rep2_freq = rep2 / rep2_cnt, rep3_freq = rep3 / rep3_cnt,
            mean_freq = mean(c(rep1_freq, rep2_freq, rep3_freq)),
            rep1_RPM = rep1 * 1e6 / rep1_cnt, rep2_RPM = rep2 * 1e6 / rep2_cnt, rep3_RPM = rep3 * 1e6 / rep3_cnt,
            mean_RPM = mean(c(rep1_RPM, rep2_RPM, rep3_RPM))
        ) %>%
        arrange(-mean_RPM)

    write.csv(x = df, file = outname, quote = FALSE, row.names = FALSE)
}

df100 <- read.csv("data/barcode_count_merged_analysis_BC_pool_100_pKI086_mix.tsv", sep = "\t", header = TRUE)
df1580 <- read.csv("data/barcode_count_merged_analysis_BC_pool_1580_pKI086_mix.tsv", sep = "\t", header = TRUE)

normalize(df1580, "data/BC1580_normalized_count_06012022.csv")
normalize(df100, "data/BC100_normalized_count_06012022.csv")

