library(tidyr)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

df.rep1 <- read.csv("data/GFP_isolation/rep1_result/barcode_counts_rep1.tsv", header = TRUE, sep = "\t")
df.rep2 <- read.csv("data/GFP_isolation/rep2_result/barcode_counts_rep2.tsv", header = TRUE, sep = "\t")

whitelist <- read_excel("merged_96lib_pDNA_gDNA_01082023.xlsx")
whitelist <- whitelist %>% mutate(clone_name2 = str_replace(clone_name, " ", "."))

# rep1
# maxtrix to long format and normalize
df2.rep1 <- df.rep1 %>%
    pivot_longer(cols = !clone_name, values_to = "raw_count") %>%
    mutate(clone_name = str_replace(clone_name, " ", ".")) %>%
    mutate(is_match = ifelse(clone_name == name, "OT", "NT"))

df.total.rep1 <- df2.rep1 %>%
    group_by(clone_name) %>%
    summarise(total = sum(raw_count))

df.heatmap.rep1 <- merge(df2.rep1, df.total.rep1, by = "clone_name") %>%
    mutate(norm_bc_freq = raw_count / total) %>%
    mutate(RPM = raw_count * 1e6 / total)

# rep2
# maxtrix to long format and normalize
df2.rep2 <- df.rep2 %>%
    pivot_longer(cols = !clone_name, values_to = "raw_count") %>%
    mutate(clone_name = str_replace(clone_name, " ", ".")) %>%
    mutate(is_match = ifelse(clone_name == name, "OT", "NT"))

df.total.rep2 <- df2.rep2 %>%
    group_by(clone_name) %>%
    summarise(total = sum(raw_count))

df.heatmap.rep2 <- merge(df2.rep2, df.total.rep2, by = "clone_name") %>%
    mutate(norm_bc_freq = raw_count / total) %>%
    mutate(RPM = raw_count * 1e6 / total)

df.merged <- merge(df.heatmap.rep1, df.heatmap.rep2, by = c("clone_name", "is_match"))



g.scatter <- ggplot(df.merged, aes(x = RPM.x, y = RPM.y)) +
    geom_point(aes(colour = is_match)) +
    theme_classic() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.text.x = element_text(colour = "black", size = FONT.SIZE),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.line = element_line(colour = "black", size = LINE.W),
        axis.ticks = element_line(colour = "black", size = LINE.W)
    ) +
    scale_color_manual(values = c("OT" = "deeppink2", "NT" = "gray")) +
    labs(x = "", y = "") +
    scale_x_continuous(limits = c(0, 1e6)) +
    scale_y_continuous(limits = c(0, 1e6)) +
    guides(colour = "none")
ggsave(plot = g.scatter, filename = "scatter_rep1_rep2_01112023.png", w = 3, h = 2.9)