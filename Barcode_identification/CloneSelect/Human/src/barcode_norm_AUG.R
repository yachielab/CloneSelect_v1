library(tidyr)
library(dplyr)
library(readxl)

FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

df <- read.csv("data/GFP_isolation/rep1_result/barcode_counts_rep1_only_atg.tsv", header = TRUE, sep = "\t")

whitelist <- read_excel("merged_96lib_pDNA_gDNA_01082023.xlsx")
whitelist <- whitelist %>% mutate(clone_name2 = str_replace(clone_name, " ", "."))

# maxtrix to long format
df2 <- df %>%
    pivot_longer(cols = !clone_name, values_to = "raw_count") %>%
    mutate(clone_name = str_replace(clone_name, " ", ".")) %>%
    mutate(is_match = ifelse(clone_name == name, "OT", "NT"))

df.total <- df2 %>%
    group_by(clone_name) %>%
    summarise(total = sum(raw_count))

df.heatmap <- merge(df2, df.total, by = "clone_name") %>%
    mutate(norm_bc_freq = raw_count / total) %>%
    # mutate(rank = as.numeric(str_replace(clone_name, "Clone.", "")))
    arrange(name, is_match) %>%
    mutate(new_order = as.numeric(str_replace(clone_name, "Clone.", "")))
df.merged <- merge(x = df.heatmap, y = whitelist %>% select(clone_name2, order, gDNA_freq, Is_used_for_isolation), by.x = "name", by.y = "clone_name2")


g.heatmap.all <- ggplot(df.merged, aes(y = reorder(clone_name, -new_order), x = reorder(name, order))) +
    geom_tile(aes(fill = norm_bc_freq), colour = NA) +
    theme_classic() +
    theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, size = FONT.SIZE), # element_blank(),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.line = element_blank(),
        panel.border = element_blank()
    ) +
    # scale_fill_gradientn(
    #    colours = c("gray98", "dodgerblue", "blue", "blue2", "darkblue"),
    #    limits = c(0, 1), na.value = "white"
    # ) +
    scale_fill_gradientn(
        colours = c("gray97", "#c4f4b9", "green2"),
        limits = c(0, 1), na.value = "white"
    ) +
    labs(x = "", y = "") +
    guides(fill = "none")
ggsave(plot = g.heatmap.all, filename = "BC_heatmap_AUG_01112023.pdf", device = cairo_pdf, w = 5.2, h = 2.6)

g.jitter <- ggplot(df.merged, aes(x = reorder(clone_name, -new_order), y = norm_bc_freq)) +
    geom_jitter(data = df.merged %>% filter(is_match == "NT"), aes(colour = is_match), width = 0.3, size = 0.7, alpha = 0.5, shape = 16) +
    geom_jitter(data = df.merged %>% filter(is_match == "OT"), aes(colour = is_match), width = 0.1, size = 1.5, shape = 16) +
    scale_y_continuous(position = "right") +
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
ggsave(plot = g.jitter, filename = "BC_jitter_AUG_01112022_.pdf", device = cairo_pdf, w = 2.21, h = 1.28)