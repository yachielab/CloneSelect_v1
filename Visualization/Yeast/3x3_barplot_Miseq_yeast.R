library(ggplot2)
library(dplyr)

df <- read.table("data/20190516_Target-AID_BY4741_3x3_Miseq.tsv", header = TRUE, sep = "\t")
df <- df %>% mutate(is_match = ifelse(gRNA == Cell, "OT", "NT"))

FONT.SIZE <- 7
FONT.SIZE.LABEL <- 8
LINE.W <- 0.232

g <- ggplot(df, aes(x = gRNA, y = ATG.conversion * 100)) +
    facet_wrap(~Cell) +
    geom_bar(aes(fill = is_match), stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("OT" = "darkblue", "NT" = "gray")) +
    theme_classic() +
    scale_x_discrete(expand = c(0, 0), labels = c("BC2" = "yBC-C1", "BC4" = "yBC-C2", "BC8" = "yBC-C3")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +
    theme(
        panel.spacing.x = unit(0.8, "mm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE.LABEL),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.text.x = element_text(colour = "black", size = FONT.SIZE, angle = 90, vjust = 0.5, hjust = 1),
    ) +
    labs(x = "", y = "") +
    guides(fill = "none")
ggsave(plot = g, filename = "./plots/3x3_RFP_miseq_08172022.pdf", w = 1.5, h = 1.9, device = cairo_pdf)