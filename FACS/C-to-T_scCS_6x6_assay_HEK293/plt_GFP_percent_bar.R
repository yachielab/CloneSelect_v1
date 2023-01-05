library(ggplot2)
library(scales)
library(dplyr)
library(lemon)


FONT.SIZE <- 8
LABEL.FONT.SIZE <- 7
LINE.W <- 0.232 # equivalent to 0.5pt in Keynote

df <- read.csv("pop_stat_03152022.csv", header = TRUE)

# barcode order
df$Cell_line <- factor(df$Cell_line, c("BC2", "BC4", "BC6", "BC8", "BC10", "BC12"))
df$gRNA <- factor(df$gRNA, c("BC2gRNA", "BC4gRNA", "BC6gRNA", "BC8gRNA", "BC10gRNA", "BC12gRNA"))

df %>% mutate(
    Cell_line2 =
        case_when(
            Cell_line == "BC2" ~ "BC-C1",
            Cell_line == "BC4" ~ "BC-C2",
            Cell_line == "BC6" ~ "BC-C3",
            Cell_line == "BC8" ~ "BC-C4",
            Cell_line == "BC10" ~ "BC-C5",
            Cell_line == "BC12" ~ "BC-C6",
        )
) -> df

df <- df %>% arrange(gRNA)

g1 <- ggplot(df %>% dplyr::filter(Target_type != "Control"), aes(y = percent, x = Cell_line2)) +
    facet_rep_wrap(~gRNA, ncol = 6, repeat.tick.labels = "left") +
    # facet_grid(~gRNA) +
    geom_bar(aes(fill = Target_type), stat = "identity", width = 0.9) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = c("Match" = "green3", "Mismatch" = "gray")) +
    labs(x = "", y = "EGFP positive cells (%)") +
    theme(
        panel.background = element_rect(fill = "transparent", colour = "transparent"),
        plot.background = element_rect(fill = "transparent", , colour = "transparent"),
        panel.spacing = unit(-0.1, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = LABEL.FONT.SIZE, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(colour = "black", size = LABEL.FONT.SIZE),
        axis.title = element_text(size = FONT.SIZE),
        axis.line = element_line(size = LINE.W, colour = "black"),
        axis.ticks = element_line(size = LINE.W, colour = "black"),
        # axis.ticks.x = element_blank()
    ) +
    guides(fill = "none")

ggsave(plot = g1, filename = "GFP_barplot_09242022.pdf", h = 2.5, w = 4.5, device = cairo_pdf)