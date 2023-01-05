library(ggplot2)
library(dplyr)

# define const for visualization
FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt

plt.bar <- function(df, ...) {
    arglist <- list(...) # additional ggplot2 object(s)
    g <- ggplot(df) +
        facet_wrap(~Cell_line) +
        geom_bar(aes(group = gRNA, fill = col, x = gRNA, y = percent), position = "dodge", stat = "summary", fun = "mean") +
        geom_point(aes(y = percent, x = gRNA), colour = "black", size = 2.0, shape = 1, stat = "identity") +
        # scale_fill_manual(values = c("NT_gRNA" = "gray", "BC2_gRNA" = "green3")) +
        scale_fill_identity() +
        scale_x_discrete(expand = c(0, 0), labels = c("", "")) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
        theme_classic() +
        theme(
            panel.background = element_blank(),
            plot.background = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5, size = FONT.SIZE),
            axis.text.y = element_text(colour = "black", size = FONT.SIZE),
            axis.title = element_text(size = FONT.SIZE),
            axis.line = element_line(size = LINE.W, colour = "black"),
            axis.ticks = element_line(size = LINE.W, colour = "black"),
        ) +
        labs(y = "EGFP positive cells (%)", x = "") +
        guides(fill = "none")

    # overwrite existing layers
    if (length(arglist) > 0) {
        return(g + arglist)
    } else {
        return(g)
    }
}

df.stat <- read.csv("pop_stat_03152022.csv", header = TRUE)

# reporter cells
df.GTG <- df.stat %>% filter(Sample_type == "Treatment")
g.bar.GTG <- plt.bar(df.GTG)

# for ATG ctrl cells
df.ATG <- df.stat %>% filter(Sample_type == "Control")
g.bar.ATG <- plt.bar(df.ATG) + scale_y_continuous(expand = c(0, 0), limits = c(0, 100))

w <- 1.1
h <- 1.7
ggsave(plot = g.bar.GTG, filename = "./plots/GFP_barplt_GTG_cells_05022022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = g.bar.ATG, filename = "./plots/GFP_barplt_ATG_cells_05022022.pdf", w = w, h = h, device = cairo_pdf)