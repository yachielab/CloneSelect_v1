library(ggplot2)
library(dplyr)

FONT.SIZE <- 7
LINE.W <- 0.232

# Barcode color map
# color_map <- c("BC2" = "#2596be", "BC4" = "#be2596", "BC6" = "#96be25")

plt.bar <- function(df, ...) {
    arglist <- list(...) # additional ggplot2 object(s)
    g <- ggplot(df, aes(x = Cell_line, y = percent)) +
        facet_wrap(~gRNA) +
        geom_bar(aes(group = Cell_line, fill = Target_type), position = "dodge", stat = "summary", fun = "mean") +
        geom_point(colour = "black", size = 2.0, shape = 1, stat = "identity") +
        scale_fill_manual(values = c("OT" = "green3", "NT" = "gray")) +
        scale_x_discrete(expand = c(0, 0), labels = c("BC2" = "BC-C1", "BC4" = "BC-C2", "BC6" = "BC-C3")) +
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

df.stat <- read.csv("GFP_activation_03172022.csv", header = TRUE)

# subset data.frame
df.lenti <- df.stat %>% dplyr::filter(Delivery_method == "Lenti", Sample_type == "Treatment")
df.tx <- df.stat %>% dplyr::filter(Delivery_method == "Transfection", Sample_type == "Treatment")

w <- 1.9
h <- 2.1
g.bar.lenti <- plt.bar(df.lenti, ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 50)))
g.bar.tx <- plt.bar(df.tx, ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 50)))
ggsave(plot = g.bar.lenti, filename = "./plots/lenti_GFP_mESC_05022022.pdf", w = w, h = h, device = cairo_pdf)
ggsave(plot = g.bar.tx, filename = "./plots/tx_GFP_mESC_05022022.pdf", w = w, h = h, device = cairo_pdf)