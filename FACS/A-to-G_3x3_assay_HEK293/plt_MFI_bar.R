library(ggplot2)
library(data.table)
library(dplyr)

FONT.SIZE <- 7
LABEL.FONT.SIZE <- 5
LINE.W <- 0.232 # 0.5pt

# Load files
exported_file <- "./processed_csv/single_cell_events_ABE_CA060_02172022.csv"
df <- fread(exported_file, header = TRUE, sep = ",")
df <- data.frame(df)
df %>%
    filter(is_matched == "OT", B525.FITC.A > GFP.cutoff) %>%
    group_by(system, rep, cell_line) %>%
    summarise(MFI = median(B525.FITC.A)) %>%
    mutate(type = "OT") -> df.OT

g.MFI <- ggplot(df.OT, aes(y = MFI, x = cell_line)) +
    facet_wrap(~system) +
    # geom_bar(aes(group = type, fill = system), position = "dodge", stat = "summary", fun = "median") +
    geom_point(
        aes(group = type, fill = cell_line),
        size = 2, shape = 21,
        position = position_dodge(width = 0.9)
    ) +
    scale_y_log10(
        breaks = c(10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        expand = c(0, 0),
        limits = c(1, 10^6)
    ) +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, size = FONT.SIZE),
        axis.text.y = element_text(colour = "black", size = FONT.SIZE),
        axis.title = element_text(size = FONT.SIZE),
        axis.line = element_line(size = LINE.W, colour = "black"),
        axis.ticks.y = element_line(size = LINE.W, colour = "black"),
        axis.ticks.x = element_blank()
    ) +
    labs(x = "", y = "Median fluorescence intensity (a.u.)") +
    guides(fill = "none")
ggsave(plot = g.MFI, filename = sprintf("plots/GFP_MFI_%s.pdf", "09242022"), w = 1.5, h = 2.0, device = cairo_pdf)