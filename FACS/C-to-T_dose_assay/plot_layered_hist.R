library(data.table)
library(ggplot2)
library(ggridges)
library(openCyto)
library(dplyr)

FONT.SIZE <- 8
LABEL.FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt
DATE <- "09242022"

exported_file <- "./processed_csv/single_cell_events_CA044.csv"
df <- fread(exported_file, header = TRUE, sep = ",")

inv.log.trans <- function(x) {
    ts <- flowWorkspace::flowjo_log_trans()
    inv <- ts[["inverse"]](x)
    return(inv)
}

gt <- gatingTemplate("CA044_manual_gates.csv")
GFP.gate <- gt@edgeData@data[["/P1|/P1/GFP"]]$gtMethod@args$gate@.Data$`Combined Events`@.Data[[1]]
GFP.cutoff <- GFP.gate@min[["FITC-A"]] %>% inv.log.trans()

# rename barcode name
supp.labs <- c("BC-C1", "BC-C2", "BC-C3")
names(supp.labs) <- c("BC2", "BC4", "BC6")

plt.hist <- function(data) {
    g <- ggplot(data, aes(y = dose_str, x = `FITC-A`, fill = stat(x))) +
        facet_wrap(~barcode, labeller = as_labeller(supp.labs)) +
        geom_density_ridges_gradient(
            na.rm = TRUE,
            colour = "white", scale = 2.5, quantile_lines = TRUE,
            quantiles = 2, vline_size = 0.5, vline_color = "gray", size = 0.2
        ) +
        scale_fill_gradientn(limits = c(0, 5), colours = c("black", "gray10", "green2", "green2", "greenyellow")) +
        geom_vline(xintercept = GFP.cutoff, linetype = "dashed", size = LINE.W, colour = "gray20") +
        scale_x_log10(
            breaks = c(10^0, 10^2, 10^4),
            labels = scales::trans_format("log10", scales::math_format(10^.x)),
            expand = c(0, 0),
            limits = c(0.5, 10^5)
        ) +
        scale_y_discrete(expand = c(0, 0)) +
        labs(x = "EGFP intensity (a.u.)", y = "") +
        theme_classic() +
        theme(
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.spacing.y = unit(0.5, "lines"),
            panel.spacing.x = unit(0.8, "lines"),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(size = LINE.W),
            axis.title = element_text(size = FONT.SIZE),
            axis.text.x = element_text(colour = "black", size = LABEL.FONT.SIZE),
            axis.text.y = element_text(colour = "black", size = LABEL.FONT.SIZE),
            strip.text.x = element_text(colour = "black", size = LABEL.FONT.SIZE),
            strip.background = element_blank()
        )
    return(g)
}

# Reorder with proper order
df$dose_str <- factor(df$dose_str, levels = c("0ng", "50ng", "100ng", "200ng", "400ng", "600ng", "800ng"))
g.hist <- plt.hist(df %>% dplyr::filter(.data[["FITC-A"]] > 1))
ggsave(plot = g.hist + guides(fill = "none"), filename = sprintf("GFP_hist_%s.pdf", DATE), w = 2.9, h = 2.45, device = cairo_pdf)