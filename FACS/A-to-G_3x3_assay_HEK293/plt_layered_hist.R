library(ggplot2)
library(data.table)
library(dplyr)
library(ggridges)
library(scales)

FONT.SIZE <- 7
LABEL.FONT.SIZE <- 5
LINE.W <- 0.232 # 0.5pt

# Load files
exported_file <- "./processed_csv/single_cell_events_ABE_CA060_02172022.csv"
df <- fread(exported_file, header = TRUE, sep = ",")
df.abe <- df %>% filter(system == "A-to-G_BC")
df.gRNA <- df %>% filter(system == "gRNA_BC")
df.CRISPRa <- df %>% filter(system == "CRISPRa_BC")

inv.log.trans <- function(x) {
    ts <- flowWorkspace::flowjo_log_trans()
    inv <- ts[["inverse"]](x)
    return(inv)
}

# Access GFP gate object
gt <- openCyto::gatingTemplate("./CA059_02062022_Gate.csv")
GFP.gate <- gt@edgeData@data[["/P1|/P1/GFP"]]$gtMethod@args$gate@.Data[["all"]]@.Data[[1]]
GFP.cutoff <- GFP.gate@min[["FL1-A"]] %>% inv.log.trans()

plt.hist <- function(data) {
    g <- ggplot(data, aes(y = rep, x = `B525-FITC-A`, fill = stat(x))) +
        # facet_wrap(cell_line ~ tx_reagent2) + # for debug
        facet_wrap(cell_line ~ tx_reagent) + # x=cell_line, y=tx_reagent (gRNA or reporter)
        geom_density_ridges_gradient(
            colour = "white", scale = 2.5, quantile_lines = TRUE,
            quantiles = 2, vline_size = 0.5, vline_color = "gray", size = 0.2
        ) +
        scale_fill_gradientn(limits = c(1, 7), colours = c("black", "gray20", "green2", "greenyellow")) +
        geom_vline(xintercept = GFP.cutoff, linetype = "dashed", size = LINE.W, colour = "gray20") +
        scale_x_log10(
            breaks = c(10^1, 10^3, 10^5, 10^7),
            labels = scales::trans_format("log10", math_format(10^.x)),
            expand = c(0, 0),
            limits = c(10, 10000000)
        ) +
        scale_y_discrete(expand = c(0, 0)) +
        labs(x = "EGFP intensity (a.u.)", y = "") +
        theme_classic() +
        theme(
            plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
            panel.spacing.y = unit(0.2, "lines"), # margin between vertical direction
            panel.spacing.x = unit(0.9, "lines"), # margen between horizontal direction
            panel.background = element_rect(fill = "transparent", colour = "transparent"),
            plot.background = element_rect(fill = "transparent", colour = "transparent"),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(colour = "black", size = LINE.W),
            axis.title = element_text(size = FONT.SIZE),
            axis.text.x = element_text(colour = "black", size = FONT.SIZE),
            axis.text.y = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            strip.background = element_blank()
        ) +
        guides(fill = "none")
    return(g)
}

# debug
# df.abe <- df.abe %>% mutate(tx_reagent2 = paste0(tx_reagent, "_gRNA"))

# p.hist.abe <- plt.hist(df.abe)
# ggsave(plot = p.hist.abe, filename = "./plots/ABE_GFP_histograms_07182022.pdf", w = fig.w, h = fig.h)
# stop("ok")

p.hist.gRNA <- plt.hist(df.gRNA)
p.hist.CRISPRa <- plt.hist(df.CRISPRa)

fig.w <- 2.7
fig.h <- 2.4

ggsave(plot = p.hist.abe, filename = "./plots/ABE_GFP_histograms_07182022.pdf", w = fig.w, h = fig.h)
ggsave(plot = p.hist.gRNA, filename = "./plots/gRNA_GFP_histograms_07182022.pdf", w = fig.w, h = fig.h)
ggsave(plot = p.hist.CRISPRa, filename = "./plots/CRISPRa_GFP_histograms_07182022.pdf", w = fig.w, h = fig.h)