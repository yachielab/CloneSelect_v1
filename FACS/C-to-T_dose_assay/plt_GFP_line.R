library(RColorBrewer)
library(ggplot2)
library(ggridges)
library(ggforce)
library(dplyr)
library(scales)
library(data.table)
library(rmarkdown)
library(openCyto)
library(flowWorkspace)
library(CytoExploreR)

# loading FCSs
cytoset <- load_cytoset_from_fcs(path = "./fcs")
gs <- GatingSet(cytoset)

# Transformation
trans.log <- cyto_transformer_log(gs, channels = c("FITC-A"), plot = FALSE)
gs.trans <- cyto_transform(gs, trans = trans.log, plot = FALSE)

cyto_details(gs.trans)$exp_id <- "CA044"
cyto_details(gs.trans)$dose_str <- c(
    rep("0ng", 6), rep("100ng", 6), rep("200ng", 6),
    rep("400ng", 6), rep("50ng", 6), rep("600ng", 6),
    rep("800ng", 6)
)
cyto_details(gs.trans)$dose_ng <- c(
    rep(0, 6), rep(100, 6), rep(200, 6),
    rep(400, 6), rep(50, 6), rep(600, 6),
    rep(800, 6)
)
cyto_details(gs.trans)$OT <- rep(c("NC", "PC"), 21)
cyto_details(gs.trans)$barcode <- rep(c("BC2", "BC2", "BC4", "BC4", "BC6", "BC6"), 7)

gt <- gatingTemplate("CA044_manual_gates.csv")
gt_gating(gt, gs.trans)

stat <- flowWorkspace::gs_pop_get_stats(gs.trans, "GFP", type = "percent")
sample.meta <- cyto_details(gs.trans)
stat.merged.df <- merge(stat, sample.meta, by.x = "sample", by.y = "name")

FONT.SIZE <- 7
LABEL.FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt
DATE <- "04192022"

g.line.plt <- ggplot(
    stat.merged.df %>% dplyr::filter(OT == "PC"),
    aes(
        x = as.numeric(dose_ng),
        y = percent, colour = barcode, group = barcode
    )
) +
    geom_line(alpha = 0.8) +
    geom_point(size = 1.8) +
    scale_x_continuous(expand = c(0.05, 0.01), limits = c(0, 800)) +
    scale_y_continuous(
        limits = c(0, 0.3), expand = c(0.01, 0),
        labels = scales::percent_format(suffix = "")
    ) +
    labs(x = "Input DNA (ng)", y = "%EGFP positive cells") +
    scale_colour_manual(values = c("BC2" = "#D53E4F", "BC4" = "#FA8E59", "BC6" = "#fee08b")) +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(colour = "black", size = FONT.SIZE),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(colour = "black", size = LINE.W),
        axis.text = element_text(colour = "black", size = LABEL.FONT.SIZE)
    )

ggsave(plot = g.line.plt + guides(colour = "none"), filename = "plots/GFP_line_09242022.pdf", w = 2.4, h = 2.0, device = cairo_pdf)