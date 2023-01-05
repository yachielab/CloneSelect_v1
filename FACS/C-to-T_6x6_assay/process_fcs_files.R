library(readxl)
library(dplyr)
library(CytoExploreR)
library(flowWorkspace)

# Load metadata from excel file
sample_sheet <- read_excel(path = "04112022_CA062_SampleSheet.xlsx")
sample_sheet <- sample_sheet %>%
    mutate(is_match = ifelse(Cell_line == Query_gRNA, "OT", "NT"))

# Load flow samples
cs <- load_cytoset_from_fcs(path = "./fcs")
gs <- GatingSet(cs)

trans.log <- cyto_transformer_log(gs, channels = c("FL1-A", "FL5-A"), plot = FALSE)
gs.trans <- cyto_transform(gs, trans = trans.log, plot = FALSE)

sample.df <- pData(gs.trans)
merged.df <- merge(sample_sheet, sample.df, by.x = "Filename", by.y = "name")

pData(gs.trans)$System <- merged.df$System
pData(gs.trans)$Rep <- merged.df$Rep
pData(gs.trans)$Cell_line <- merged.df$Cell_line
pData(gs.trans)$Query_gRNA <- merged.df$Query_gRNA
pData(gs.trans)$Cell_line2 <- merged.df$Cell_line2
pData(gs.trans)$Query_gRNA2 <- merged.df$Query_gRNA2
pData(gs.trans)$Assay_group <- merged.df$Assay_group
pData(gs.trans)$is_match <- merged.df$is_match

# Manual gating
is_interactive_gate <- FALSE
if (is_interactive_gate) {
    # Live cell gate
    cyto_gate_draw(
        gs.trans,
        parent = "root",
        alias = "P1",
        channels = c("FSC-A", "SSC-A"),
        gatingTemplate = "CA062_GatingTemplate_04112022.csv",
        type = "polygon",
        xlim = c(-10, 1 * 10^6),
        ylim = c(-10, 1 * 10^6),
    )

    # Single-cell gate
    cyto_gate_draw(
        gs.trans,
        parent = "P1",
        alias = "P2",
        channels = c("FSC-A", "FSC-Width"),
        gatingTemplate = "CA062_GatingTemplate_04112022.csv",
        xlim = c(-10, 1 * 10^6),
        ylim = c(-10, 6 * 10^3),
        type = "polygon"
    )

    # GFP gate
    # cyto_gate_draw(
    cyto_gate_edit(
        gs.trans["06-CtoT-A4.fcs"],
        parent = "P2",
        alias = "GFP",
        channels = c("FL1-A", "FSC-A"),
        gatingTemplate = "CA062_GatingTemplate_04112022.csv",
        xlim = c(-1e2, 3 * 10^6),
        ylim = c(-1e2, 1 * 10^6),
        type = "interval"
    )
}

# Apply gate for all samples
gt <- gatingTemplate("CA062_GatingTemplate_04112022.csv")
if ("/P1/P2/GFP" %in% gs_get_pop_paths(gs.trans)) {
    gs_pop_remove(gs.trans, "P1") # remove existing gates
    gt_gating(gt, gs.trans)
} else {
    gt_gating(gt, gs.trans)
}

# percent EGFP+
GFP.stat <- gs_pop_get_stats(gs.trans, node = "GFP", type = "percent") %>% mutate(percent = percent * 100)

# mediam EGFP expr
GFP.MFI <- gs_pop_get_stats(gs.trans, nodes = c("GFP"), type = pop.MFI, inverse.transform = TRUE)

df.stat <- merge(pData(gs.trans), c(GFP.stat, GFP.MFI), by.x = "name", by.y = "sample")

# export data
write.csv(x = df.stat, file = "GFP_stat_CA062_04112022.csv", quote = FALSE, row.names = FALSE)

# save_gs(gs.trans, path = "CS062_gs", verbose = TRUE)