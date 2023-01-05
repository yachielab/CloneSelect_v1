library(CytoExploreR)
library(flowWorkspace)
library(dplyr)
library(data.table)

# Loading dataset
cs <- load_cytoset_from_fcs(path = "./fcs")
gs <- GatingSet(cs)
sampleNames(gs)
sampleNames(gs) %>% length() # sample numbers

# Transformation
trans.log <- cyto_transformer_log(gs, channels = c("FL1-A", "FL5-A"), plot = FALSE)
gs.trans <- cyto_transform(gs, trans = trans.log, plot = FALSE)

# Manual gate
gt <- gatingTemplate("./CA059_02062022_Gate.csv")
gt_gating(gt, gs.trans)

# Add experimenal metadata
cyto_details(gs.trans)$rep <- rep(c("Rep1", "Rep2", "Rep3"), 27)
cyto_details(gs.trans)$exp_id <- c(rep("CA060", 27), rep("CA059", 54))
cyto_details(gs.trans)$system <- c(rep("gRNA_BC", 27), rep("A-to-G_BC", 27), rep("CRISPRa_BC", 27))

lenti_plasmid <- c(
    rep("CS-254", 9),
    rep("CS-255", 9),
    rep("CS-256", 9),
    rep("CS-121", 9),
    rep("CS-246", 9),
    rep("CS-247", 9),
    rep("CS-251", 9),
    rep("CS-252", 9),
    rep("CS-253", 9)
)
cyto_details(gs.trans)$lenti_plasmid <- lenti_plasmid

cell_line <- c(
    rep("T3", 9),
    rep("T4", 9),
    rep("T5", 9),
    rep("T3", 9),
    rep("T4", 9),
    rep("T5", 9),
    rep("T3", 9),
    rep("T4", 9),
    rep("T5", 9)
)
cyto_details(gs.trans)$cell_line <- cell_line

# Tx plasmid reagent of either gRNA or CRISPRa plasmid (query plasmid seq)
cyto_details(gs.trans)$tx_reagent <- rep(rep(c(rep("T3", 3), rep("T4", 3), rep("T5", 3)), 3), 3)

is_matched <- c(
    rep("OT", 3), rep("NT", 3), rep("NT", 3),
    rep("NT", 3), rep("OT", 3), rep("NT", 3),
    rep("NT", 3), rep("NT", 3), rep("OT", 3)
)
cyto_details(gs.trans)$is_matched <- rep(is_matched, 3)

# merge stats and metadata
sample.names <- sampleNames(gs.trans)
sample.stats <- gs_pop_get_stats(gs.trans, type = "percent", nodes = c("GFP")) %>% mutate(percent = percent * 100)
sample.meta <- cyto_details(gs.trans)
stat.df <- merge(sample.stats, sample.meta, by.x = "sample", by.y = "name")

# export single-cell event data to csv
gs.exp.df <- data.table()
for (sample in sample.names) {
    message(sprintf("Processing: %s", sample))
    df <- gs_get_singlecell_expression(gs.trans, c("P1"),
        threshold = FALSE, marginal = FALSE, inverse.transform = TRUE,
        other.markers = c("SSC-A", "FL1-A", "FL5-A")
    )[[sample]] %>%
        data.table() %>%
        mutate(sample = sample)
    gs.exp.df <- rbindlist(list(gs.exp.df, df))
}
exp.df <- merge(gs.exp.df, sample.meta, by.x = "sample", by.y = "name")
write.csv(sample.meta, quote = FALSE, row.names = FALSE, file = "./processed_csv/sample_meta_02172022.csv")
write.csv(exp.df, quote = FALSE, row.names = FALSE, file = "./processed_csv/single_cell_events_ABE_CA060_02172022.csv")

# export GatingSet object to the disk
gs_save_dir <- "./gs_achieve/gs_trans_CA060"
if (file.exists(gs_save_dir)) {
    message(sprintf("Removed old objects %s", gs_save_dir))
    unlink(gs_save_dir, recursive = TRUE)
} else {
    message(sprintf("Created %s", gs_save_dir))
    dir.create(gs_save_dir)
}
save_gs(gs.trans, path = file.path(gs_save_dir), verbose = TRUE)

