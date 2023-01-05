library(flowWorkspace)
library(CytoExploreR, include.only = "cyto_details")
library(tidyr)
library(dplyr)

inv.log.trans <- function(x) {
    ts <- flowWorkspace::flowjo_log_trans()
    inv <- ts[["inverse"]](x)
    return(inv)
}
log.trans <- function(x) {
    ts <- flowWorkspace::flowjo_log_trans()
    trs <- ts[["transform"]](x)
    return(trs)
}

compute_dynamic_gate <- function(gs, inv.bins, parent_gate = "P1") {
    # We assume dynamic gate will be created for P1 population.
    i <- 1

    df.stat <- data.frame()
    for (min.cutoff in inv.bins) {
        message(sprintf("Processing of %f, %d of %d...", min.cutoff, i, length(inv.bins)))

        # Create dynamic FL1-A gate
        rg_gfp <- flowCore::rectangleGate(
            "FL1-A" = c(min.cutoff, Inf), "FSC-A" = c(-Inf, Inf),
            filterId = "GFP_dynamic"
        )
        gates <- flowWorkspace::gs_get_pop_paths(gs)
        if ("/P1/GFP_dynamic" %in% gates) {
            flowWorkspace::gs_pop_remove("GFP_dynamic", gs = gs)
            flowWorkspace::gs_pop_add(gs, rg_gfp, parent = parent_gate)
            flowWorkspace::recompute(gs)
        } else {
            flowWorkspace::gs_pop_add(gs, rg_gfp, parent = parent_gate)
            flowWorkspace::recompute(gs)
        }
        stat <- flowWorkspace::gs_pop_get_stats(gs, "GFP_dynamic", type = "percent") %>%
            dplyr::mutate(percent, cutoff = inv.log.trans(min.cutoff))
        df.stat <- rbind(df.stat, stat)
        i <- i + 1
    }
    return(df.stat)
}

# load data
gs.trans <- load_gs("./gs_achieve/gs_trans_CA060", verbose = TRUE)

gs.positive <- subset(gs.trans, is_matched == "OT")
gs.negative <- subset(gs.trans, is_matched == "NT")

bins <- seq(0, 1e5, by = 120)
bins.trans <- log.trans(bins)

# apply dynamic gates
df.positive <- compute_dynamic_gate(gs.positive, bins.trans[1:800])
df.negative <- compute_dynamic_gate(gs.negative, bins.trans[1:800])

# merge data.frame
sample.meta.pos <- cyto_details(gs.positive)
sample.meta.neg <- cyto_details(gs.negative)
stat.df.pos <- merge(df.positive, sample.meta.pos, by.x = "sample", by.y = "name")

# Group two expected negative cell lines and calculate mean GFP percent
stat.df.neg <- merge(df.negative, sample.meta.neg, by.x = "sample", by.y = "name") %>%
    dplyr::group_by(cutoff, lenti_plasmid, cell_line, system, is_matched, rep, exp_id) %>%
    dplyr::summarise(mean.percent.neg = mean(percent))

df.merged <- merge(stat.df.pos, stat.df.neg,
    by = c("cutoff", "rep", "cell_line", "system")
)
df.merged <- df.merged %>% rename(activation = percent, error = mean.percent.neg)

# calculate average values across replicates
df.merged %>%
    dplyr::group_by(cell_line, system, cutoff) %>%
    dplyr::summarise(mean.activation = mean(activation), mean.error = mean(error)) %>%
    data.frame() %>%
    # Add extrapolated values
    # add_row(system = "A-to-G_BC", cell_line = c("T3", "T4", "T5"), mean.error = 0, mean.activation = 0) %>%
    # add_row(system = "gRNA_BC", cell_line = c("T3", "T4", "T5"), mean.error = 0, mean.activation = 0) %>%
    # add_row(system = "CRISPRa_BC", cell_line = c("T3", "T4", "T5"), mean.error = 0, mean.activation = 0) %>%
    # add_row(system = "A-to-G_BC", cell_line = c("T3", "T4", "T5"), mean.error = 1, mean.activation = 1) %>%
    # add_row(system = "gRNA_BC", cell_line = c("T3", "T4", "T5"), mean.error = 1, mean.activation = 1) %>%
    # add_row(system = "CRISPRa_BC", cell_line = c("T3", "T4", "T5"), mean.error = 1, mean.activation = 1) %>%
    mutate(col_str = paste(system, cell_line, sep = "")) %>%
    arrange(
        cell_line, system, -cutoff
    ) -> df.mean

# output to csv
message("Exported df.mean and df.merged to csv...")
write.csv(df.merged, quote = FALSE, row.names = FALSE, file = sprintf("./processed_csv/ABE_error_activation_each_rep_%s.csv", "03072022"))
write.csv(df.mean, quote = FALSE, row.names = FALSE, file = sprintf("./processed_csv/ABE_error_activation_mean_%s.csv", "03072022"))