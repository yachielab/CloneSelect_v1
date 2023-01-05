library(flowWorkspace)
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
    df.stat <- data.frame()
    i <- 1
    for (min.cutoff in inv.bins) {
        message(sprintf("Processing of %f, %d of %d...", min.cutoff, i, length(inv.bins)))

        # Create dynamic FL1-A gate
        rg_gfp <- flowCore::rectangleGate(
            "FITC-A" = c(min.cutoff, Inf), "FSC-A" = c(-Inf, Inf),
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

# load pre-processed GatingSet object from local disk
gs.trans <- load_gs("./gs_trans", verbose = TRUE)

gs.positive <- subset(gs.trans, OT == "PC")
gs.negative <- subset(gs.trans, OT == "NC")

bins <- seq(-100, 1e4, by = 5) # FIXME: need to be fixed - negative GFP values should be considered...
bins.trans <- log.trans(bins)

# apply dynamic gates
df.positive <- compute_dynamic_gate(gs.positive, bins.trans[1:200]) # 1000
df.negative <- compute_dynamic_gate(gs.negative, bins.trans[1:200])

# merge data.frame
sample.meta.pos <- pData(gs.positive)
sample.meta.neg <- pData(gs.negative)
stat.df.pos <- merge(df.positive, sample.meta.pos, by.x = "sample", by.y = "name")
stat.df.neg <- merge(df.negative, sample.meta.neg, by.x = "sample", by.y = "name")

df.merged <- merge(stat.df.pos, stat.df.neg,
    by = c("cutoff", "dose_str", "barcode", "exp_id", "dose_ng", "dose", "pop")
)

df.merged <- df.merged %>%
    select(-c(sample.x, sample.y, OT.x, OT.y)) %>%
    rename(activation = percent.x, error = percent.y) %>%
    arrange(barcode, dose_str, -cutoff)

write.csv(df.merged, quote = FALSE, row.names = FALSE, file = sprintf("dose_error_activation_mean_%s.csv", "04192022"))s