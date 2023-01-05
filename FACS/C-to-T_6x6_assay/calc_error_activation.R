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

compute_dynamic_gate <- function(gs, inv.bins, parent_gate = "P2") {
    df.stat <- data.frame()
    i <- 1
    for (min.cutoff in inv.bins) {
        message(sprintf("Processing of %f, %d of %d...", min.cutoff, i, length(inv.bins)))

        # Create dynamic FL1-A gate
        rg_gfp <- flowCore::rectangleGate(
            "FL1-A" = c(min.cutoff, Inf), "FSC-A" = c(-Inf, Inf),
            filterId = "GFP_dynamic"
        )
        gates <- flowWorkspace::gs_get_pop_paths(gs)
        if ("/P1/P2/GFP_dynamic" %in% gates) {
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

# load pre-processed GatingSet object
gs.trans <- load_gs("./CS062_gs", verbose = TRUE)

gs.positive <- subset(gs.trans, is_match == "OT")
gs.negative <- subset(gs.trans, is_match == "NT")

bins <- seq(0, 1e5, by = 100) # by=50
bins.trans <- log.trans(bins)

# apply dynamic gates
df.positive <- compute_dynamic_gate(gs.positive, bins.trans[1:1000]) # [1:500]
df.negative <- compute_dynamic_gate(gs.negative, bins.trans[1:1000])

# merge data.frame
sample.meta.pos <- pData(gs.positive)
sample.meta.neg <- pData(gs.negative)
stat.df.pos <- merge(df.positive, sample.meta.pos, by.x = "sample", by.y = "name")

# Group two expected negative cell lines and calculate mean GFP percent
stat.df.neg <- merge(df.negative, sample.meta.neg, by.x = "sample", by.y = "name") %>%
    dplyr::group_by(is_match, Rep, System, cutoff, Cell_line, Cell_line2, Assay_group) %>%
    dplyr::summarize(mean.percent.neg = mean(percent))

df.merged <- merge(stat.df.pos, stat.df.neg,
    by = c("cutoff", "Rep", "Cell_line", "System", "Cell_line2", "Assay_group")
)
df.merged <- df.merged %>% rename(activation = percent, error = mean.percent.neg)

# calculate average values across replicates
BCs <- c("BC2", "BC4", "BC6", "BC8", "BC10", "BC12")
df.mean <- df.merged %>%
    dplyr::group_by(Cell_line, System, cutoff) %>%
    dplyr::summarise(mean.activation = mean(activation), mean.error = mean(error)) %>%
    # data.frame() %>%
    # Add extrapolated (0, 0) values
    # add_row(System = "C-to-T", Cell_line = BCs, mean.error = 0, mean.activation = 0) %>%
    # add_row(System = "gRNA", Cell_line = BCs, mean.error = 0, mean.activation = 0) %>%
    # add_row(System = "CRISPRa", Cell_line = BCs, mean.error = 0, mean.activation = 0) %>%
    # Add extrapolated (1, 1) values
    # add_row(System = "C-to-T", Cell_line = BCs, mean.error = 1, mean.activation = 1) %>%
    # add_row(System = "gRNA", Cell_line = BCs, mean.error = 1, mean.activation = 1) %>%
    # add_row(System = "CRISPRa", Cell_line = BCs, mean.error = 1, mean.activation = 1) %>%
    mutate(col_str = paste(System, Cell_line, sep = "")) %>%
    arrange(
        Cell_line, System, -cutoff
    )

write.csv(df.mean, quote = FALSE, row.names = FALSE, file = sprintf("6x6_error_activation_mean_%s.csv", "04132022"))