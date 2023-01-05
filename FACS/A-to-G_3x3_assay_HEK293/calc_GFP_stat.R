library(flowWorkspace)

gs.trans <- load_gs("./gs_achieve/gs_trans_CA060", verbose = TRUE)

sample.names <- gs_pop_get_stats(gs.trans, nodes = c("root"))$sample
sample.stats <- gs_pop_get_stats(gs.trans, type = "percent", nodes = c("GFP")) %>% mutate(percent = percent * 100)
sample.meta <- pData(gs.trans)
stat.df <- merge(sample.stats, sample.meta, by.x = "sample", by.y = "name")

# mediam EGFP expr
GFP.MFI <- gs_pop_get_stats(gs.trans, nodes = c("GFP"), type = pop.MFI, inverse.transform = TRUE)

df <- merge(GFP.MFI, stat.df, by.x = "sample", by.y = "sample")

# export data
write.csv(x = df, file = "./processed_csv/GFP_stat_ABE3x3_04112022.csv", quote = FALSE, row.names = FALSE)