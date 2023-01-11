library(ggplot2)
library(dplyr)
library(tidyr)

df <- read.table("data_summary_01302020.tsv", header = TRUE, sep = "\t")
df <- df %>% mutate(gRNA_type = case_when(Cell_line == gRNA ~ "OT", Cell_line != gRNA ~ "NT"))

df %>%
    group_by(gRNA_type) %>%
    summarise(value = list(EGFP_int)) %>%
    spread(gRNA_type, value) %>%
    mutate(p_value = t.test(unlist(OT), unlist(NT), var.equal = FALSE)$p.value)