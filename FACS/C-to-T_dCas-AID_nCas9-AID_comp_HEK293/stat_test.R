library(dplyr)
library(tidyr)

df <- read.csv("GFP_stats.csv")
pval <- df %>%
    filter(AID_type != "ctrl") %>%
    group_by(AID_type, Target_type) %>%
    summarise(value = list(percent)) %>%
    spread(Target_type, value) %>%
    mutate(p_value = t.test(unlist(OT), unlist(NT))$p.value)
print(pval)