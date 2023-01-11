library(dplyr)
library(tidyr)

df <- read.table("11072021_CA053_summary.tsv", sep = "\t", header = TRUE)
df$System <- factor(df$System, levels = c("CtoT_BC", "CRISPRa_BC", "gRNA_BC"))

aid.mfi <- df %>% filter(System == "CtoT_BC", gRNA_type == "OT_gRNA")
grna.mfi <- df %>% filter(System == "gRNA_BC", gRNA_type == "OT_gRNA")
ca.mfi <- df %>% filter(System == "CRISPRa_BC", gRNA_type == "OT_gRNA")

# C-to-T vs gRNA BC, MFI
wilcox.test(aid.mfi$Median_GFP_exp, grna.mfi$Median_GFP_exp)

# C-to-T vs CRIPSRa BC, MFI
wilcox.test(aid.mfi$Median_GFP_exp, ca.mfi$Median_GFP_exp)

## comparison with OT vs NT
p.val.df <- df %>%
    filter(System != "ctrl") %>%
    group_by(System, gRNA_type) %>%
    summarise(value = list(Percent_FITC)) %>%
    spread(gRNA_type, value) %>%
    mutate(p_value = t.test(unlist(OT_gRNA), unlist(NT_gRNA), var.equal = FALSE)$p.value)