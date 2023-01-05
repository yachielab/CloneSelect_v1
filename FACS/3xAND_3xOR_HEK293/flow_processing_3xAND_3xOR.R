library(ggplot2)
library(scales)
library(dplyr)
library(ggforce)

FONT.SIZE <- 7
LINE.W <- 0.232

## Loading ABE001 data
df1 <- read.csv("exported_data/ABE001_exported/119_01.fcs_spill_applied_P1.txt", sep = "\t", header = TRUE, skip = 1)
df2 <- read.csv("exported_data/ABE001_exported/119_02.fcs_spill_applied_P1.txt", sep = "\t", header = TRUE, skip = 1)
df3 <- read.csv("exported_data/ABE001_exported/119_03.fcs_spill_applied_P1.txt", sep = "\t", header = TRUE, skip = 1)
df4 <- read.csv("exported_data/ABE001_exported/119_04.fcs_spill_applied_P1.txt", sep = "\t", header = TRUE, skip = 1)
df5 <- read.csv("exported_data/ABE001_exported/119_05.fcs_spill_applied_P1.txt", sep = "\t", header = TRUE, skip = 1)
df6 <- read.csv("exported_data/ABE001_exported/119_06.fcs_spill_applied_P1.txt", sep = "\t", header = TRUE, skip = 1)
df7 <- read.csv("exported_data/ABE001_exported/119_07.fcs_spill_applied_P1.txt", sep = "\t", header = TRUE, skip = 1)
df1$sample <- "1"
df2$sample <- "2"
df3$sample <- "3"
df4$sample <- "1+2"
df5$sample <- "2+3"
df6$sample <- "1+3"
df7$sample <- "1+2+3"
df1$order <- 1
df2$order <- 2
df3$order <- 3
df4$order <- 4
df5$order <- 5
df6$order <- 6
df7$order <- 7
df.ABE <- rbind(df1, df2, df3, df4, df5, df6, df7)

## ABE 3-AND
## fold-change of %EGFP cells
.lab <- c("1", "2", "3", "1+2", "2+3", "1+3", "1+2+3")
.order <- c(1, 2, 3, 4, 5, 6, 7)
.fc.val <- c(0, 0, 0, 0, 1, 1, 45)
.gfp.percent <- c(0.00, 0.00, 0.00, 0.00, 0.01, 0.01, 0.45)
fc.ABE <- data.frame(.lab, .fc.val, .gfp.percent, .order)

## Loading CA032 (AID) data
df1 <- read.csv("exported_data/CA032_exported/1.fcs_file_internal_comp_Ungated.txt", sep = "\t", header = TRUE, skip = 1)
df2 <- read.csv("exported_data/CA032_exported/2.fcs_file_internal_comp_Ungated.txt", sep = "\t", header = TRUE, skip = 1)
df3 <- read.csv("exported_data/CA032_exported/3.fcs_file_internal_comp_Ungated.txt", sep = "\t", header = TRUE, skip = 1)
df4 <- read.csv("exported_data/CA032_exported/4.fcs_file_internal_comp_Ungated.txt", sep = "\t", header = TRUE, skip = 1)
df5 <- read.csv("exported_data/CA032_exported/5.fcs_file_internal_comp_Ungated.txt", sep = "\t", header = TRUE, skip = 1)
df6 <- read.csv("exported_data/CA032_exported/6.fcs_file_internal_comp_Ungated.txt", sep = "\t", header = TRUE, skip = 1)
df7 <- read.csv("exported_data/CA032_exported/7.fcs_file_internal_comp_Ungated.txt", sep = "\t", header = TRUE, skip = 1)
df1$sample <- "1"
df2$sample <- "2"
df3$sample <- "3"
df4$sample <- "1+2"
df5$sample <- "2+3"
df6$sample <- "1+3"
df7$sample <- "1+2+3"
df1$order <- 1
df2$order <- 2
df3$order <- 3
df4$order <- 4
df5$order <- 5
df6$order <- 6
df7$order <- 7
df.AID <- rbind(df1, df2, df3, df4, df5, df6, df7)

## AID 3-OR
## fold-change of %EGFP cells
.fc.val <- c(0, 22, 247, 12, 121, 188, 96)
.gfp.percent <- c(0.00, 0.18, 2.42, 0.11, 1.20, 1.84, 0.92)
fc.AID <- data.frame(.lab, .fc.val, .gfp.percent, .order)

stop("ok")

## Fold-change bar plot for ABE AND gate
g.fc.plot.ABE <- ggplot(fc.ABE, aes(x = reorder(.lab, .order), y = .fc.val)) +
    geom_bar(stat = "identity", fill = "gray80", position = position_dodge(width = 0.8)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 50), breaks = c(0, 25, 50)) +
    theme(
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE)
    ) +
    labs(y = "", x = "")
# ggsave(plot=g.fc.plot.ABE, filename='./plots/ABE_3AND_fc_barplot.pdf', w=3, h=1)

g.percent.plot.ABE <- ggplot(fc.ABE, aes(x = reorder(.lab, .order), y = .gfp.percent)) +
    geom_bar(stat = "identity", fill = "green3", position = position_dodge(width = 0.8)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0, 0.5, 1.0)) +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE)
    ) +
    labs(y = "%EGFP\npositive cells", x = "")
ggsave(plot = g.percent.plot.ABE, filename = "./plots/ABE_3AND_GFP_barplot.pdf", w = 2.65, h = 1)

## Fold-change bar plot for Target-AID OR gate
g.fc.plot.AID <- ggplot(fc.AID, aes(x = reorder(.lab, .order), y = .fc.val)) +
    geom_bar(stat = "identity", fill = "gray80", position = position_dodge(width = 0.8)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 300), breaks = c(0, 100, 200, 300)) +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE)
    ) +
    labs(y = "", x = "")
# ggsave(plot=g.fc.plot.AID, filename='./plots/AID_3OR_fc_barplot.pdf', w=3, h=1)

g.percent.plot.AID <- ggplot(fc.AID, aes(x = reorder(.lab, .order), y = .gfp.percent)) +
    geom_bar(stat = "identity", fill = "green3", position = position_dodge(width = 0.8)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 5), breaks = c(0, 2.5, 5.0)) +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE)
    ) +
    labs(x = "", y = "%EGFP\npositive cells")
ggsave(plot = g.percent.plot.AID, filename = "./plots/AID_3OR_GFP_barplot.pdf", w = 2.65, h = 1)

df.ABE <- df.ABE %>% filter(FITC.A > 1)
col <- ifelse(df.ABE$FITC.A > 91.61, "green3", "gray80") ## Need to be checked
g.jitter.ABE <- ggplot(df.ABE %>% filter(FITC.A > 1), aes(y = FITC.A, x = reorder(sample, order))) +
    # geom_jitter(colour='green3', position=position_jitter(0.4), size=0.5) +
    geom_violin(size = 0.2, colour = "gray50", trim = FALSE) +
    geom_sina(colour = col, shape = 16, size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = 91.61, colour = "gray40", linetype = "dashed", size = 0.2) +
    scale_y_log10(expand = c(0, 0), limits = c(1, 10000), labels = trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() +
    theme(
        ## plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = FONT.SIZE)
    ) +
    labs(x = "", y = "EGFP intensity (a.u.)")
ggsave(plot = g.jitter.ABE, filename = "./plots/ABE_3AND_sina.pdf", w = 2.55, h = 1.7)

df.AID <- df.AID %>% filter(FITC.A > 1)
col <- ifelse(df.AID$FITC.A > 91.61, "green3", "gray80") ## Need to be checked
g.jitter.AID <- ggplot(df.AID %>% filter(FITC.A > 1), aes(y = FITC.A, x = reorder(sample, order))) +
    geom_violin(size = 0.2, colour = "gray50", trim = FALSE) +
    geom_sina(colour = col, shape = 16, size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = 91.61, colour = "gray40", linetype = "dashed", size = 0.2) +
    # geom_jitter(colour='green3', position=position_jitter(0.4), size=0.5) +
    ## geom_crossbar(xintercept=85.39, width=0.1) +
    # scale_y_continuous(limits=c(0, 3000)) +
    scale_y_log10(expand = c(0, 0), limits = c(1, 10000), labels = trans_format("log10", scales::math_format(10^.x))) +
    theme_classic() +
    theme(
        ## plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(colour = "black", size = LINE.W),
        axis.line = element_line(size = LINE.W),
        axis.title = element_text(size = FONT.SIZE),
        axis.text = element_text(colour = "black", size = FONT.SIZE),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = FONT.SIZE)
    ) +
    labs(x = "", y = "EGFP intensity (a.u.)")
ggsave(plot = g.jitter.AID, filename = "./plots/AID_3OR_sina.pdf", w = 2.55, h = 1.7)