library(ggplot2)
library(scales)
library(dplyr)
library(ggridges)
library(stringr)

DATE <- Sys.Date()

breaks_log10 <- function(x) {
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    10^(seq.int(low, high))
}

breaks_5log10 <- function(x) {
    low <- floor(log10(min(x) / 5))
    high <- ceiling(log10(max(x) / 5))
    5 * 10^(seq.int(low, high))
}

FONT.SIZE <- 7
LABEL.FONT.SIZE <- 7
LINE.W <- 0.232 # 0.5pt
BC <- "V4-BC2"
Target_AID <- "A"
gRNA_BC <- "B"
CRISPR_BC <- "C"

generate.data.BC2 <- function(system, ...) {
    if (str_detect(system, "CtoT_BC")) {
        well_id <- "A"
    } else if (str_detect(system, "gRNA_BC")) {
        well_id <- "B"
    } else if (str_detect(system, "CRISPRa_BC")) {
        well_id <- "C"
    } else {
        stop("Error found:")
    }

    df_1 <- read.csv(sprintf("P1_exported_data/CA053_%s1.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_2 <- read.csv(sprintf("P1_exported_data/CA053_%s2.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_3 <- read.csv(sprintf("P1_exported_data/CA053_%s3.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_1$rep <- "Rep1"
    df_2$rep <- "Rep2"
    df_3$rep <- "Rep3"
    df_merged_on <- rbind(df_1, df_2, df_3)
    df_merged_on$cond <- "OT gRNA"

    df_4 <- read.csv(sprintf("P1_exported_data/CA053_%s4.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_5 <- read.csv(sprintf("P1_exported_data/CA053_%s5.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_6 <- read.csv(sprintf("P1_exported_data/CA053_%s6.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_4$rep <- "Rep1"
    df_5$rep <- "Rep2"
    df_6$rep <- "Rep3"
    df_merged_off <- rbind(df_4, df_5, df_6)
    df_merged_off$cond <- "NT gRNA"
    df_all <- rbind(df_merged_on, df_merged_off)
    df_all$system <- system
    return(df_all)
}

generate.data.BC4 <- function(system, ...) {
    if (str_detect(system, "CtoT_BC")) {
        well_id <- "A"
    } else if (str_detect(system, "gRNA_BC")) {
        well_id <- "B"
    } else if (str_detect(system, "CRISPRa_BC")) {
        well_id <- "C"
    } else {
        stop("Error found:")
    }

    df_1 <- read.csv(sprintf("P1_exported_data/CA053_%s7.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_2 <- read.csv(sprintf("P1_exported_data/CA053_%s8.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_3 <- read.csv(sprintf("P1_exported_data/CA053_%s9.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_1$rep <- "Rep1"
    df_2$rep <- "Rep2"
    df_3$rep <- "Rep3"
    df_merged_on <- rbind(df_1, df_2, df_3)
    df_merged_on$cond <- "NT gRNA"

    df_4 <- read.csv(sprintf("P1_exported_data/CA053_%s10.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_5 <- read.csv(sprintf("P1_exported_data/CA053_%s11.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_6 <- read.csv(sprintf("P1_exported_data/CA053_%s12.fcs_file_internal_comp_P1.txt", well_id), sep = "\t", skip = 1, header = TRUE)
    df_4$rep <- "Rep1"
    df_5$rep <- "Rep2"
    df_6$rep <- "Rep3"
    df_merged_off <- rbind(df_4, df_5, df_6)
    df_merged_off$cond <- "OT gRNA"
    df_all <- rbind(df_merged_on, df_merged_off)
    df_all$system <- system
    return(df_all)
}

plt.density <- function(x, ...) {
    g <- ggplot(x %>% filter(`B525.FITC.A_.FL1.A.` > 1.0), aes(y = reorder(rep, desc(rep)), x = `B525.FITC.A_.FL1.A.`, fill = stat(x))) +
        facet_wrap(~cond, ncol = 1) +
        geom_density_ridges_gradient(
            colour = "white", size = LINE.W, quantile_lines = TRUE, quantiles = 2,
            vline_size = 0.5, vline_color = "gray", scale = 2.0
        ) +
        scale_fill_gradientn(limits = c(1, 7), colours = c("black", "gray20", "green2", "greenyellow")) +
        scale_x_log10(
            breaks = c(10^1, 10^3, 10^5, 10^7),
            labels = scales::trans_format("log10", math_format(10^.x)),
            expand = c(0, 0),
            limits = c(10, 10000000)
        ) +
        scale_y_discrete(expand = c(0, 0)) +
        labs(x = "", y = "") +
        theme_classic() +
        theme(
            plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
            panel.spacing.y = unit(0.2, "lines"), # space between top and bottom
            panel.background = element_rect(fill = "transparent", colour = "transparent"),
            plot.background = element_rect(fill = "transparent", colour = "transparent"),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(colour = "black", size = LINE.W),
            axis.title = element_text(size = FONT.SIZE),
            axis.text.x = element_text(colour = "black", size = LABEL.FONT.SIZE),
            axis.text.y = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            strip.background = element_blank()
        ) +
        guides(fill = "none")
    return(g)
}

fig.size.w <- 1.0
fig.size.h <- 1.7

### BC2 data
df1 <- generate.data.BC2("CtoT_BC")
df2 <- generate.data.BC2("gRNA_BC")
df3 <- generate.data.BC2("CRISPRa_BC")
ggsave(plt.density(df1), filename = sprintf("./plt_histogram/BC2_CtoTBC_%s.pdf", "04052022"), w = fig.size.w, h = fig.size.h)
ggsave(plt.density(df2), filename = sprintf("./plt_histogram/BC2_gRNABC_%s.pdf", "04052022"), w = fig.size.w, h = fig.size.h)
ggsave(plt.density(df3), filename = sprintf("./plt_histogram/BC2_CRISPRaBC_%s.pdf", "04052022"), w = fig.size.w, h = fig.size.h)

### BC4 data
df4 <- generate.data.BC4("CtoT_BC")
df5 <- generate.data.BC4("gRNA_BC")
df6 <- generate.data.BC4("CRISPRa_BC")
ggsave(plt.density(df4), filename = sprintf("./plt_histogram/BC4_CtoTBC_%s.pdf", "04052022"), w = fig.size.w, h = fig.size.h)
ggsave(plt.density(df5), filename = sprintf("./plt_histogram/BC4_gRNABC_%s.pdf", "04052022"), w = fig.size.w, h = fig.size.h)
ggsave(plt.density(df6), filename = sprintf("./plt_histogram/BC4_CRISPRaBC_%s.pdf", "04052022"), w = fig.size.w, h = fig.size.h)