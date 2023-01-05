library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scales)
# library(jcolors)

FONT.SIZE <- 7
FONT.SIZE.TITLE <- 8
LINE.W <- 0.232 # 0.5pt in Keynote

# Load data
CBE.df <- read.csv("data/6x6_error_activation_mean_04132022.tsv", sep = ",", header = TRUE)
ABE.df <- read.csv("data/ABE_error_activation_mean_03072022.csv", sep = ",", header = TRUE)

# Use same col names
colnames(CBE.df) <- c("cell_line", "system", "cutoff", "mean.activation", "mean.error", "col_str")

# Merge them
merged.df <- rbind(CBE.df, ABE.df)

# Sanitize different label name in two datasets
clean.df <- merged.df %>%
    mutate(system = replace(system, system == "C-to-T", "C-to-T")) %>%
    mutate(system = replace(system, system == "A-to-G_BC", "A-to-G")) %>%
    mutate(system = replace(system, system == "CRISPRa", "High-copy CRISPRa")) %>%
    mutate(system = replace(system, system == "gRNA_BC", "High-copy CRISPRa")) %>%
    mutate(system = replace(system, system == "CRISPRa_BC", "Low-copy CRISPRa")) %>%
    mutate(system = replace(system, system == "gRNA", "Low-copy CRISPRa")) %>%
    mutate(
        bc_type =
            case_when(
                cell_line == "T3" ~ "A-to-G BC",
                cell_line == "T4" ~ "A-to-G BC",
                cell_line == "T5" ~ "A-to-G BC",
                TRUE ~ "C-to-T BC"
            )
    )

clean.df$system <- factor(clean.df$system, levels = c("C-to-T", "A-to-G", "Low-copy CRISPRa", "High-copy CRISPRa"))


# get some stats
clean.df %>%
    filter(system == "A-to-G") %>%
    mutate(mean.error = mean.error * 100, mean.activation = mean.activation * 100) %>%
    filter(mean.error <= 0.5) %>%
    group_by(cell_line) %>%
    summarise(max = max(mean.activation)) %>%
    data.frame()

stop("ok...")


plt.ROC <- function(df) {
    g <- ggplot(df, aes(y = mean.activation, x = mean.error)) +
        geom_abline(intercept = 0, linetype = "dashed", colour = "black", size = LINE.W) +
        geom_line(aes(group = col_str, colour = bc_type), size = 0.5) +
        coord_cartesian(xlim = c(0, 0.01), ylim = c(0, 0.4)) +
        scale_colour_manual(values = c("C-to-T BC" = "#E05A99", "A-to-G BC" = "#87A5E0")) +
        scale_y_continuous(
            expand = c(0, 0),
            limits = c(0, 0.4),
            label = scales::percent_format(suffix = "")
        ) +
        scale_x_continuous(
            expand = c(0.005, 0),
            limits = c(0, 0.1),
            breaks = c(0, 0.005, 0.01),
            label = scales::percent_format(suffix = "")
        ) +
        labs(x = "", y = "") +
        theme_classic() +
        theme(
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "lines"),
            panel.background = element_rect(fill = "transparent", colour = "transparent"),
            plot.background = element_rect(fill = "transparent", , colour = "transparent"),
            plot.title = element_text(colour = "black", size = FONT.SIZE),
            strip.background = element_rect(colour = NA, fill = "gray90"),
            panel.border = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_line(colour = "black", size = LINE.W),
            axis.line = element_line(colour = "black", size = LINE.W),
            axis.text = element_text(colour = "black", size = FONT.SIZE)
        ) +
        guides(colour = "none")
    return(g)
}

# Iterate over activation systems
g1 <- plt.ROC(clean.df %>% filter(system == "C-to-T"))
g2 <- plt.ROC(clean.df %>% filter(system == "A-to-G"))
g3 <- plt.ROC(clean.df %>% filter(system == "Low-copy CRISPRa", bc_type == "C-to-T BC"))
g4 <- plt.ROC(clean.df %>% filter(system == "Low-copy CRISPRa", bc_type == "A-to-G BC"))

ggsave(plot = g1, filename = "./plots/08192022_ROC_C-to-T.pdf", w = 1.3, h = 1.3, device = cairo_pdf)
ggsave(plot = g2, filename = "./plots/08192022_ROC_A-to-G.pdf", w = 1.3, h = 1.3, device = cairo_pdf)
ggsave(plot = g3, filename = "./plots/08192022_ROC_C-to-T_CaLow.pdf", w = 1.3, h = 1.3, device = cairo_pdf)
ggsave(plot = g4, filename = "./plots/08192022_ROC_A-to-G_CaLow.pdf", w = 1.3, h = 1.3, device = cairo_pdf)