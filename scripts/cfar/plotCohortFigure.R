#' Generate the cohort summary figure. 
#' @description
#' Generate the dot plots showing the viral load, CD4 counts and repertoire
#' clonality across time
#' @return Cohort summary figure
library(tidyverse)
library(LymphoSeq2)
library(patchwork)
library(tidygraph)
library(lubridate)
library(ggraph)
library(ggpubr)
library(ggalluvial)
library(readxl)
library(knitr)
library(scales)
library(optparse)
library(DT)
library(ComplexHeatmap)
set.seed(12357) # Set the seed for reproducibility
option_list <- list(
    make_option(c("-a", "--amino_table"),
        type = "character",
        help = "Path to the sampled amino acid table",
        dest = "amino_table"),
    make_option(c("-m", "--meta_table"),
        type = "character",
        help = "Path to the study meta data in RDA format",
        dest = "meta_table"),
    make_option(c("-r", "--repertoire_summary_table"),
        type = "character",
        help = "Path to the repertoire summary tables",
        dest = "rep_table")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list, add_help_option=FALSE), 
    print_help_and_exit = TRUE) # Parse the arguments

load(parser$amino_table)
load(parser$meta_table)
summary_table <- read_tsv(parser$rep_table)

comp_palette <- c("Post-ART" = "#f1a340", "Pre-ART" = "#998ec3", "Normal" = "#000000")
comp_one <- list(c("Pre-ART", "Post-ART"), c("Pre-ART", "Normal"), c("Post-ART", "Normal"))
comp_two <- list(c("Pre-ART", "Post-ART"))

cc_table <- study_cd4_count %>% 
  rename(cd4_count = result) %>% 
  mutate(resultsRel = cd4_count, 
    patientID = forcats::fct_rev(patientID),
    yearsRel = years_rel_to_art,
    time_point = dplyr::if_else(years_rel_to_art <= 0, "Pre-ART", "Post-ART"))

## Calculate days relative ART treatment for ImmunoSeq and Viral load table
vl_table <- study_viral_load %>%
  rename(viral_load = result) %>% 
  dplyr::mutate(patientID = forcats::fct_rev(patientID),
  yearsRel = years_rel_to_art,
  time_point = dplyr::if_else(years_rel_to_art <= 0, "Pre-ART", "Post-ART"))


itable <- summary_table %>%
  filter(repertoire_id %in% amino_table$repertoire_id & str_detect(repertoire_id, "CFAR")) %>%
  left_join(study_sample_table, by = "repertoire_id") %>% 
  left_join(study_patient_table, by = "patientID") %>% 
  dplyr::mutate(patientID = forcats::fct_rev(patientID),
    yearsRel = years_rel_to_art, 
    time_point = dplyr::if_else(years_rel_to_art <= 0, "Pre-ART", "Post-ART"))

#Calculate ranking

rank_table <- vl_table %>% 
  filter(time_point == "Post-ART") %>%
  mutate(vl_supression = if_else(viral_load < 200, TRUE, FALSE)) %>% 
  group_by(patientID, vl_supression) %>% 
  summarise(count = n()) %>% 
  pivot_wider(names_from = vl_supression, values_from = count, values_fill = 0) %>% 
  mutate(supression_rate = `TRUE`/(`FALSE` + `TRUE`)) %>% 
  select(supression_rate)

# Setting min and max date range
date_lim <- c(-12.0, 12.0)      

# Plot viral load
plot_vl <- vl_table %>%
  left_join(rank_table, by = "patientID") %>%
  mutate(patientID = forcats::fct_reorder(patientID, supression_rate)) %>%
  ggplot(aes(x = yearsRel, y = patientID)) +
  geom_point(aes(color = time_point, size = log10(viral_load ))) +
  geom_segment(x=-12, xend=12, aes(y=patientID, yend=patientID)) +
  theme_classic(base_size = 20) +
  xlab("Viral load (cells/mL)\nDays relative to ART initiation") + 
  ylab("Patient ID") +
  labs(size = "Log10(VL)") +
  scale_color_manual(values = c("Pre-ART" = "#998ec3", "Post-ART" = "#f1a340")) +
  scale_x_continuous(position = "top", limits = date_lim, breaks = seq(-12, 12, 4)) + 
  guides(color = "none") + 
  scale_size_binned(range = c(2.0, 12.0), trans = "exp") +
  theme(legend.position="bottom", 
    legend.box="vertical", 
    legend.margin=margin(),
    legend.key.width = unit(2, "cm"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text( face = "bold"))

# Plot CD4 count variation
plot_cc <- cc_table %>% 
  left_join(rank_table, by = "patientID") %>%
  mutate(patientID = forcats::fct_reorder(patientID, supression_rate)) %>%
  ggplot(aes(x = yearsRel, y = patientID)) +
  geom_point(aes(color = time_point, size = resultsRel)) +
  geom_segment(x=-12, xend=12, aes(y=patientID, yend=patientID)) +
  theme_classic(base_size = 20) +
  xlab("CD4 count (cell/uL)\nDays relative to ART initiation") + 
  ylab("Patient ID") +
  labs(color = "Treatment status", size = "CD4 Count") +
  scale_x_continuous(position = "top", limits = date_lim, breaks = seq(-12, 12, 4)) + 
  scale_size_binned_area(max_size = 12) +
  scale_color_manual(values = c("Pre-ART" = "#998ec3", "Post-ART" = "#f1a340")) +
  guides(color = "none") +
  theme(legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank())

# Plot clonality variation 
plot_cl<-  itable %>%
  left_join(rank_table, by = "patientID") %>%
  mutate(patientID = forcats::fct_reorder(patientID, supression_rate)) %>%
  ggplot(aes(x = yearsRel, y = patientID)) +
  geom_point(aes(color = time_point, size = clonality)) +
  geom_segment(x=-12, xend=12, aes(y=patientID, yend=patientID)) +
  theme_classic(base_size = 20) +
  xlab("T-cell repertoire sequencing\nDays relative to ART initiations") + 
  ylab("Patient ID") +
  labs(color = "Treatment status", size = "Clonality of \n TCR repertoire") +
  scale_x_continuous(position = "top", limits = date_lim, breaks = seq(-12, 12, 4)) + 
  scale_size_binned_area(max_size = 12) +
  scale_color_manual(values = c("Pre-ART" = "#998ec3", "Post-ART" = "#f1a340")) +
  guides(color = "none") +
  theme(legend.position = "bottom",
    legend.box="vertical", 
    legend.margin=margin(),
    legend.key.width = unit(2, "cm"),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank())

# Compose plots
figure_one <- patchwork::wrap_plots(plot_vl, 
    plot_cc, 
    plot_cl, 
    nrow = 1, 
    widths = c(8,8,8)) + 
  patchwork::plot_layout()

ggsave("Figure1_CFAR_summary.pdf", figure_one, device = "pdf", width = 24,
    height = 14, units = "in", limitsize = FALSE)