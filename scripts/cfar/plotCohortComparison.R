#' Generate the cohort comparison figure. 
#' @description
#' Generate the box plots comparing the control and PLHIV cohort
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
data_summary <- read_tsv(parser$rep_table)

readAverageMetrics <- function(rda_path) {
  load(rda_path)
  return(summary_table)
}

readTranslationTable <- function(table_path) {
  load(table_path)
  return(translation_efficacy)
} 

averaged_metrics <- list.files(path = "./", 
    pattern = "summary.rda", 
    full.names = TRUE) %>% 
  map(readAverageMetrics) %>% 
  bind_rows() %>% 
  select(-simpson_index, -inverse_simpson, -chao_estimate, -hill_estimate,
    -kemp_estimate) %>% 
  group_by(repertoire_id) %>% 
  summarize(total_sequences = mean(total_sequences),
    unique_productive_sequences = mean(unique_productive_sequences),
    total_count = mean(total_count),
    clonality = mean(clonality),
    gini_coefficient = mean(gini_coefficient),
    top_productive_sequence = mean(top_productive_sequence),
    translation_efficacy = mean(translation_efficacy)) %>% 
  left_join(study_sample_table, by = "repertoire_id") %>% 
  select(patientID, time_point, repertoire_id, everything())

comp_palette <- c("Post-ART" = "#f1a340", "Pre-ART" = "#998ec3", "Normal" = "#000000")
comp_one <- list(c("Pre-ART", "Post-ART"), c("Pre-ART", "Normal"), c("Post-ART", "Normal"))
comp_two <- list(c("Pre-ART", "Post-ART"))

# Modify the study viral load table
study_viral_load <- study_viral_load %>% 
  select(patientID, years_rel_to_art, result) %>% 
  rename(viral_load = result) %>%
  left_join(study_patient_table, by = "patientID") %>% 
  mutate(time_from_art = years_rel_to_art,
    time_point = if_else(time_from_art <= 0, "Pre-ART", "Post-ART"))
# Modify the study CD4 table
study_cd4_count <- study_cd4_count %>% 
  select(patientID, years_rel_to_art, result) %>% 
  rename(cd4_count = result) %>%
  left_join(study_patient_table, by = "patientID") %>% 
  mutate(time_from_art = years_rel_to_art,
    time_point = if_else(time_from_art <= 0, "Pre-ART", "Post-ART"))

viral_load <- study_viral_load %>% 
  mutate(time_point = as_factor(time_point),
    time_point = fct_relevel(time_point, c("Pre-ART", "Post-ART"))) %>%
  ggboxplot(x = "time_point", y = "viral_load", color = "time_point", palette = comp_palette, add = "jitter") +
  stat_compare_means(comparisons = comp_two) +
  theme_classic(base_size = 20) +
  scale_y_log10() + 
  labs(color = "Time point", x = "Time point", y = "Log10(Viral load \n (Copies/mL))") + 
  annotation_logticks() 

cdfour_count <- study_cd4_count %>% 
  mutate(time_point = as_factor(time_point),
    time_point = fct_relevel(time_point, c("Pre-ART", "Post-ART"))) %>%
  ggboxplot(x = "time_point", y = "cd4_count", color = "time_point", palette = comp_palette, add = "jitter") +
  stat_compare_means(comparisons = comp_two) +
  theme_classic(base_size = 20) +
  labs(color = "Time point", x = "Time point", y = "CD4 Count \n (Cells/uL)") 

total_sequences <- ggboxplot(averaged_metrics, x = "time_point", y = "total_sequences", color = "time_point", palette = comp_palette, add = "jitter") + 
  stat_compare_means(comparisons = comp_one) +
  theme_classic(base_size = 20) +
  ylim(0, 6000) +
  labs(color = "Time point", x = "Time point", y = "Total sequences") +
  theme(axis.title.x = element_blank(),
    axis.text.x = element_blank())

unique_productive <- ggboxplot(averaged_metrics, x = "time_point", y = "unique_productive_sequences", color = "time_point", palette = comp_palette, add = "jitter") + 
  stat_compare_means(comparisons = comp_one) +
  theme_classic(base_size = 20) +
  ylim(0, 6000) +
  labs(color = "Time point", x = "Time point", y = "Unique productive \n sequences") +
  theme(axis.title.x = element_blank(),
    axis.text.x = element_blank())

clonality <- ggboxplot(averaged_metrics, x = "time_point", y = "clonality", color = "time_point", palette = comp_palette, add = "jitter") + 
  stat_compare_means(comparisons = comp_one) +
  theme_classic(base_size = 20) +
  ylim(0, 1.0) +
  labs(color = "Time point", x = "Time point", y = "Clonality") +
  theme(axis.title.x = element_blank(),
    axis.text.x = element_blank())

gini_coefficient <- ggboxplot(averaged_metrics, x = "time_point", y = "gini_coefficient", color = "time_point", palette = comp_palette, add = "jitter") + 
  stat_compare_means(comparisons = comp_one) +
  theme_classic(base_size = 20) +
  labs(color = "Time point", x = "Time point", y = "Gini coefficient") +
  theme(axis.title.x = element_blank(),
    axis.text.x = element_blank())

top_productive <- ggboxplot(averaged_metrics, x = "time_point", y = "top_productive_sequence", color = "time_point", palette = comp_palette, add = "jitter") + 
  stat_compare_means(comparisons = comp_one) +
  theme_classic(base_size = 20) +
  labs(color = "Time point", x = "Time point", y = "Percentage of repertoire \n occupied by the top clone") 

te_table <- list.files("./", pattern = "translation_efficacy.rda", full.names = TRUE) %>% 
  map(readTranslationTable) %>% 
  bind_rows() %>% 
  left_join(study_sample_table, by = "repertoire_id") %>% 
  group_by(time_point, repertoire_id) %>% summarize(translation_efficacy = mean(translation_efficacy))

translation_efficacy <- ggboxplot(te_table, x = "time_point", y = "translation_efficacy", color = "time_point", palette = comp_palette, add = "jitter") + 
  stat_compare_means(comparisons = comp_one) +
  theme_classic(base_size = 20) +
  labs(color = "Time point", x = "Time point", y = "Translation efficacy") 

figure_two <- (viral_load + cdfour_count) / (total_sequences + unique_productive) / (clonality + gini_coefficient) / (top_productive + translation_efficacy) +
  patchwork::plot_layout(guides = 'collect')


ggsave("Figure2_CFAR_comparison.pdf", figure_two, device = "pdf", width = 14,
    height = 16, units = "in", limitsize = FALSE)