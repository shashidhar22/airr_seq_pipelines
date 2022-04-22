#' Generate the cohort network figure. 
#' @description
#' Generate the nextwork plots comparing the control and PLHIV cohort
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
library(circlize)
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

comp_palette <- c("Post-ART" = "#f1a340", "Pre-ART" = "#998ec3", "Normal" = "#000000")
comp_one <- list(c("Pre-ART", "Post-ART"), c("Pre-ART", "Normal"), c("Post-ART", "Normal"))
comp_two <- list(c("Pre-ART", "Post-ART"))

hiv_matrix <- amino_table %>% 
  filter(str_detect(repertoire_id, "CFAR")) %>% 
  scoringMatrix(mode = "Similarity")
norm_matrix <- amino_table %>% 
  filter(!str_detect(repertoire_id, "CFAR")) %>% 
  scoringMatrix(mode = "Similarity")
hiv_break_names <- study_sample_table %>% 
  filter(repertoire_id %in% colnames(hiv_matrix)) %>% 
  arrange(patientID, years_rel_to_art) %>% 
  pull(patientID) 
hiv_break_points <- study_sample_table %>% 
  filter(repertoire_id %in% colnames(hiv_matrix)) %>% 
  arrange(patientID, years_rel_to_art) %>% 
  pull(repertoire_id) %>% 
  as.character()
hiv_row_names <- study_sample_table %>% 
  filter(repertoire_id %in% colnames(hiv_matrix)) %>% 
  arrange(patientID, years_rel_to_art) %>% 
  pull(patientID)
names(hiv_row_names) <- hiv_break_points
names(hiv_break_names) <- hiv_break_points
row_orders <- study_sample_table %>% 
  filter(repertoire_id %in% colnames(hiv_matrix)) %>% 
  arrange(patientID, years_rel_to_art) %>% 
  pull(repertoire_id)
hiv_map <- Heatmap(hiv_matrix, cluster_rows = FALSE, cluster_columns = FALSE, row_title = NULL, column_title = NULL, show_row_names = FALSE, row_order = hiv_break_points, column_order = hiv_break_points, show_column_names = FALSE, name = "Similarity\nScore", col = colorRamp2(c(0, 0.2, 1), c("#f0f9e8", "#7bccc4", "#0868ac")))
pdf(file = "Figure5A_CFAR_PLHIV.pdf", width = 6, height = 6)
hiv_map
dev.off()

norm_break_names <- study_sample_table %>% 
  filter(repertoire_id %in% colnames(norm_matrix)) %>% 
  arrange(patientID, years_rel_to_art) %>% 
  pull(patientID) 
norm_break_points <- study_sample_table %>% 
  filter(repertoire_id %in% colnames(norm_matrix)) %>% 
  arrange(patientID, years_rel_to_art) %>% 
  pull(repertoire_id)
norm_row_names <- study_sample_table %>% 
  filter(repertoire_id %in% colnames(norm_matrix)) %>% 
  arrange(patientID, years_rel_to_art) %>% 
  pull(time_point)
names(norm_row_names) <-norm_break_points
names(norm_break_names) <- norm_break_points
norm_map <- Heatmap(norm_matrix, cluster_rows = FALSE, cluster_columns = FALSE, row_title = NULL, column_title = NULL, show_row_names = FALSE, show_column_names = FALSE, name = "Similarity\nScores", col = colorRamp2(c(0, 0.2, 1), c("#f0f9e8", "#7bccc4", "#0868ac")))
pdf(file = "Figure5B_CFAR_Control.pdf", width = 6, height = 6)
norm_map
dev.off()

hiv_names <- rownames(hiv_matrix) 
samples_names <- study_sample_table %>% 
  pull(repertoire_id) 
patient_names <- study_sample_table %>% 
  pull(patientID)
names(patient_names) <- samples_names
hiv_table <- hiv_matrix %>% 
  as_tibble() %>% 
  mutate(rep_one = hiv_names) %>% 
  pivot_longer(-rep_one, names_to = "rep_two", values_to = "similarity_score") %>% 
  filter(rep_one != rep_two & patient_names[rep_one] == patient_names[rep_two]) %>% 
  mutate(group = "Within Individual",
         names = if_else(rep_one < rep_two, paste(rep_one, rep_two, sep =  "_"), paste(rep_two, rep_one, sep = "_")),
         pat_one = patient_names[rep_one],
         pat_two = patient_names[rep_two]) %>% 
  distinct(names, .keep_all = TRUE)

norm_names <- rownames(norm_matrix) 
norm_table <- norm_matrix %>%
  as_tibble() %>% 
  mutate(rep_one = norm_names) %>% 
  pivot_longer(-rep_one, names_to = "rep_two", values_to = "similarity_score") %>% 
  filter(rep_one != rep_two ) %>% 
  mutate(group = "Normal",
         names = if_else(rep_one < rep_two, paste(rep_one, rep_two, sep =  "_"), paste(rep_two, rep_one, sep = "_")),
         pat_one = patient_names[rep_one],
         pat_two = patient_names[rep_two]) %>% 
  distinct(names, .keep_all = TRUE)


across_table <- hiv_matrix %>% 
  as_tibble() %>% 
  mutate(rep_one = hiv_names) %>% 
  pivot_longer(-rep_one, names_to = "rep_two", values_to = "similarity_score") %>% 
  filter(rep_one != rep_two & patient_names[rep_one] != patient_names[rep_two]) %>% 
  mutate(group = "Across Individual",
         names = if_else(rep_one < rep_two, paste(rep_one, rep_two, sep =  "_"), paste(rep_two, rep_one, sep = "_")),
         pat_one = patient_names[rep_one],
         pat_two = patient_names[rep_two]) %>% 
  distinct(names, .keep_all = TRUE)

comp_three <- list(c("Within Individual", "Across Individual"), c("Across Individual", "Normal"), c("Within Individual", "Normal"))
comp_palette_three <- c("Within Individual" = "#39b6b9", "Across Individual" = "#89c348", "Normal" = "#000000")
sim_table <- bind_rows(hiv_table, norm_table, across_table)
overlap_map <- sim_table %>% 
    mutate(group = as_factor(group),
           group = fct_relevel(group, c("Within Individual", "Across Individual", "Normal"))) %>% 
    ggboxplot( x = "group", y = "similarity_score", color = "group", palette = comp_palette_three, add = "jitter") + 
    stat_compare_means(comparisons = comp_three) +
    labs(color = "Time point", x = "Time point", y = "Similarity") +
    theme_classic(base_size = 20)  

ggsave("Figure5C_CFAR_similarity.pdf", overlap_map, device = "pdf", width = 12,
    height = 4, units = "in", limitsize = FALSE)

