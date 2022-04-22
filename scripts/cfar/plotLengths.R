#' Generate the cohort Length figure. 
#' @description
#' Generate the length plots comparing the control and PLHIV cohort
#' @return Cohort length figure
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

annotate_table <- amino_table %>% 
    left_join(study_sample_table, by = "repertoire_id") %>% 
    left_join(study_patient_table, by = "patientID") %>%
    mutate(lengths = nchar(junction_aa)) %>% 
    select(repertoire_id, time_group, time_point, junction_aa, lengths) 

length_dist <- annotate_table %>% 
    group_by(repertoire_id, lengths) %>% 
    summarise(count = n()) %>%
    ungroup() %>% 
    group_by(repertoire_id) %>%
    mutate(percent = (count *100)/sum(count))

length_matrix <- length_dist %>% 
    pivot_wider(id_cols = lengths, names_from = repertoire_id, 
        values_from = percent) %>%
    arrange(lengths) %>% 
    select(-lengths) %>%
    as.matrix()

length_rownames <- length_dist %>% 
    arrange(lengths) %>%
    pull(lengths) %>% 
    unique()

rownames(length_matrix) <- length_rownames

column_breaks <- annotate_table %>% 
    select(repertoire_id, time_point) %>% 
    distinct() %>%
    select(time_point) %>% 
    pull(time_point)

library(circlize)
length_heatmap <- draw(Heatmap(length_matrix, 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE, 
                          column_split = column_breaks, 
                          col = colorRamp2(c(0, 20, 40), c("#f0f9e8", "#7bccc4", "#0868ac")),
                          show_column_names = FALSE, 
                          name = "Percentage\nsequences",
                          heatmap_legend_param = list(
                              legend_direction = "horizontal", 
                              legend_width = unit(6, "cm"))),
                      heatmap_legend_side = "bottom")

pdf(file = "Figure6_CFAR_length.pdf", width = 14, height = 5)
length_heatmap
dev.off()