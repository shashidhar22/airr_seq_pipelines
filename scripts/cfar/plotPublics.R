#' Generate the public sequences comparison figure. 
#' @description
#' Generate the box plots comparing the occurence of public sequences in control and PLHIV cohort
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

comp_palette <- c("Post-ART" = "#f1a340", "Pre-ART" = "#998ec3", "Normal" = "#000000")
comp_one <- list(c("Pre-ART", "Post-ART"), c("Pre-ART", "Normal"), c("Post-ART", "Normal"))
comp_two <- list(c("Pre-ART", "Post-ART"))

loadPublicSum <- function(public_path) {
    public_table <- read_tsv(public_path)
    repertoire_id <- public_path %>% 
        basename() %>% 
        str_remove("_public_summary.tsv")
    public_table <- public_table %>% 
        mutate(repertoire_id = repertoire_id)
    return(public_table)
}

public_table <- list.files(path = "./", 
        pattern = "_public_summary.tsv", full.names = TRUE) %>%
    purrr::map(loadPublicSum) %>%
    bind_rows()

annotate_table <- amino_table %>% 
    left_join(study_sample_table, by = "repertoire_id") %>% 
    left_join(study_patient_table, by = "patientID") %>%
    mutate(lengths = nchar(junction_aa)) %>% 
    select(repertoire_id, time_group, time_point, junction_aa, lengths) 

public_table <- annotate_table %>% 
    select(repertoire_id, time_group, time_point) %>% 
    distinct() %>% 
    inner_join(public_table, by = "repertoire_id")

ebv_plot <- public_table %>%
    filter(pathology == "EBV") %>% 
    ggboxplot(x = "time_point", y = "public_frequency", color = "time_point", palette = comp_palette, add = "jitter") + 
    stat_compare_means(comparisons = comp_one) +
    theme_classic(base_size = 20) +
    labs(color = "Time point", x = "Time point", y = "Frequency of annotated sequences") +
    yscale("log10", .format = TRUE) +
    ylim(0, 10e-3) +
    ggtitle("EBV")


cmv_plot <- public_table %>%
    filter(pathology == "CMV") %>% 
    ggboxplot(x = "time_point", y = "public_frequency", color = "time_point", palette = comp_palette, add = "jitter") + 
    stat_compare_means(comparisons = comp_one) +
    labs(color = "Time point", x = "Time point", y = "Frequency of annotated sequences") +
    yscale("log10", .format = TRUE) +
    ylim(0, 10e-3) +
    ggtitle("CMV") +
    theme_classic(base_size = 20) +
    theme(axis.title.y = element_blank())


hiv_plot <- public_table %>%
    filter(pathology == "HIV-1") %>% 
    ggboxplot(x = "time_point", y = "public_frequency", color = "time_point", palette = comp_palette, add = "jitter") + 
    stat_compare_means(comparisons = comp_one) +
    labs(color = "Time point", x = "Time point", y = "Frequency of annotated sequences") +
    yscale("log10", .format = TRUE) +
    ylim(0, 10e-3) +
    ggtitle("HIV") +
    theme_classic(base_size = 20) +
    theme(axis.title.y = element_blank())

public_plot <- ebv_plot + cmv_plot + hiv_plot + plot_layout(guides = "collect")

ggsave("SupFig1_CFAR_antigen_specific.pdf", public_plot, device = "pdf", width = 14,
    height = 5, units = "in", limitsize = FALSE)