#' Merge sample tables
#' @description:
#' This script collects and merges the following view from LymphoSeq
#' * Raw study table
#' * Nucleotide level aggregation of the study table
#' * Amino acid level aggregation of the study table
#' * Study summary table
#' @import tidyverse 
library(tidyverse)
library(optparse)
library(readxl)
library(LymphoSeq2)
set.seed(12357) # Set the seed for reproducibility
# Accept arguments
option_list <- list(
    make_option(c("-s", "--study_id"), 
                type="character",  
                help="Study ID",
                dest="study_id")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
# Load and merge study tables. Here a function is need because RDA stores the 
# variable names as well, hence it cannot be merged directly.
load_study <- function(study_file) {
    load(study_file)
    return(study_table)
}
study_table <- list.files(pattern = "_study.rda") %>% 
    purrr::map(load_study) %>%
    dplyr::bind_rows()
save(study_table, file = "study_table.rda")
# Load and merge nucleotide level aggregation tables.
load_nucleotide <- function(study_file) {
    load(study_file)
    return(nucleotide_table)
}
nucleotide_table <- list.files(pattern = "_nucleotide.rda") %>% 
    purrr::map(load_nucleotide) %>%
    dplyr::bind_rows()
save(nucleotide_table, file = "nucleotide_table.rda")
# Load and merge amino acid level aggregation tables.
load_amino_acid <- function(study_file) {
    load(study_file)
    return(amino_table)
}
amino_table <- list.files(pattern = "_amino.rda") %>% 
    purrr::map(load_amino_acid) %>%
    dplyr::bind_rows()
save(amino_table, file = "amino_table.rda")
# Load and merge study summary tables.
load_summary <- function(study_file) {
    load(study_file)
    return(summary_table)
}
summary_table <- list.files(pattern = "_summary.rda") %>%
    purrr::map(load_summary) %>%
    dplyr::bind_rows()
save(summary_table, file = "summary_table.rda")

# Save object for LymphoSeq Shiny input
out_file <- paste(parser$study_id, "LS.rda", sep = "_")
save(study_table, nucleotide_table, amino_table, summary_table,
    file = out_file)