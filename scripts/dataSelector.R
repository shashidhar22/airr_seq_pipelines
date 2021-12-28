#' Select sample
#' @description:
#' This script randomly selects a specified number of normal sample to include 
#' in the data analysis. It also removes samples with total number of sequences 
#' less than a specified threshold. To do this the script takes in two user 
#' inputs.
#' * The number of normal samples to select
#' * The minimum number of sequences per sample
#' @param norm_count The number of normal samples to select
#' @param min_seq The minimum number of sequences per sample
#' @import tidyverse 
library(tidyverse)
library(optparse)
library(readxl)
library(LymphoSeq2)
set.seed(12357) # Set the seed for reproducibility
option_list <- list(
    optparse::make_option(c("-r", "--study_table"), 
        type="character",  
        help="Path to RDA file with study data",
        dest="study_table"),
    optparse::make_option(c("-n", "--nucleotide_table"), 
        type="character",  
        help="Path to RDA file with nucleotide level aggregation data",
        dest="nucleotide_table"),   
    optparse::make_option(c("-a", "--amino_table"), 
        type="character",  
        help="Path to RDA file with amino acid level aggregation data",
        dest="amino_table"), 
    optparse::make_option(c("-s", "--summary_table"), 
        type="character",  
        help="Path to RDA file with study summary table",
        dest="summary_table"),   
    optparse::make_option(c("-m", "--meta_table"),  
        type = "character",
        help = "Path to the meta data RDA object",
        dest = "meta_table"),
    optparse::make_option(c("-c", "--sample_count"), 
        type="integer",  
        help="Number of samples to be included in the study",
        dest="study_count"),
    optparse::make_option(c("-t", "--min_seq"),
        type="integer",
        help="Minimum number of sequences per sample",
        dest="min_seq"),
    optparse::make_option(c("-i", "--study_id"),        
        type = "character",
        help = "Study ID",
        dest = "study_id")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE) # Parse the arguments

load(parser$study_table) # Load the study data
load(parser$nucleotide_table) # Load the nucleotide level aggregation data
load(parser$amino_table) # Load the amino acid level aggregation data
load(parser$summary_table) # Load the study summary table  
#' Filter out samples that do not meet the minimum sequence threshold
summary_table <- summary_table %>% 
    filter(total_sequences >= parser$min_seq)
#' Sample datasets
#' 
#' Randomly select the number of samples from the normal samples
#' 
#' @param sample_count The number of samples to select
#' @return A tibble containing a subset of sample names based on group of data
#' and the number of samples specified by the user
sample_datasets <- function(sample_table, sample_count) {
    set.seed(12357) # Set the seed for reproducibility
    sample_group <- sample_table %>% # Get the group name of the chunk
        dplyr::pull(group) %>%
        unique()
    # If group == "Dean", then select sample_count number of samples
    if (sample_group == "Dean") {
        sample_table <- sample_table %>%
            dplyr::sample_n(size = sample_count)
    } 
    return(sample_table)
}
# Sample dataset based on group of the dataset. All CFAR samples are to be
# selected. Dean samples are to be selected based on the value of min_seq
# parameter.
summary_table <- summary_table %>%
    dplyr::mutate(group = dplyr::if_else(
        stringr::str_detect(repertoire_id, "CFAR"), "CFAR", "Dean")) %>%
    dplyr::group_by(group) %>%
    dplyr::group_split() %>%
    purrr::map(~sample_datasets(.x, parser$study_count)) %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(file_name = str_c(repertoire_id, "tsv", sep = "."))

# Write the sample data summary to a file for future reference
output_summary <- paste(parser$study_id, "data_summary.tsv", sep = "_")
write_tsv(summary_table, output_summary)

# Subsample the AIRR data tables to include only the selected samples
study_table <- study_table %>%
    filter(repertoire_id %in% summary_table$repertoire_id)
nucleotide_table <- nucleotide_table %>%
    filter(repertoire_id %in% summary_table$repertoire_id)
amino_table <- amino_table %>%
    filter(repertoire_id %in% summary_table$repertoire_id)

# Write the sample data to a file for future reference
save(study_table, file=paste(parser$study_id, "study.rda", sep = "_"))
save(nucleotide_table, file=paste(parser$study_id, "nucleotide.rda", sep = "_"))
save(amino_table, file=paste(parser$study_id, "amino.rda", sep = "_"))




