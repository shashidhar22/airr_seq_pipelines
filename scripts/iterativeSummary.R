#' Iterative calculation of summary statistics
#' @description
#' Given variation in sequencing depth, we need to sample the data in order to 
#' estimate the summary statistics. This script does this iteratively.
#' @param study_table The sampled study table
#' @param iterations The number of iterations
#' @param minimum_sequences The minimum number of samples per sample
#' @param meta_table The metadata table
#' @return The averaged summary statistics
library(tidyverse)
library(LymphoSeq2)
library(readxl)
library(optparse)
set.seed(12357) # Set the seed for reproducibility
option_list <- list(
    optparse::make_option(c("-r", "--study_table"), 
        type="character",  
        help="Path to TSV file with study data",
        dest="study_table"),
    optparse::make_option(c("-i", "--iterations"),
        type="numeric",
        help="Number of iterations",
        dest="iterations"),
    optparse::make_option(c("-c", "--minimum_sequences"),
        type="numeric",
        help="Minimum number of sequences",
        dest="min_sequences"),
    optparse::make_option(c("-m", "--meta_table"),
        type="character",
        help="Path to RDA file with metadata",
        dest="meta_table")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE) # Parse the arguments

# Load the study table
study_table <- read_tsv(parser$study_table)
nucleotide_table <- productiveSeq(study_table, aggregate = "junction")

#' Random sample data
#' @description
#' The data from each sample will be sampled down to the minimum number of 
#' sequences sequences and diversity will be calculated for each sample. This 
#' process will be repeated 100 times and the results will be averaged.
#' @param study_table The raw study table table
#' @param minimum_sequences The minimum number of sequences per sample
#' @return Tibble with the diversity and number of sequences per sample
randomSample <- function(study_table,  min_sequences) {
    top_aminos <- study_table %>% 
        productiveSeq() %>% 
        topSeqs(top = 100) %>% 
        pull(junction_aa) %>% 
        unique()
    study_table <- study_table %>% 
        uncount(duplicate_count) %>% 
        sample_n(size = min_sequences, replace = FALSE) %>% 
        dplyr::group_by(repertoire_id, junction, junction_aa, v_call, j_call, 
            d_call, v_family, j_family, d_family, reading_frame) %>% 
        dplyr::summarize(duplicate_count = n()) %>% 
        dplyr::ungroup() %>%
        dplyr::mutate(duplicate_frequency = duplicate_count/sum(duplicate_count)) %>%
        dplyr::select(repertoire_id, everything())
    repertoire <- study_table %>% 
        pull(repertoire_id) %>%
        unique()
    summary_table <- clonality(study_table)
    translation_efficacy <- study_table %>% 
        productiveSeq(aggregate = "junction") %>%
        filter(junction_aa %in% top_aminos) %>%
        group_by(junction_aa) %>% 
        summarize(translation_efficacy = n()) %>% 
        pull(translation_efficacy) %>% 
        mean()
    summary_table <- summary_table %>%
        mutate(translation_efficacy = translation_efficacy)
    print(summary_table)
    return(summary_table)
}
#' Iterative sampling
#' @description
#' Iteratively samples each repertoire as specified by the user and 
# calculates the diversity of each sample
bootstrapSummary <- function(study_table, iterations, min_sequences) {
    summary_table <- rerun(randomSample(study_table, min_sequences), .n = iterations) %>% 
        bind_rows() 
    
    #    group_by(repertoire_id) %>% 
    #    summarize(clonality = mean(clonality), 
    #        gini = mean(giniCoefficient), 
    #        unique_sequences = mean(unique_sequences),
    #        total_sequences = mean(total_sequences),
    #        translation_efficacy = mean(translation_efficacy)) 
    return(summary_table)
}

# Get averaged summary statistics table
summary_table <- bootstrapSummary(study_table, parser$iterations, 
    parser$min_sequences)

# Get translation efficacy tables
top_aminos <- study_table %>% 
        productiveSeq() %>% 
        topSeqs(top = 100) %>% 
        pull(junction_aa) %>% 
        unique()
translation_efficacy <- study_table %>%
    productiveSeq(aggregate = "junction") %>%
    filter(junction_aa %in% top_aminos) %>%
    group_by(repertoire_id, junction_aa) %>% 
    summarize(translation_efficacy = n()) 

# Save the summary table
repertoire_id = summary_table %>% 
    pull(repertoire_id) %>%
    unique()
out_file = paste(repertoire_id, "_avg_summary.rda", sep = "")
save(summary_table, file = out_file)

# Save the translation efficacy table
out_file = paste(repertoire_id, "_translation_efficacy.rda", sep = "")
save(translation_efficacy, file = out_file)
