#' Load datasets
#' @description 
#' Load the selected datasets into tables using LymphoSeq. Multiple views are 
#' generated from the data and stored as RDA objects that are passed to the
#' study data merging module.
#' 
#' @param data_path Path to ImmunoSeq data
#' @param sample_name Tuple containing sample name and file name
#' @param sampling_rate The depth to which the data is sampled
#' @param n_iterations Number of iterations to perform
library(LymphoSeq2)
library(tidyverse)
set.seed(12357) # Set the seed for reproducibility
# Accept user input
option_list <- list(
    optparse::make_option(c("-d", "--data_path"), 
                type="character",  
                help="Path to ImmunoSEQ file",
                dest="data_path")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE) # Parse the arguments

# Read a ImmunoSEQ file and generate multiple views
study_table <- LymphoSeq2::readImmunoSeq(parser$data_path) %>% 
    filter(nchar(junction_aa) >= 6 & 
        str_sub(junction_aa, start = 1, end = 1 ) == "C" & 
        str_sub(junction_aa, start = -1L, end = -1L) == "F")
nucleotide_table <- LymphoSeq2::productiveSeq(study_table, 
    aggregate = "junction")
amino_table <- LymphoSeq2::productiveSeq(study_table,
    aggregate = "junction_aa")
summary_table <- LymphoSeq2::clonality(study_table)
mean_tcount <- nucleotide_table %>% 
    group_by(repertoire_id, junction_aa) %>% 
    summarize(tcount = n()) %>% ungroup() %>% 
    group_by(repertoire_id) %>% 
    summarize(mean_tcount = mean(tcount))
summary_table <- left_join(summary_table, 
    mean_tcount, 
    by = "repertoire_id")
# Define output file names
sample_name <- summary_table %>% 
    pull(repertoire_id) %>% 
    unique()
study_file <- paste(sample_name, "_study.rda", sep = "")
nucleotide_file <- paste(sample_name, "_nucleotide.rda", sep = "")
amino_file <- paste(sample_name, "_amino.rda", sep = "")
summary_file <- paste(sample_name, "_summary.rda", sep = "")
# Write the study data to a file
save(study_table, file = study_file)
save(nucleotide_table, file = nucleotide_file)
save(amino_table, file = amino_file)
save(summary_table, file = summary_file)