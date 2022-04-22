#' Load datasets
#' @description 
#' Load the selected datasets into tables using LymphoSeq. Multiple views are 
#' generated from the data and stored as RDA objects that are passed to the
#' study data merging module.
#' 
#' @param data_path Path to ImmunoSeq data
#' @param min_sequences The minimum number of sequences needed in a sample
library(LymphoSeq2)
library(tidyverse)
set.seed(12357) # Set the seed for reproducibility
# Accept user input
option_list <- list(
    optparse::make_option(c("-m", "--min_sequences"), 
                type="integer",  
                help="Path to ImmunoSEQ file",
                dest="min_seq")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE) # Parse the arguments
# Get list of all data files in the study
file_list <- base::list.files(path = "./", pattern = "tsv", full.names = TRUE)
# Function that checks if the number of sequences in each file is greater than
# equal to the minimum require number of sequences
getValidFiles <- function(airr_path) {
    read_immuno_seq <- function(airr_path, ...) {
        tryCatch(LymphoSeq2::readImmunoSeq(airr_path), error = function(c) {
            c$message <- paste0(c$message, " (in ", airr_path, ")")
            stop(c)
        })
    }
    study_table <- read_immuno_seq(airr_path)
    repertoire_id <- study_table %>% 
        dplyr::pull(repertoire_id) %>% 
        base::unique()
    if (base::sum(study_table$duplicate_count) >= parser$min_seq) {
        filtered_data <- tibble::tibble(sample = repertoire_id,
            valid = TRUE, 
            number_of_sequences = base::sum(study_table$duplicate_count),
            file_name = airr_path)
    } else {
        filtered_data <- tibble::tibble(sample = repertoire_id,
            valid = FALSE, 
            number_of_sequences = base::sum(study_table$duplicate_count),
            file_name = airr_path)
    }
    return(filtered_data)
}

validity_table <- file_list %>% 
    purrr::map(getValidFiles) %>% 
    dplyr::bind_rows()

valid_file_list <- validity_table %>% 
    dplyr::filter(valid == TRUE) %>% 
    dplyr::pull(file_name) 
# Copy all the valid files to a separate folder for nextflow to emit
base::file.copy(valid_file_list, "filtered")
# Write a table describing the files that were accepted and those that weren't
readr::write_csv(validity_table, file = "Sequence_count_summary.csv")

