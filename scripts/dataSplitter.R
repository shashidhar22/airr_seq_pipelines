#' Split data table
#' @description 
#' In order to completely leverage the power of nextflow, we want to have a 
#' function that can split the data table into multiple smaller tables for all
#' views generated by the pipeline.
#' @param study_table The raw AIRR seq data from the study
#' @param amino_table The amino acid data from the study
#' @param nuc_table The nucleotide data from the study
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
        dest="amino_table")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE) # Parse the arguments

load(parser$study_table) # Load the study data
load(parser$nucleotide_table) # Load the nucleotide data
load(parser$amino_table) # Load the amino acid data

#' Split and write
#' @description
#' Split the data table into multiple smaller tables for all views generated by
#' the pipeline.
#' @param data_table AIRR Seq data table
#' @param data_type The type of data table to split
split_data <- function(data_table, data_type) {
    # Get sample name
    sample_name <- data_table %>% 
        pull(repertoire_id) %>% 
        unique() 
    if (data_type == "study") {
        # Write sample study table
        out_file <- paste(sample_name, "_study.tsv", sep="")
        write_tsv(data_table, out_file)
    } else if (data_type == "nucleotide") {
        # Write sample nucleotide table
        out_file <- paste(sample_name, "_nucleotide.tsv", sep="")
        write_tsv(data_table, out_file)
    } else if (data_type == "amino") {
        # Write sample amino acid table
        out_file <- paste(sample_name, "_amino.tsv", sep="")
        write_tsv(data_table, out_file)
    }
}

# Split and write per sample study table
study_table %>% 
group_by(repertoire_id) %>%
group_split() %>% 
map(split_data, "study") 

# Split and write per sample nucleotide table
nucleotide_table %>%
group_by(repertoire_id) %>%
group_split() %>%
map(split_data, "nucleotide")

# Split and write per sample amino acid table
amino_table %>%
group_by(repertoire_id) %>%
group_split() %>%
map(split_data, "amino")


