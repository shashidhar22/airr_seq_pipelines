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
        str_sub(junction_aa, start = -1L, end = -1L) == "F") %>%
    mutate(v_call = if_else(str_detect(v_call,"\\*"), 
            str_c(
                str_c("TRBV", as.numeric(str_extract(v_call, "\\d+"), sep = "")), 
                    as.numeric(str_extract(v_call,"\\d+$")), sep = "-"),
            str_c("TRBV", as.numeric(str_extract(v_call, "\\d+")), sep = "")),
        j_call = if_else(str_detect(j_call,"\\*"),
            str_c(
                str_c("TRBJ", as.numeric(str_extract(j_call, "\\d+"), sep = "")), 
                    as.numeric(str_extract(j_call,"\\d+$")), sep = "-"),
            str_c("TRBJ", as.numeric(str_extract(j_call, "\\d+")), sep = "")),
        d_call = if_else(str_detect(d_call,"\\*"),
            str_c(
                str_c("TRBD", as.numeric(str_extract(d_call, "\\d+"), sep = "")), 
                    as.numeric(str_extract(d_call,"\\d+$")), sep = "-"),
            str_c("TRBD", as.numeric(str_extract(d_call, "\\d+")), sep = "")))

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
