#' Find public sequences
#' @description
#' There are multiple databases that curate the antigenic specificity of 
#' TCR sequences. This script will search the McPAS and VDJdb databases for
#' the prevalence of public sequences in the repertoire. The function will take
#' a tibble of productive sequences and return a tibble with two additional
#' columns for the antigenic specificity of the sequences. The first column
#' will be for the McPAS database and the second for the VDJdb database. 
library(LymphoSeq2)
library(tidyverse)     
library(readxl)
library(optparse)
# Accept arguments
option_list <- list(
    make_option(c("-d", "--database"),
        type = "character",
        help = "Comma separated list of paths to the database containing the public sequences",
        dest="database"),
    make_option(c("-m", "--mode"),
        type = "character",
        help = "Mode of lookup. Strict: Match VDJ gene. Lenient: Match just junction seuqences",
        dest = "mode")
)
parser <- parse_args(OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE)

#' Load databases
#' @description 
#' Load public databases into a table containing the TRB nucleotide and amino
#' acid sequences, the VDJ gene calls and the associated pathology. The function
#' will filter out only Human public sequences from VDJDb
load_database <- function(db_path) {
    db_table <- read_tsv(db_path) %>% 
        filter(gene == "TRB" & species == "HomoSapiens") %>% 
        select(cdr3, v.segm, j.segm, antigen.species) %>% 
        rename(junction_aa = cdr3, 
            v_call = v.segm, j_call = j.segm,
            pathology = antigen.species) %>% 
        mutate(pathology = recode(pathology, `M. tuberculosis` = "M.Tuberculosis",
            "Psoriatic arthritis" = "Psoriatic Arthritis",
            "Hepatitis C virus (HCV)" = "Hepatitis C virus"))
    return(db_table)   
}
# Load public databases and amino acid tables
db_table <- load_database(parser$database)
amino_file <- list.files(pattern = "_amino.rda")
load(amino_file)
# Merge tables
if ( parser$mode == "strict" ) {
    amino_table <-  left_join(amino_table, db_table, 
        by = c("junction_aa", "v_call", "j_call"))
} else if ( parser$mode == "lenient") {
    # In lenient mode we just want to merge the table based on the junction
    # sequences and the pathology, that's why we remove v_call and j_call
    db_pathology <- db_table %>% 
        select(junction_aa, pathology) 
    amino_table <- left_join(amino_table, db_pathology, by = "junction_aa")
}

# Get a summary table of all the public sequences found in the sample and the
# cumulative frequency
public_summary <- amino_table %>% 
    group_by(pathology) %>% 
    summarize(public_count = n(),
        public_frequency = sum(duplicate_frequency))
# Output summary table and amino table with the pathology listed
sample_name <- str_remove(amino_file, "_amino.rda")
summary_file <- paste(sample_name, "public_summary.tsv", sep = "_")
annotated_file <- paste(sample_name, "public_annotated.rda", sep = "_")
write_tsv(public_summary, summary_file)
save(amino_table, file = annotated_file)