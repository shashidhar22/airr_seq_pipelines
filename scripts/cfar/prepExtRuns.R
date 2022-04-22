#' Prepare HLA information in the required format for external tools
#' @description
#' Tools like GLIPH2, DeepTCR and tcrdist require the HLA information to be 
#' consistent and in a particular format. In our current study, the CFAR cohort 
#' and the Dean cohort have different formatting of the HLA information. These 
#' need to be made consistent and output in the required format.
#' @param meta_data The RDA file containing all the metadata information 
#' @param study_id The study ID of the cohort
#' @return Data frames with the HLA information in the required format
library(tidyverse)
library(LymphoSeq2)
library(readxl)
library(optparse)
set.seed(12357) # Set the seed for reproducibility
option_list <- list(
    optparse::make_option(c("-a", "--amino_table"),
        type = "character",
        help = "Path to the sampled amino acid table",
        dest = "amino_table"),
    optparse::make_option(c("-m", "--meta_table"),
        type="character",
        help="Path to RDA file with metadata",
        dest="meta_table"),
    optparse::make_option(c("-n", "--min_seuqences"),
        type="integer",
        help="Minimum number of sequences to be considered",
        dest="min_sequences"),
    optparse::make_option(c("-s", "--study_id"),
        type="character",
        help="Study ID of the cohort",
        dest="study_id")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE) # Parse the arguments

# Load the RDA files
load(parser$meta_table)
load(parser$amino_table)
# Format the amino acid table to IGMT format VDJ gene names
amino_table <- amino_table %>%
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
# Add the patient ID to the amino acid table and generate GLIPH ID
amino_table <- amino_table %>% 
    left_join(study_sample_table, by = "repertoire_id") %>%
    mutate(status = if_else(time_point == "Normal", "Normal", "HIV"))
# Now that the amino table has the patient information, we can generate a table
# that has the information of the repertoire id and the patient id. This can 
# be used to transform the HLA table to contain the repertoire_id as the first
# column.
sample_table <- amino_table %>% 
    select(repertoire_id, patientID) %>%
    distinct()
# Format metadata for GLIPH runs
gliph_hla_table <- study_hla_table %>% 
    right_join(sample_table, by = "patientID") %>%
    select(-patientID, -years_rel_to_art, -timePoint, -time_group, -sex, -age_at_collection) %>%
    pivot_longer(-repertoire_id, names_to = "name", values_to = "hla") %>%
    mutate(hla = str_extract(hla, "\\w+\\*\\d+")) %>% 
    pivot_wider(id_cols = repertoire_id, names_from = name, values_from = hla,
        values_fn = first)

gliph_out_file <- paste(parser$study_id, "_gliph_hla.tsv", sep = "")    
write_tsv(gliph_hla_table, gliph_out_file, col_names = FALSE)
# Sample the minimum number of sequence per repertoire to create the GLIPH 
# input file
# Gliph One: HIV vs Normal
gliph_amino_table_one <- amino_table %>% 
    mutate(gliph_id = str_c(repertoire_id, status, sep = ":")) %>%
    group_by(repertoire_id) %>%
    arrange(desc(duplicate_count)) %>% 
    slice_head(n = parser$min_sequences) %>%
    ungroup()
# Gliph Two: Pre-ART vs Post-ART vs Normal
gliph_amino_table_two <- amino_table %>% 
    mutate(gliph_id = str_c(repertoire_id, time_point, sep = ":")) %>%
    group_by(repertoire_id) %>%
    arrange(desc(duplicate_count)) %>% 
    slice_head(n = parser$min_sequences) %>%
    ungroup()

gliph_one_amino_file <- paste(parser$study_id, "_gliph_one_amino_table.rds", sep = "")
saveRDS(gliph_amino_table_one, file = gliph_one_amino_file)

gliph_two_amino_file <- paste(parser$study_id, "_gliph_two_amino_table.rds", sep = "")
saveRDS(gliph_amino_table_two, file = gliph_two_amino_file)

# Format metadata for DeepTCR runs
deeptcr_hla_table <- study_hla_table %>% 
    right_join(sample_table, by = "patientID") %>% 
    select(-patientID, -years_rel_to_art, -timePoint, -time_group, -sex, -age_at_collection) %>%
    pivot_longer(-repertoire_id, names_to = "name", values_to = "hla") %>%
    pivot_wider(id_cols = repertoire_id, names_from = name, values_from = hla,
        values_fn = first) %>%
    mutate(repertoire_id = str_c(repertoire_id, "tsv", sep = "."))
# Save deep TCR HLA files
deeptcr_out_file <- paste(parser$study_id, "_deep_tcr_hla.tsv", sep = "")
write_tsv(deeptcr_hla_table, deeptcr_out_file, col_names = FALSE)
# Write the output into different folder based on status and timePoint
# Normal list
sample_table <- sample_table %>% 
    left_join(study_sample_table, by = c("repertoire_id", "patientID")) %>%
    mutate(status = if_else(time_point == "Normal", "Normal", "HIV"))
normal_list <- sample_table %>% 
    filter(status == "Normal") %>%
    select(repertoire_id) %>%
    mutate(repertoire_id = str_c(repertoire_id, "tsv", sep = ".")) %>% 
    pull(repertoire_id)
hiv_list <- sample_table %>% 
    filter(status == "HIV") %>%
    select(repertoire_id) %>%
    mutate(repertoire_id = str_c(repertoire_id, "tsv", sep = ".")) %>% 
    pull(repertoire_id)
preart_list <- sample_table %>% 
    filter(time_point == "Pre-ART") %>%
    select(repertoire_id) %>%
    mutate(repertoire_id = str_c(repertoire_id, "tsv", sep = ".")) %>% 
    pull(repertoire_id)
postart_list <- sample_table %>%
    filter(time_point == "Post-ART") %>%
    select(repertoire_id) %>%
    mutate(repertoire_id = str_c(repertoire_id, "tsv", sep = ".")) %>% 
    pull(repertoire_id)
# Write the output into different folder based on status and timePoint
# DeepTCR 1
dir.create("deeptcr_one")
dir.create("deeptcr_one/normal")
dir.create("deeptcr_one/hiv")
file.copy(deeptcr_out_file, "deeptcr_one/")
file.copy(normal_list, "deeptcr_one/normal")
file.copy(hiv_list, "deeptcr_one/hiv")
#DeepTCR 2
dir.create("deeptcr_two")
dir.create("deeptcr_two/normal")
dir.create("deeptcr_two/preart")
dir.create("deeptcr_two/postart")
file.copy(deeptcr_out_file, "deeptcr_two/")
file.copy(normal_list, "deeptcr_two/normal")
file.copy(preart_list, "deeptcr_two/preart")
file.copy(postart_list, "deeptcr_two/postart")