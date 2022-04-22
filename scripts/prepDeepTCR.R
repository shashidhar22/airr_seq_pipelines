library(tidyverse)     
library(lubridate)
library(readxl)
library(optparse)
# Accept arguments
option_list <- list(
    optparse::make_option(c("-a", "--amino_data"), 
                type="character",  
                help="Path to annotated amino acid collapsed RDA file",
                dest="apath"),
    optparse::make_option(c("-m", "--hla_data"), 
            type="character",  
            help="Path to HLA table, the first column should correspond to sample name",
            dest="hpath"),
    optparse::make_option(c("-s", "--study_id"), 
                type="character",  
                help="Study ID",
                dest="sid")            
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE)
gliph_amino_table <- readRDS(parser$apath)

atable <- gliph_amino_table %>% 
    mutate(cdrA = NA) %>%
    select(junction_aa, v_call, j_call, cdrA, gliph_id, duplicate_frequency) 

htable <- read_tsv(parser$hpath) # HLA table
sid <- parser$sid

hla_table <- atable %>% 
    select(repertoire_id, patient_id) %>% 
    distinct() %>% 
    left_join(htable, by = "patient_id") %>% 
    select(-patient_id) %>% 
    mutate(repertoire_id = str_c(repertoire_id, "tsv", sep = ".")) %>% 

out_file <- paste(parser$sid, "HLA.csv", sep = "_")
write_csv(hla_table, file = out_file, col_names = FALSE)