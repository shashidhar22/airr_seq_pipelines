library(LymphoSeq2)
library(tidyverse)     
library(lubridate)
library(readxl)
library(patchwork)
library(lubridate)
library(readxl)
library(viridis)
library(fossil)
library(optparse)
# Accept arguments
option_list <- list(
    make_option(c("-a", "--amino_data"), 
                type="character",  
                help="Path to productive sequence RDA file",
                dest="adata"),
    make_option(c("-m", "--meta_data"), 
                type="character",  
                help="Path to metadta RDA file",
                dest="meta")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
# Load metadata and productive amino acid table
load(parser$meta)
load(parser$adata)
# Merge with metadata
atable <- dplyr::left_join(atable, sample_table, 
                            by = "repertoire_id") %>%
            dplyr::left_join(patient_table, by = "patientID")
# Filter out samples that show consistent response
cresponse <- c(1673, 1838, 2027, 2312, 1334, 1402, 1423,
                1508, 1632, 1669, 1707, 1711, 1887, 2003,
                2034)
atable <- atable %>% 
            dplyr::filter(patientID %in% cresponse)
htable <- atable %>% 
            dplyr::filter(prevalence == 0)
lptable <- atable %>%
            dplyr::filter(prevalence < 100 & prevalence > 0)
ctable <- atable %>% dplyr::filter(prevalence == 100)
# Save all high, low and 
save(htable, file = "CFAR_specific.rda")
save(lptable, file = "CFAR_prevalent.rda")
save(ctable, file = "CFAR_common.rda")         