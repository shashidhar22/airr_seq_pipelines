# Load all libraries
library(tidyverse)     
library(lubridate)
library(readxl)
library(optparse)
library(devtools)
devtools::load_all("~/projects/LymphoSeq2")
# Accept arguments
option_list <- list(
    optparse::make_option(c("-t", "--tcr_data"), 
                        type = "character",  
                        help = "Path to TCR data folder",
                        dest = "tdata")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                            print_help_and_exit = TRUE)
# Load data using LymphoSeq, tranform it to get productive sequences at amino
# acid level, nucleotide level and repertoire summary table.
stable <- LymphoSeq2::readImmunoSeq(parser$tdata) 
atable <- LymphoSeq2::productiveSeq(stable) 
ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")
itable <- LymphoSeq2::clonality(stable)
# Pull repertoire name and use for naming the files
repertoire_id <- stable %>% 
                dplyr::pull(repertoire_id) %>% 
                unique()
# Save all repertoire related tables into RDA format
save(stable, file = paste(repertoire_id, "stable.rda", sep = "_"))
save(atable, file = paste(repertoire_id, "atable.rda", sep = "_"))
save(ntable, file = paste(repertoire_id, "ntable.rda", sep = "_"))
save(itable, file = paste(repertoire_id, "itable.rda", sep = "_"))
