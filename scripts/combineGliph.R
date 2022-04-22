library(tidyverse)     
library(lubridate)
library(readxl)
library(optparse)
# Accept arguments
option_list <- list(
    optparse::make_option(c("-i", "--input_path"), 
                type="character",  
                help="Path to length grouped gliph output",
                dest="ipath"),
    optparse::make_option(c("-p", "--patient_id"), 
                type="character",  
                help="Patient ID",
                dest="ptid")                
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                                print_help_and_exit = TRUE)
# Load all gliph tables and combine outputs
gfiles <- list.files(pattern = "*.tsv")
gtable <- gfiles %>% 
            purrr::map(readr::read_tsv) %>% 
            dplyr::bind_rows()
out_path <- paste(parser$ptid, "tsv", sep = ".")
write_tsv(gtable, out_path)