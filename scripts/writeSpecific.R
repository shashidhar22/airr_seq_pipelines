library(LymphoSeq2)
library(tidyverse)     
library(lubridate)
library(readxl)
library(patchwork)
library(lubridate)
library(readxl)
library(viridis)
library(fossil)
library(tidygraph)
library(ggraph)
library(optparse)
# Accept arguments
option_list <- list(
    make_option(c("-c", "--cfar_hiv"), 
                type="character",  
                help="Path to RDA file contain CFAR study specific sequence data",
                dest="chiv")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
# Load HIV specific amino acid sequences
load(parser$chiv)
# Split by repertoire and write to file
split_write <- function(htable) {
    repertoire <- htable %>%
                dplyr::pull(repertoire_id) %>%
                unique()
    out_file <- paste(repertoire, "specific.tsv", sep="_")
    write_tsv(htable, out_file)
}
htable %>% 
dplyr::group_by(repertoire_id) %>%
dplyr::group_split() %>%
purrr::map(split_write)