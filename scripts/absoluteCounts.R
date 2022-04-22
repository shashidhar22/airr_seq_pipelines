library(tidyverse)
library(ggpubr)
library(lubridate)
library(patchwork)
library(coin)
library(LymphoSeq2)

option_list <- list(
    optparse::make_option(c("-m", "--meta_data"), 
                type="character",  
                help="Path to meta data RDA file",
                dest="mdata"),
    optparse::make_option(c("-a", "--amino_data"), 
                type="character",  
                help="Path to amino data RDA file",
                dest="adata"),
    optparse::make_option(c("-s", "--study_data"), 
                type="character",  
                help="Path to study data RDA file",
                dest="sdata"),
    optparse::make_option(c("-i", "--immuno_data"), 
                type="character",  
                help="Path to ImmunoSeq summary data RDA file",
                dest="idata")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                               print_help_and_exit = TRUE)
# Load data
load(parser$mdata)