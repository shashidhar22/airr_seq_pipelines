library(tidyverse)     
library(lubridate)
library(readxl)
library(optparse)
# Accept arguments
option_list <- list(
    optparse::make_option(c("-a", "--amino_data"), 
        type="character",  
        help="Path to amino acid RDA file",
        dest="adata"),
    optparse::make_option(c("-p", "--public_path"), 
            type="character",  
            help="Path to public tcr databases",
            dest="ppath")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE)
load(parser$adata)
# Function definitions
readPublic <- function(file_path) {
    ptable <- read_csv(file_path)
    sample_name <- basename(file_path) %>% str_extract("^\\w+")
    if ("aminoAcid.beta" %in% colnames(ptable)) {
        ptable <- ptable %>%
            pivot_longer(cols= c(aminoAcid.beta),
                names_to = "Chain",
                values_to = "aminoAcid") %>%
            filter(Chain == "aminoAcid.beta") %>%
            mutate(antigen = sample_name) %>%
            select(antigen, aminoAcid)
    } else if ("CDR3" %in% colnames(ptable)){
        ptable <- ptable %>%
            mutate(antigen = sample_name) %>%
            rename(aminoAcid = CDR3) %>%
            select(antigen, aminoAcid)
    } else {
        ptable <- ptable %>%
            mutate(antigen = sample_name) %>%
            select(antigen, aminoAcid)              
    }
    ptable <- ptable %>% drop_na()
    return(ptable)
}

# Load public TRB data
public_list <- list.files(parser$ppath, full.names = TRUE, all.files = FALSE, 
    recursive = FALSE, pattern = ".csv", 
    include.dirs = FALSE)
pbtable <- public_list %>% 
    purrr::map(readPublic) %>% 
    dplyr::bind_rows() %>% 
    dplyr::rename(junction_aa = aminoAcid)
mcppath <- paste(parser$ppath, "publicDB/McPAS-TCR.csv", sep = "/") 
mptable <- readr::read_csv(mcppath) %>% 
    dplyr::rename(junction_aa = CDR3.beta.aa, antigen = Pathology) %>% 
    dplyr::select(antigen, junction_aa)

antigen_table <- c(EBV = "Epstein Barr virus (EBV)", 
    CMV = "Cytomegalovirus (CMV)",
    flu = "Influenza",
    HSV = "Herpes simplex virus 1 (HSV1)",
    `Human herpes virus 1` = "Herpes simplex virus 1 (HSV1)",
    `Hepatitis C virus` = "Hepatitis C virus (HCV)",
    HIV = "Human immunodeficiency virus (HIV)",
    KS = "Kaposi Sarcoma",
    `Psoriatic arthritis` = "Psoriatic Arthritis",
    `HIV-1` = "Human immunodeficiency virus (HIV)",
    InfluenzaA = "Influenza",
    HPV = "Human papilloma virus",
    HCV = "Hepatitis C virus (HCV)",
    `HSV-2` = "Herpes simplex virus 2 (HSV2)",
    LCMV = "Lymphocytic choriomeningitis virus (LCMV)",
    YFV = "Yellow fever virus")

filter_species <- c("MusMusculus", "HomoSapiens", "GallusGallus")
vdjpath <- paste(parser$ppath, "publicDB/vdjdb_full.txt", sep = "/") 
vdjtable <- readr::read_delim(vdjpath, , delim="\t" ) %>%
    dplyr::select(cdr3.beta, antigen.species) %>%
    dplyr::rename(junction_aa = cdr3.beta, antigen = antigen.species) %>%
    dplyr::filter(!antigen %in% filter_species)

pbtable <- dplyr::bind_rows(pbtable, mptable) %>%
    dplyr::bind_rows(vdjtable) %>%
    dplyr::mutate(antigen = recode_factor(antigen, !!!antigen_table)) %>%
    dplyr::distinct_all() 

patable <- left_join(atable, pbtable, by = "junction_aa") 

save(patable, file = "patable_cfar.rda")
