#' Prep control databases
#' @description 
#' Read the Emerson dataset and parse out all the relevant data and metadata
#' to generate the reference database
#' 
#' @param data_path Path to ImmunoSeq data
library(LymphoSeq2)
library(tidyverse)
set.seed(12357) # Set the seed for reproducibility
# Accept user input
option_list <- list(
    optparse::make_option(c("-d", "--data_path"), 
                type="integer",  
                help="Path to ImmunoSEQ file",
                dest="data_path")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
    print_help_and_exit = TRUE) # Parse the arguments

read_controls <- function(control_path) {
    control_data <- read_tsv(control_path, show_col_types = FALSE)
    control_data <- control_data %>% 
        dplyr::mutate(age = str_extract(sample_catalog_tags, "Age:\\d+ Years"),
            sex = str_extract(sample_catalog_tags, "Biological Sex:\\w+,"),
            disease_status = str_extract(sample_catalog_tags, "Disease state:\\w+,"),
            ethnic_group = str_extract(sample_catalog_tags, "Ethnic Group:[:alnum:]+[:blank:]?[:alnum:]+"),
            hla = str_extract_all(sample_catalog_tags, "HLA MHC class I:HLA-[A,B]\\*\\d+,"),
            virus = str_extract(sample_catalog_tags, "Virus Diseases:[:alnum:]+[:blank:]?[:graph:]+")) %>%
        dplyr::select(sample_name, age, sex, disease_status, ethnic_group, hla, 
            virus, bio_identity, templates) %>% 
        dplyr::mutate(bio_identity = stringr::str_replace_all(bio_identity, "\\+", ",")) %>%
        tidyr::separate(bio_identity, into = c("junction_aa", "v_gene", "j_gene"),
            sep = ",") %>% 
        dplyr::mutate(v_call =  str_c(str_c("TRBV", as.numeric(str_extract(v_gene, "\\d+"), sep = "")), as.numeric(str_extract(v_gene,"\\d+$")), sep = "-"),
                j_call = str_c(str_c("TRBJ", as.numeric(str_extract(j_gene, "\\d+"), sep = "")), as.numeric(str_extract(j_gene,"\\d+$")), sep = "-"),
                age = as.integer(str_extract(age, "\\d+")),
                sex = stringr::str_remove(sex, "Biological Sex:"),
                sex = stringr::str_remove(sex, ","),
                disease_status = stringr::str_remove(disease_status, "Disease state:"),
                disease_status = stringr::str_remove(disease_status, ","),
                ethnic_group = stringr::str_remove(ethnic_group, "Ethnic Group:"),
                ethnic_group = stringr::str_remove(ethnic_group, ","),
                virus = stringr::str_remove(virus, "Virus Diseases:"),
                virus = stringr::str_remove(virus, ",")) %>% 
        dplyr::filter(junction_aa != "X")
    return(control_data)                
}

control_table <- read_controls(parser$data_path)
out_file <- control_table %>% 
    dplyr::pull(sample_name) %>% 
    base::unique() %>% 
    stringr::str_c("rda", sep = ".")
save(control_table, file = out_file)

