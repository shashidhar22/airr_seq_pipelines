library(devtools)
library(tidyverse)     
library(lubridate)
library(readxl)
library(optparse)
devtools::load_all("~/projects/LymphoSeq2")
# Accept arguments
option_list <- list(
    make_option(c("-t", "--tcr_data"), 
                type="character",  
                help="Path to TCR data folder",
                dest="tdata"),
    make_option(c("-m", "--meta_data"), 
                type="character",  
                help="Path to metadta file",
                dest="meta")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
# The metadta is stored across multiple sheets
sheet_list <- c('CFAR_Sample_info', "CFAR_HLA_info", 
            "CFAR_Patient_info", "CFAR_Viral_load", 
            "CFAR_CD4_count", "CFAR_Coinfections")
# Load and process immunoSeq data
stable <- LymphoSeq2::readImmunoSeq(parser$tdata)
atable <- LymphoSeq2::productiveSeq(stable, prevalence = TRUE)
ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")
# Load in the metadata tables 
sample_table <- readxl::read_excel(parser$meta, sheet = sheet_list[1])
hla_table <- readxl::read_excel(parser$meta, sheet = sheet_list[2])
patient_table <- readxl::read_excel(parser$meta, sheet = sheet_list[3]) 
vl_table <- readxl::read_excel(parser$meta, sheet = sheet_list[4])
cc_table <- readxl::read_excel(parser$meta, sheet = sheet_list[5])
infection_table <- readxl::read_excel(parser$meta, sheet = sheet_list[6])
# Transform meta data

sample_table <- sample_table %>%
                dplyr::filter(!is.na(dateTCR)) %>%
                dplyr::mutate(patientID = as.integer(patientID),
                              patientID = forcats::as_factor(patientID),
                              dateTCR = ymd(dateTCR),
                              time_point = forcats::as_factor(timePoint),
                              repertoire_id = forcats::as_factor(sampleID)) %>%
                dplyr::select(-sampleID, -timePoint)
patient_table <- patient_table %>%
                 dplyr::filter(!is.na(dateART)) %>%
                 dplyr::mutate(patientID = as.integer(patientID),
                               patientID = forcats::as_factor(patientID),
                               dateART = lubridate::as_date(dateART),
                               sex = forcats::as_factor(sex),
                               dateMinEvidenceHIV = ymd(dateMinEvidenceHIV),
                               birthYear = lubridate::as_date(birthYear),
                               age_at_art = (birthYear %--% dateART) %/% years(1))
hla_table <- hla_table %>%
             dplyr::select(-dateART) %>%
             tidyr::pivot_longer(-patientID, 
                                 names_to = "allele",
                                 values_to = "value") %>%
             dplyr::mutate(patientID = as.integer(patientID),
                           patientID = forcats::as_factor(patientID),
                           HLA = str_remove(allele, '(allele\\d|Allele\\d)'),
                           HLA = str_to_upper(HLA),
                           HLA = paste("HLA", HLA, sep = " "), 
                           allele = str_extract(allele, '\\d$'),
                           value = str_remove(value, '\\w+\\*'),
                           value = str_replace(value, "Not_Present", "NA"))
vl_table <- vl_table %>%
            dplyr::mutate(patientID = as.integer(patientID),
                         patientID = forcats::as_factor(patientID),
                         dateCollection = lubridate::as_date(dateCollection)) %>%
            dplyr::left_join(patient_table, by = "patientID") %>%
            dplyr::mutate(time_point = dplyr::if_else(dateCollection <= dateART, "Pre-ART", "Post-ART"))
cc_table <- cc_table %>%
            dplyr::mutate(patientID = as.integer(patientID),
                        patientID = forcats::as_factor(patientID),
                        collectionDate = lubridate::as_date(collectionDate)) %>%
            dplyr::rename(dateCollection = collectionDate) %>%
            dplyr::left_join(patient_table, by = "patientID") %>%
            dplyr::mutate(time_point = dplyr::if_else(dateCollection <= dateART, "Pre-ART", "Post-ART"))

infection_table <- infection_table %>%
                    dplyr::mutate(patientID = as.integer(patientID),
                                patientID = forcats::as_factor(patientID))
# Get clonality table
itable <- LymphoSeq2::clonality(stable) %>%
          dplyr::left_join(sample_table, by = "repertoire_id") %>%
          dplyr::left_join(patient_table, by = "patientID")
# Merge metadata with ImmunoSeq data
atable <- atable %>%
          dplyr::left_join(sample_table, by = "repertoire_id") %>%
          dplyr::left_join(patient_table, by = "patientID")
stable <- stable %>%
          dplyr::left_join(sample_table, by = "repertoire_id") %>%
          dplyr::left_join(patient_table, by = "patientID")
ntable <- ntable %>%
          dplyr::left_join(sample_table, by = "repertoire_id") %>%
          dplyr::left_join(patient_table, by = "patientID")


# Create RDS objects
save(stable, file = "stable_cfar.rda")
save(atable, file = "atable_cfar.rda")
save(ntable, file = "ntable_cfar.rda")
save(itable, file = "itable_cfar.rda")
save(sample_table, patient_table, vl_table, cc_table, hla_table, infection_table, file = "mframe_cfar.rda")