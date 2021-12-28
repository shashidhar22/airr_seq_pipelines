# Author : Shashidhar Ravishankar
# Date: 09/23/2021
# Title : Format Metadata
# Description : This script formats the metadata for the cfar and dean data set 

# Load required packages
library(tidyverse)
library(readxl)
library(lubridate)
library(forcats)
# Accept arguments from nextflow
option_list <- list(
    optparse::make_option(c("-c", "--cfar_path"), 
        type="character",  
        help="Path to CFAR metadata",
        dest="cpath"),
    optparse::make_option(c("-n", "--norm_path"), 
        type="character",  
        help="Path to Dean metadata",
        dest="npath")                
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list), 
                                print_help_and_exit = TRUE)

# Get the HLA information from Dean metadata
dean_metadata <- read_tsv(parser$npath)
# The metadata available in the Dean metadata contains the HLA information in a
# text field with all the other metadata. The code below will first extract the
# HLA string, seprate them into different rows, collect the first two alleles
# for each HLA-A and HLA-B and then pivot the table again to generate the final
# metadata table in the same format as that from CFAR.
dean_hla <- dean_metadata %>%
    filter(str_detect(sample_tags, "HLA")) %>%
    select(sample_name, sample_tags) %>%
    mutate(hla_info = str_extract_all(sample_tags, "HLA-\\w\\*\\d+")) %>%
    rowwise() %>%
    group_by(sample_name) %>%
    summarize(hla_type = toString(unique(hla_info[[1]]))) %>%
    filter(hla_type != "") %>%
    separate_rows(hla_type, sep = ",") %>%
    group_by() %>% 
    mutate(allele = str_remove(hla_type, "HLA-"),
        allele_name = str_extract(allele, "\\w+")) %>%
    ungroup() %>%
    group_by(sample_name, allele_name) %>%
    mutate(allele_number = row_number()) %>%
    filter(allele_number < 3) %>%
    mutate(allele_name = str_c(allele_name, allele_number, sep = "")) %>%
    pivot_wider(id_cols = sample_name, values_from = allele, 
        names_from = allele_name) %>%
    mutate(B2 = if_else(is.na(B2), B1, B2), A2 = if_else(is.na(A2), A1, A2)) %>%
    rename(aAllele1 = A1, aAllele2 = A2, bAllele1 = B1, bAllele2 = B2, 
        patientID = sample_name) %>%
    mutate(aAllele1 = str_trim(aAllele1), aAllele2 = str_trim(aAllele2), 
        bAllele1 = str_trim(bAllele1), bAllele2 = str_trim(bAllele2),
        sampleID = patientID, timePoint = "Normal")
        
# Get the HLA information from CFAR metadata
cfar_metadata <- read_excel(parser$cpath, 
        sheet = "CFAR_Sample_info",
        col_types = c("text", "text", "date", "text", "text")) %>%
    filter(str_detect(sampleID, "CFAR")) 

cfar_hla <- read_excel(parser$cpath, sheet = "CFAR_HLA_info",
        col_types = c("text", "date", "text", "text", "text", "text", "text", 
        "text", "text", "text", "text", "text", "text", "text", "text", "text", 
        "text", "text", "text", "text")) %>%
    right_join(cfar_metadata, by = "patientID") 
# Merge the HLA information from Dean and CFAR metadata
study_hla_table <- bind_rows(cfar_hla, dean_hla) %>% 
    mutate(dateART = as_date(dateART), dateTCR = as_date(dateTCR))

# The metadata will be used at many different steps of the analysis for various
# purposes. The code below will generate multiple views of the metadata that can
# be used for different purposes.

# Let's first combine the Sample information table and patient information table
# for the Dean and CFAR data sets.
cfar_sample_table <- read_excel(parser$cpath, sheet = "CFAR_Sample_info",
        col_types = c("text", "text", "date", "text", "text")) %>%
    filter(str_detect(sampleID, "CFAR")) 
dean_sample_table <- read_tsv(parser$npath) %>% 
    select(sample_name) %>%
    rename(sampleID = sample_name) %>%
    mutate(patientID = sampleID,
        timePoint = "Normal",
        time_group = "Normal")
study_sample_table <- bind_rows(cfar_sample_table, dean_sample_table) %>%
    mutate(patientID = as_factor(patientID),
        dateTCR = as_date(dateTCR),
        time_point = as_factor(timePoint),
        repertoire_id = as_factor(sampleID)) %>%
    select(-sampleID, -timePoint)

# Let's merge Patient information table next
cfar_patient_table <- read_excel(parser$cpath, sheet = "CFAR_Patient_info",
        col_types = c("text", "text", "text", "text", "date", "text", "text", "text", "text", "date")) %>% 
        mutate(birthYear = ymd(birthYear, truncated = 2L),
            deathYear = ymd(deathYear, truncated = 2L),
            dateART = as_date(dateART),
            dateMinEvidenceHIV = as_date(dateMinEvidenceHIV))
dean_patient_table <- read_tsv(parser$npath) %>%
    select(sample_name) %>%
    rename(patientID = sample_name) %>%
    mutate(HIV_status = "Negative")
study_patient_table <- bind_rows(cfar_patient_table, dean_patient_table) %>%
    mutate(patientID = as_factor(patientID),
        dateART = as_date(dateART),
        sex = as_factor(sex),
        dateMinEvidenceHIV = ymd(dateMinEvidenceHIV),
        age_at_art = (birthYear %--% dateART) %/% years(1))

# Let's format the Viral load, CD4 counts and infection tables 
study_viral_load <- read_excel(parser$cpath, sheet = "CFAR_Viral_load") %>%
    mutate(patientID = as_factor(patientID),
        dateCollection = as_date(dateCollection)) %>%
    rename(date_collection = dateCollection) %>%
    left_join(study_patient_table, by = "patientID") %>%
    mutate(time_point = if_else(date_collection <= dateART, 
        "Pre-ART", "Post-ART"))

study_cd4_count <- read_excel(parser$cpath, sheet = "CFAR_CD4_count",
        col_types = c("text", "date", "text", "numeric", "text")) %>%
    mutate(patientID = as_factor(patientID),
        collectionDate = as_date(collectionDate)) %>%
    rename(date_collection = collectionDate) %>%
    left_join(study_patient_table, by = "patientID") %>%
    mutate(time_point = if_else(date_collection <= dateART, 
        "Pre-ART", "Post-ART"))

study_infections <- read_excel(parser$cpath, sheet = "CFAR_Coinfections",
        col_types = c("text", "date", "date", "date", "date", "date", "date",
            "date", "date", "date", "date", "date", "date", "date", "date")) %>%
    mutate(patientID = as_factor(patientID),
            across(EBV:Esophagus_cancer, as_date)) 

# Store all these output to an RDA object
out_file <- "Study_meta_frame.rda"
save(study_sample_table, study_patient_table, study_hla_table, study_viral_load,
    study_cd4_count, study_infections, file = out_file)
    


