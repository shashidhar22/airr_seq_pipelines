library(LymphoSeq2)
library(tidyverse)     
library(lubridate)
library(readxl)
library(patchwork)
library(viridis)
library(fossil)
library(optparse)
# Accept arguments
option_list <- list(
    make_option(c("-p", "--pub_path"), 
                type="character",  
                help="Path to public sequences",
                dest="ppath")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
vdjdb_path <- paste(parser$ppath, "vdjdb_full.txt", sep = "/")
mcpas_path <- paste(parser$ppath, "McPAS-TCR.csv", sep = "/")
cmv2_path <- paste(parser$ppath, "CMV.publicTRB.2.csv", sep = "/")
cmv_path <- paste(parser$ppath, "CMV.publicTRB.csv", sep = "/")
ebv_path <- paste(parser$ppath, "EBV.publicTRB.csv", sep = "/")
flu_path <- paste(parser$ppath, "flu.M1.58.66.A2.TRB.csv", sep = "/")
hiv_path <- paste(parser$ppath, "HIV.publicTRB.csv", sep = "/")
hsv_path <- paste(parser$ppath, "HSV.publicTRB.csv", sep = "/")
ks_path <- paste(parser$ppath, "KS.publicTRB.csv", sep = "/")
vdjdb_frame <- readr::read_tsv(vdjdb_path) %>%
               dplyr::filter(species == "HomoSapiens" & vdjdb.score > 0) %>%
               dplyr::select(cdr3.alpha, cdr3.beta, mhc.a, mhc.b, mhc.class, antigen.epitope, antigen.species, vdjdb.score) %>%
               tidyr::pivot_longer(cols = cdr3.alpha:cdr3.beta,
                                   names_to = "CDR3_sequence_type",
                                   values_to = "CDR3_sequence") %>%
               dplyr::filter(!is.na(CDR3_sequence)) %>%
               dplyr::rename(CDR3 = CDR3_sequence, 
                             seq_type = CDR3_sequence_type,
                             hla_allele1 = mhc.a,
                             hla_allele2 = mhc.b,
                             hla_class = mhc.class,
                             epitope = antigen.epitope,
                             antigen = antigen.species,
                             score = vdjdb.score) %>%
               dplyr::mutate(seq_type = recode(seq_type, 
                                           cdr3.beta = "TRB", 
                                           cdr3.alpha = "TRA"))
mcpas_frame <- readr::read_csv(mcpas_path) %>%
               dplyr::filter(Species == "Human") %>%
               dplyr::select(CDR3.alpha.aa, CDR3.beta.aa, MHC, Pathology, Epitope.peptide, T.Cell.Type) %>%
               tidyr::pivot_longer(cols = CDR3.alpha.aa:CDR3.beta.aa,
                                   names_to = "seq_type",
                                   values_to = "CDR3") %>%
               dplyr::filter(!is.na(CDR3)) %>%
               dplyr::rename(hla_allele1 = MHC,
                             antigen = Pathology,
                             cell_type = T.Cell.Type,
                             epitope = Epitope.peptide) %>%
               dplyr::mutate(seq_type = recode(seq_type,
                                            CDR3.alpha.aa = "TRA",
                                            CDR3.beta.aa = "TRB"))
cmv2_frame <- readr::read_csv(cmv2_path) %>%
              dplyr::select(CDR3) %>%
              dplyr::mutate(antigen = "CMV",
                            seq_type = "TRB")
cmv_frame <- readr::read_csv(cmv_path) %>%
             dplyr::select(aminoAcid.alpha, aminoAcid.beta, peptide, mhcRestrictingAllele) %>%
             dplyr::mutate(antigen = "CMV") %>%
             tidyr::pivot_longer(cols = aminoAcid.alpha:aminoAcid.beta,
                                 names_to = "seq_type",
                                 values_to = "CDR3") %>%
             dplyr::filter(!is.na(CDR3)) %>%
             dplyr::rename(epitope = peptide,
                           hla_allele1 = mhcRestrictingAllele) %>%
             dplyr::mutate(seq_type = recode(seq_type,
                                             aminoAcid.beta = "TRB",
                                             aminoAcid.alpha = "TRA"))
ebv_frame <- readr::read_csv(ebv_path) %>%
             dplyr::select(aminoAcid.alpha, aminoAcid.beta, peptide, mhcRestrictingAllele, altRestrictingAllele) %>%
             dplyr::mutate(antigen = "EBV") %>%
             tidyr::pivot_longer(cols = aminoAcid.alpha:aminoAcid.beta,
                                 names_to = "seq_type",
                                 values_to = "CDR3") %>%
             dplyr::filter(!is.na(CDR3)) %>%
             dplyr::rename(epitope = peptide,
                           hla_allele1 = mhcRestrictingAllele,
                           hla_allele2 = altRestrictingAllele) %>%
             dplyr::mutate(seq_type = recode(seq_type,
                                             aminoAcid.beta = "TRB",
                                             aminoAcid.alpha = "TRA"))
flu_frame <- readr::read_csv(flu_path) %>%
             dplyr::select(aminoAcid.alpha, aminoAcid.beta, peptide, mhcRestricitingAllele, altRestrictingAllele) %>%
             dplyr::mutate(antigen = "Influenza") %>%
             tidyr::pivot_longer(cols = aminoAcid.alpha:aminoAcid.beta,
                                 names_to = "seq_type",
                                 values_to = "CDR3") %>%
             dplyr::filter(!is.na(CDR3)) %>%
             dplyr::rename(epitope = peptide,
                           hla_allele1 = mhcRestricitingAllele,
                           hla_allele2 = altRestrictingAllele) %>%
             dplyr::mutate(seq_type = recode(seq_type,
                                             aminoAcid.beta = "TRB",
                                             aminoAcid.alpha = "TRA"))
hiv_frame <- readr::read_csv(hiv_path) %>%
             dplyr::select(aminoAcid.beta, peptide, mhcRestrictingAllele, altRestrictingAllele) %>%
             dplyr::mutate(antigen = "HIV") %>%
             tidyr::pivot_longer(cols = aminoAcid.beta,
                                 names_to = "seq_type",
                                 values_to = "CDR3") %>%
             dplyr::filter(!is.na(CDR3)) %>%
             dplyr::rename(epitope = peptide,
                           hla_allele1 = mhcRestrictingAllele,
                           hla_allele2 = altRestrictingAllele) %>%
             dplyr::mutate(seq_type = recode(seq_type,
                                             aminoAcid.beta = "TRB"))
hsv_frame <- readr::read_csv(hsv_path) %>%
             dplyr::select(aminoAcid.alpha, aminoAcid.beta, peptide, mhcRestrictingAllele, altRestrictingAllele) %>%
             dplyr::mutate(antigen = "HSV") %>%
             tidyr::pivot_longer(cols = aminoAcid.alpha:aminoAcid.beta,
                                 names_to = "seq_type",
                                 values_to = "CDR3") %>%
             dplyr::filter(!is.na(CDR3)) %>%
             dplyr::rename(epitope = peptide,
                           hla_allele1 = mhcRestrictingAllele,
                           hla_allele2 = altRestrictingAllele) %>%
             dplyr::mutate(seq_type = recode(seq_type,
                                             aminoAcid.beta = "TRB",
                                             aminoAcid.alpha = "TRA"))
ks_frame <- readr::read_csv(ks_path) %>%
             dplyr::select(aminoAcid, specificity, mhcRestrictingAllele, certainty) %>%
             dplyr::mutate(antigen = "KS") %>%
             dplyr::rename(CDR3 = aminoAcid,
                           hla_allele1 = mhcRestrictingAllele,
                           score = certainty) %>%
             dplyr::mutate(seq_type = "TRB")

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
                   YFV = "Yellow fever virus",
                   `COVID-19` = "SARS-CoV-2",
                   "HTLV-1 (Chronic" = "HTLV-1")
public_frame <- bind_rows(vdjdb_frame, mcpas_frame) %>%
                bind_rows(cmv2_frame) %>%
                bind_rows(cmv_frame) %>%
                bind_rows(ebv_frame) %>%
                bind_rows(flu_frame) %>%
                bind_rows(hiv_frame) %>%
                bind_rows(hsv_frame) %>%
                bind_rows(ks_frame) %>%
                dplyr::mutate(antigen = recode_factor(antigen, !!!antigen_table))

save(public_frame, file = "Public_sequences.rda")