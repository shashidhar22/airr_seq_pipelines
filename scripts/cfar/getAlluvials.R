#' Generate the alluvial plot for the top 50 sequences and color them by ART 
#' status
#' @description
#' Generate the alluvial plots for the top 50 sequences from each individual and
#' color sequences by whether they were found before or after ART initiation
#' @param amino_table Path to the amino acid table for the study
#' @return Alluvial plots for each individual in the study
library(tidyverse)
library(LymphoSeq2)
library(readxl)
library(optparse)
library(igraph)
library(lubridate)
library(patchwork)
library(tidygraph)
library(ggraph)
set.seed(12357) # Set the seed for reproducibility
option_list <- list(
    make_option(c("-a", "--amino_table"),
        type = "character",
        help = "Path to the sampled amino acid table",
        dest = "amino_table"),
    make_option(c("-m", "--meta_table"),
        type = "character",
        help = "Path to the study meta data in RDA format",
        dest = "meta_table"),
    make_option(c("-g", "--gliph_table"),
        type = "character",
        help = "Path to the GLIPH2 clusters table",
        dest = "gliph_table"),
    make_option(c("-h", "--gliph_hla_table"),
        type = "character",
        help = "Path to the GLIPH2 predicted HLA types",
        dest = "gliph_hla_table"),
    make_option(c("-c", "--mcpas_table"),
        type = "character",
        help = "Path to the McPAS-TCR database",
        dest = "mcpas_table"),
    make_option(c("-v", "--vdjdb_table"),
        type = "character",
        help = "Path to the VDJDb database",
        dest = "vdjdb_table")
)
parser <- optparse::parse_args(optparse::OptionParser(option_list=option_list, add_help_option=FALSE), 
    print_help_and_exit = TRUE) # Parse the arguments
# Load meta data
load(parser$meta_table)
# Read the amino acid table
load(parser$amino_table)
# Read GLIPH2 table and predicted HLA groups
study_gliph_table <- read_csv(parser$gliph_table) %>%
    select(pattern, TcRb, Sample) %>% 
    separate(Sample, c("repertoire_id", "status"), ":")
study_gliph_hla_table <- read_csv(parser$gliph_hla_table)
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
# Plot sankey diagrams 
plotSankey <- function(sample_table, vtable, cdtable, pbtable, tag_table, hla_table, cfar_publics) {
    set.seed(12345)
    ptid <- sample_table %>% 
        pull(patientID) %>% 
        unique() %>% 
        as.character()
    # Get hla vector 
    hla_vector <- hla_table %>% 
        filter(patientID == ptid) %>% 
        select(aAllele1:DRB345allele2) %>% 
        distinct() %>% 
        slice(1) %>% 
        as.character() %>% 
        str_extract("\\w+\\*\\d+")
    # Create pbfactor 
    public_pal <- c("EBV" = "#1f78b4", "CMV" = "#33a02c", "HIV" = "#e31a1c", "Other" = "#b15928")
    public_like_pal <- c("EBV_like" = "#a6cee3", "CMV_like" = "#b2df8a", "HIV_like" = "#fb9a99", "Other_like" = "#ffff99")

    pbtable <- pbtable %>% 
        filter(hla %in% hla_vector) %>% 
        mutate(palette = public_pal[pathogen])
    pbnames <- pbtable$junction_aa
    pbfactor <- pbtable$pathogen
    names(pbfactor) <- pbnames
    pbpal <- pbtable$palette
    names(pbpal) <- pbnames

    # Create public like CDR3 like (pblfactor)
    pbltable <- tag_table %>% 
        filter(patientID == ptid) %>% 
        mutate(pathogen = str_c(pathogen, "like", sep = "_"),
            palette = public_like_pal[pathogen]) 
    pblnames <- pbltable$junction_aa
    pblfactor <- pbltable$pathogen 
    names(pblfactor) <- pblnames
    pblpal <- pbltable$palette
    names(pblpal) <- pblnames

    # Filter samples
    sample_table <- sample_table %>%
        dplyr::filter(patientID == ptid)
    vtable <- vtable %>%
        dplyr::filter(patientID == ptid)
    cdtable <- cdtable %>%
        dplyr::filter(patientID == ptid)
    # Get top sequences
    sample_table <- sample_table %>%
        LymphoSeq2::cloneTrack() 
    # Collect repertoire_id information and response variable
    sample_name <- sample_table %>% 
        dplyr::pull(repertoire_id) %>% 
        unique()
    patient_name <- sample_table %>%
                    dplyr::pull(patientID) %>%
                    unique()
    #date_art <- sample_table %>%
    #            dplyr::pull(dateART) %>%
    #            unique()
    # Get viral load 
    vtable <- vtable %>% 
        dplyr::mutate(result = dplyr::case_when(viral_load == 0 ~ 1,
                TRUE ~ as.numeric(viral_load)),
            value_type = "Viral Load",
            days_rel = years_rel_to_art)
    # Get CD4 count
    cdtable <- cdtable %>%
        mutate(value_type = "CD4 count",
            days_rel = years_rel_to_art)
    # Get infection table
    cval <- c("ART", "EBV", "HBV", "HCV", "HSV1", "HSV2", "Chlamydia",
        "Gonnorrhea", "Syphilis", "Kaposi's sarcoma", "Anal cancer", "Colorectal cancer",
        "Skin cancer: squamous cell carcinoma", "Esophagus cancer", "Cervical cancer")
    cpal <- c("#e31a1c", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", 
        "#33a02c", "#ff7f00", "#6a3d9a", "#b15928", "#140063", "#7a6fa6", 
        "#20154d", "#037b80", "#33edf5")
    cdict <- c("ART", "EBV", "HBV", "HCV", "HSV1", "HSV2", "Chlamydia", 
        "Gonnorrhea", "Syphillis", "Kaposi_sarcoma", "Anal_cancer", "Colorectal_cancer",
        "Skin_cancer_squamous_cell_carcinoma", "Esophagus_cancer", "Cervical_cancer")
    names(cpal) <- cval
    names(cval) <- cdict
    # Gather all dates and create a date range for the samples
    date_lim <- c(-12, 12)
    # Color sequences by occurrence
    fcommon <- sample_table %>% 
        dplyr::filter(seen > 1 & duplicate_count > 1) %>%
        dplyr::group_by(junction_aa) %>%
        dplyr::summarize(dateEarliest = min(years_rel_to_art)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(foundStatus = dplyr::if_else(dateEarliest <= 0, "Pre-ART", "Post-ART"),
            foundStatus = dplyr::if_else(junction_aa %in% names(pbfactor), pbfactor[junction_aa], foundStatus),
            foundStatus = dplyr::if_else(junction_aa %in% names(pblfactor), pblfactor[junction_aa], foundStatus),
            foundStatus = dplyr::if_else(junction_aa %in% cfar_publics$junction_aa, "CFAR", foundStatus))
    # Now that we have identified the status for each sequence and the earliest date of occurrence. Let pull up the to 50 sequences
    atable <- sample_table
    sample_table <- sample_table %>%
        topSeqs(top = 50)
    fcommon <- fcommon %>%
        filter(junction_aa %in% sample_table$junction_aa)
    all_seq <- unique(sample_table$junction_aa)
    pre_seq <- fcommon %>% filter(foundStatus == "Pre-ART") %>% dplyr::pull(junction_aa)
    pre_pal <- if (length(pre_seq) > 0) p1(length(pre_seq)) else NULL
    pos_seq <- fcommon %>% filter(foundStatus == "Post-ART") %>% dplyr::pull(junction_aa)
    pos_pal <- if (length(pos_seq) > 0) p2(length(pos_seq)) else NULL
    pub_seq <- fcommon %>% filter(foundStatus %in% pbfactor) %>% dplyr::pull(junction_aa)
    pub_pal <- if (length(pub_seq) > 0) as.vector(pbpal[pub_seq]) else NULL
    pub_like_seq <- fcommon %>% filter(foundStatus %in% pblfactor) %>% dplyr::pull(junction_aa)
    pub_like_pal <- if (length(pub_like_seq) > 0) as.vector(pblpal[pub_like_seq]) else NULL
    cfar_seq <- fcommon %>% filter(foundStatus == "CFAR") %>% dplyr::pull(junction_aa)
    cfar_pal <- rep("#8e0152", length(cfar_seq))
    fnd_seq <- c(pre_seq, pos_seq, pub_seq, pub_like_seq, cfar_seq)
    mis_seq <- setdiff(all_seq, fnd_seq)
    mis_pal <- p3(length(mis_seq))
    tot_seq <- c(pre_seq, pos_seq, pub_seq, pub_like_seq, cfar_seq, mis_seq)
    tot_pal <- c(pre_pal, pos_pal, pub_pal, pub_like_pal, cfar_pal, mis_pal)
    names(tot_pal) <- tot_seq    
    pub_lab <- c(pbfactor[pub_seq], pblfactor[pub_like_seq])
    pub_sequences <- c(pub_seq, pub_like_seq)
    # Generate the alluvial plot to trace the frequency variation of top n Amino acid sequences
    sample_limits <- sample_table %>% 
        dplyr::mutate(yearsRel = forcats::as_factor(yearsRel)) %>% 
        dplyr::pull(yearsRel) %>%
        unique()
    sample_names <- sample_table %>%
        dplyr::pull(repertoire_id) %>%
        unique()
    names(sample_limits) <- sample_names
    sample_table <- sample_table %>% 
        dplyr::mutate(repertoire_id = yearsRel)
    sankey <- plotTrack(sample_table, alist = tot_seq, apal = tot_pal) +
        ggalluvial::geom_stratum(aes(y=duplicate_frequency)) + 
        ggplot2::geom_vline(xintercept = 0, color="red", size=1) + 
        ggplot2::theme_classic(base_size = 16) + 
        ggplot2::scale_x_continuous(limits = date_lim,
            breaks = seq(date_lim[1], date_lim[2], 2)) +
        ggplot2::xlab("Years from ART initiation") +
        ggplot2::ylab("Frequency") +
        ggplot2::theme(legend.position="none") +
        ggtitle(ptid)
    # Generate sample plot layout using Patchwork
    layout <- "
    AAAAAAAAA
    BBBBBBBBB
    CCCCCCCCC
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    DDDDDDDDD
    "
    sample_sankey <-  sankey + 
        plot_annotation(title = patient_name)
    out_path <- paste(patient_name, "alluvial.pdf", sep = "_")
    ggsave(out_path, sample_sankey, width = 8, height = 12, device = "pdf", units = "in")
    return(sankey)
}
# Get public tables 
getPublicTable <- function(atable, pbtable, htable) {
    # Amino table 
    ptid <- atable %>% 
        pull(patientID) %>%
        unique()
    # Get hla vector 
    hla_vector <- htable %>% 
        filter(patientID == ptid) %>% 
        select(aAllele1:DRB345allele2) %>% 
        distinct() %>% 
        slice(1) %>% 
        as.character() %>% 
        str_extract("\\w+\\*\\d+")
    
    atable_concise <- atable %>% 
        select(repertoire_id, junction_aa, duplicate_count, duplicate_frequency, patientID) 
    # Create pbfactor 
    pbtable <- pbtable %>% 
        filter(hla %in% hla_vector) %>% 
        mutate(patientID = ptid) %>% 
        inner_join(atable_concise, by = c("patientID", "junction_aa"))
    hla_table <- htable %>% 
        select(patientID, aAllele1:DRB345allele2) %>% 
        distinct()
    pbtable <- left_join(pbtable, hla_table, by = "patientID") 
    return(pbtable)
    
}

## Transform CD4 counts to be relative to time point nearest to date of treament

cc_table <- study_cd4_count %>% 
    rename(cd4_count = result) %>% 
    mutate(resultsRel = cd4_count, 
        patientID = forcats::fct_rev(patientID),
        yearsRel = years_rel_to_art,
        time_point = dplyr::if_else(yearsRel <= 0, "Pre-ART", "Post-ART"))

## Calculate days relative ART treatment for ImmunoSeq and Viral load table
vl_table <- study_viral_load %>%
    rename(viral_load = result) %>% 
    mutate(patientID = forcats::fct_rev(patientID),
        yearsRel = years_rel_to_art,
        time_point = dplyr::if_else(yearsRel <= 0, "Pre-ART", "Post-ART"))
# Load public TRB data

mptable <- readr::read_csv(parser$mcpas_table) %>% 
    filter(Species == "Human") %>% 
    rename(junction_aa = CDR3.beta.aa, antigen = Pathology, hla = MHC) %>%
    mutate(hla = str_remove(hla, "HLA-")) %>% 
    dplyr::select(antigen, junction_aa, hla)

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

filter_species <- c("HomoSapiens")
vdjtable <- readr::read_delim(parser$vdjdb_table, , delim="\t" ) %>%
    filter(species %in% filter_species) %>% 
    select(cdr3, antigen.species, mhc.a) %>%
    rename(junction_aa = cdr3, antigen = antigen.species, hla = mhc.a) %>%
    mutate(hla = str_remove(hla, "HLA-"))

pbtable <- mptable %>%
    bind_rows(vdjtable) %>%
    mutate(antigen = recode_factor(antigen, !!!antigen_table)) %>%
    distinct_all() %>%
    mutate(antigen = as.character(antigen), 
        pathogen = case_when(antigen == "Human immunodeficiency virus (HIV)" ~ "HIV", 
        antigen == "Epstein Barr virus (EBV)" ~ "EBV", 
        antigen == "Cytomegalovirus (CMV)" ~ "CMV", 
        TRUE ~ "Other"))  %>% 
    mutate(hla = str_extract(hla, "\\w+\\*\\d+"))

pbnames <- pbtable$junction_aa
pbfactor <- pbtable$pathogen
names(pbfactor) <- pbnames


atable <- amino_table %>% 
    left_join(study_sample_table, by = "repertoire_id") %>% 
    left_join(study_patient_table, by = "patientID") %>% 
    mutate(yearsRel = years_rel_to_art)
p1 <- colorRampPalette(c("#998ec3", "#cab2d6"))
p2 <- colorRampPalette(c("#f1a340", "#fdbf6f"))
p3 <- colorRampPalette(c("#737373", "#525252"))

study_public_table <- atable %>% 
    filter(str_detect(repertoire_id, "CFAR")) %>% 
    topSeqs(top = 50) %>%
    group_by(patientID) %>% 
    group_split() %>% 
    map(~getPublicTable(.x, pbtable, study_hla_table)) %>% 
    bind_rows()


study_public_gliph_groups <- study_public_table %>% 
    left_join(study_gliph_table, 
        by = c("repertoire_id", "junction_aa" = "TcRb")) %>% 
    select(antigen, junction_aa, hla, pathogen, patientID, repertoire_id, 
        duplicate_count, pattern, status)

public_gliph_table <- study_public_gliph_groups %>% 
    select(-patientID, -repertoire_id, -duplicate_count) %>% 
    distinct()

public_gliph_pattern_table <- public_gliph_table %>% 
    select(-junction_aa) %>% 
    distinct()

verified_publics_table <- study_gliph_hla_table %>% 
    select(pattern:Pvalue)  %>% 
    mutate(pattern = str_remove(pattern, "[gl]")) %>% 
    right_join(public_gliph_table, by = "pattern") %>% 
    filter(`'HLA allele with lowest Fisher Score'` == hla)

verified_public_gliph_table <- verified_publics_table %>% 
    select(-junction_aa) %>% 
    distinct()

top_amino_table <- amino_table %>% 
    topSeqs(top = 50) 

top_amino_gliph_table <- top_amino_table %>% 
    left_join(study_gliph_table, 
        by = c("repertoire_id", "junction_aa" = "TcRb")) %>% 
    filter(pattern %in% verified_public_gliph_table$pattern & 
        !(junction_aa %in% verified_publics_table$junction_aa)) %>% 
    left_join(verified_public_gliph_table, by = "pattern") %>% 
    left_join(study_sample_table, by = "repertoire_id") %>% 
    select(repertoire_id, patientID, junction_aa, pattern, hla, pathogen, duplicate_count, duplicate_frequency) %>% 
    left_join(study_hla_table, by = c("patientID", "repertoire_id" = "sampleID")) %>% 
    select(-timePoint, -time_group) %>% 
    pivot_longer(cols = aAllele1:DRB345allele2, 
        names_to = "Allele", values_to = "Allele_call") %>% 
    mutate(Allele_call = str_extract(Allele_call, "\\w+\\*\\d+")) %>% 
    filter(hla == Allele_call)

cfar_publics <- atable %>% 
    filter(str_detect(repertoire_id, "CFAR")) %>% 
    select(patientID, junction_aa) %>% 
    distinct() %>% 
    group_by(junction_aa) %>% 
    summarize(seen = n()) %>% 
    filter(seen > 1) %>% 
    filter(!(junction_aa %in% pbtable$junction_aa) & !(junction_aa %in% top_amino_gliph_table$junction_aa))


cfar_hla_table <- study_hla_table %>% 
    select(patientID, aAllele1:DRB345allele2) %>% 
    distinct()

cfar_gliph_table <- study_gliph_table %>% 
    left_join(study_sample_table, by = c("repertoire_id")) %>% 
    filter(status == "HIV") %>% 
    select(patientID, pattern, TcRb) %>% 
    distinct()

cfar_gliph_hla_table <- study_gliph_hla_table %>% 
    filter(Pvalue <= 0.05) %>% 
    rename(hla = `'HLA allele with lowest Fisher Score'`) %>% 
    mutate(pattern = str_remove(pattern, "[gl]")) %>% 
    select(pattern, hla, Pvalue)


cfar_public_table <- atable %>% 
    filter(junction_aa %in% cfar_publics$junction_aa) %>% 
    select(patientID, junction_aa) %>% 
    distinct() %>% 
    left_join(cfar_hla_table, by = "patientID") %>% 
    pivot_longer(aAllele1:DRB345allele2, names_to = "allele_name", values_to = "allele_call") %>% 
    mutate(allele_call = str_extract(allele_call, "\\w+\\*\\d+")) %>% 
    left_join(cfar_gliph_table, by = c("patientID", "junction_aa" = "TcRb")) %>% 
    left_join(cfar_gliph_hla_table, by = "pattern") %>% 
    filter(allele_call == hla) %>% 
    group_by(junction_aa, allele_call) %>% 
    summarise(nsamples = length(unique(patientID))) %>% 
    filter(nsamples > 1) %>% 
    pivot_wider(id_cols = junction_aa, names_from = "allele_call", values_from = "nsamples")



cfar_hootie_list <- tibble(junction_aa = c("CASSLSRGLLNEQFF", "CASSIGARGYEQYF",
    "CAWSVLKGETQYF", "CASSFELAPYEQYF", "CASSQDPGDSGNTIYF", "CASSFWDTIYF", 
    "CASSYSGGQGSSTDTQYF", "CASSTGNREKLFF", "CASSKQGDTGELFF", "CASSLAGWDTGELFF",
    "CASSLVLEGNEQFF", "CASSITGGGNTIYF", "CASSFYRGGGEQYF", "CASSIVTGAYNEQFF",
    "CASSSDTSPLHF", "CAWSEGGEAFF", "CASSPAGGRYEKLFF", "CASSIVGGEYEQYF"))

cfar_public_table <- bind_rows(cfar_public_table, cfar_hootie_list)

sankey_plot <- atable %>% 
    filter(str_detect(repertoire_id, "CFAR")) %>% 
    group_by(patientID) %>% 
    group_split() %>% 
    map(~plotSankey(.x, vl_table, cc_table, pbtable, top_amino_gliph_table, study_hla_table, cfar_public_table)) %>%
    wrap_plots(nrow=6, ncol=5)

out_path <- "CFAR_alluvial.pdf"
ggsave(out_path, sankey_plot, width = 30, height = 36, device = "pdf", units = "in", limitsize=FALSE)

