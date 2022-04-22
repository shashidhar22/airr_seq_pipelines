
load(parser$amino_table)
load(parser$meta_data)
study_public_table <- parser$public_path
study_gliph_table <- parser$gliph_path
study_gliph_hla_table <- parser$gliph_hla_path


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

study_public_table <- atable %>% 
    filter(str_detect(repertoire_id, "CFAR")) %>% 
    topSeqs(top = 50) %>%
    group_by(patientID) %>% 
    group_split() %>% 
    map(~getPublicTable(.x, pbtable, study_hla_table)) %>% 
    bind_rows()



study_gliph_table <- study_gliph_table %>% 
    select(pattern, TcRb, Sample) %>% 
    separate(Sample, c("repertoire_id", "status"), ":")

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
    select(-dateART, -dateTCR, -timePoint, -time_group) %>% 
    pivot_longer(cols = aAllele1:DRB345allele2, 
        names_to = "Allele", values_to = "Allele_call") %>% 
    mutate(Allele_call = str_extract(Allele_call, "\\w+\\*\\d+")) %>% 
    filter(hla == Allele_call)

top_amino_table_unrestricted <- top_amino_table %>% 
    left_join(study_gliph_table, 
        by = c("repertoire_id", "junction_aa" = "TcRb")) %>% 
    filter(pattern %in% public_gliph_table$pattern & 
        !(junction_aa %in% public_gliph_table$junction_aa)) %>% 
    left_join(public_gliph_pattern_table, by = "pattern") %>% 
    left_join(study_sample_table, by = "repertoire_id") %>% 
    select(repertoire_id, patientID, junction_aa, pattern, hla, pathogen, duplicate_count, duplicate_frequency) %>% 
    left_join(study_hla_table, by = c("patientID", "repertoire_id" = "sampleID")) %>% 
    select(-dateART, -dateTCR, -timePoint, -time_group) %>% 
    pivot_longer(cols = aAllele1:DRB345allele2, 
        names_to = "Allele", values_to = "Allele_call") %>% 
    mutate(Allele_call = str_extract(Allele_call, "\\w+\\*\\d+")) %>% 
    filter(hla == Allele_call)
