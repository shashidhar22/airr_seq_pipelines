getRarefiedCounts <- function(stable) {
    #Get the smallest repertoire size from the cohort
    rep_size <- stable %>% 
                dplyr::group_by(repertoire_id) %>% 
                dplyr::summarise(rep_size = sum(duplicate_count)) 
    #Uncount all junction sequences and randomly sample them to smallest rep size, finally group by unique rearrangement and add count and frequency of rearrangement
    rstable <- stable %>% 
            tidyr::uncount(duplicate_count) %>% 
            dplyr::group_by(repertoire_id) %>% 
            dplyr::sample_n(size = min(rep_size$rep_size)) %>% 
            dplyr::group_by(repertoire_id, junction, junction_aa, v_call, j_call, d_call, v_family, d_family, j_family) %>% 
            dplyr::summarise(duplicate_count = n()) %>% 
            dplyr::ungroup() %>% 
            dplyr::group_by(repertoire_id) %>% 
            dplyr::mutate(duplicate_frequency = duplicate_count/ sum(duplicate_count)) %>% 
            dplyr::ungroup()
    #Get clonality table where
    ritable <- LymphoSeq2::clonality(rstable) %>% 
               dplyr::select(repertoire_id, total_sequences) %>% 
               dplyr::rename(rarefied_count = total_sequences)
    itable <- LymphoSeq2::clonality(stable) %>% 
              dplyr::left_join(ritable, by = "repertoire_id") %>% 
    return(itable)
}