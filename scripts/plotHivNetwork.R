library(LymphoSeq2)
library(tidyverse)     
library(lubridate)
library(readxl)
library(patchwork)
library(viridis)
library(fossil)
library(tidygraph)
library(ggraph)
library(optparse)

# Accept arguments
option_list <- list(
        make_option(c("-i", "--hiv_data"), 
                    type="character",  
                    help="RDA file name for HIV specific sequences",
                    dest="hdata"),
        make_option(c("-m", "--meta_data"), 
                    type="character",  
                    help="RDA file name for meta data object",
                    dest="meta")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
# Load metadata and productive amino acid table
load(parser$meta)
htable <- read_tsv(parser$hdata)
# Plotter function for graph data
plotGraph <- function(htable) {
    sample <- htable %>%
              pull(repertoire_id) %>%
              unique()
    patient <- htable %>%
                pull(patientID) %>%
                unique()
    date <- htable %>%
            pull(dateTCR) %>%
            unique() %>%
            as.character() %>%
            str_replace_all("-","_")
    date_time <- htable %>%
                 pull(dateTCR) %>%
                 unique()
    out_path <- paste(patient, date, "specificity", "network.pdf", sep = "_")
    htable <- htable %>%
              dplyr::mutate(cdr_length = nchar(junction_aa)) %>%
              dplyr::filter(cdr_length >= 13 & cdr_length <= 16)
    net_table <- tidyr::expand_grid(amino_one = htable$junction_aa, 
                                    amino_two = htable$junction_aa) %>% 
                dplyr::filter_all(all_vars(!is.na(.))) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(editDist = adist(amino_one, amino_two, fixed=TRUE, partial=FALSE)) %>%
                dplyr::filter(editDist == 1) 
    ngraph <- tidygraph::as_tbl_graph(net_table) %>% 
            dplyr::mutate(relatedness = centrality_degree(mode = 'in'), 
                            community = group_components(), 
                            nodecolor = if_else(community <= 10L, community, 11L)) %>% 
            dplyr::group_by(community) %>% 
            dplyr::mutate(occupency = n()) %>%
            dplyr::ungroup() 
    nplot <- ggraph(ngraph, layout = 'drl') + 
            geom_edge_fan(aes(alpha = stat(index)), show.legend = FALSE) + 
            geom_node_point(aes(size = relatedness, colour = community)) + 
            theme_graph(foreground = 'steelblue', fg_text_colour = 'white') + 
            theme(legend.position = "none") +
            scale_color_viridis_c()
    ggsave(nplot,
            filename = out_path, 
            height = 15,
            width = 15,
            device = "pdf",
            units = "in",
            limitsize = FALSE)
    sgraph <- ngraph %>%
              as_tibble() %>%
              mutate(repertoire_id = sample,
                     patientID = patient,
                     dateTCR = date_time)
    return(sgraph)
}

ngraph <- htable %>% 
          group_by(repertoire_id) %>%
          group_split() %>%
          purrr::map(plotGraph) %>%
          bind_rows()
sample_name <- basename(parser$hdata) %>%
               tools::file_path_sans_ext() %>%
               stringr::str_remove("_specific")
out_file <- paste(sample_name, "NT.tsv", sep="_")
write_tsv(ngraph, out_file)