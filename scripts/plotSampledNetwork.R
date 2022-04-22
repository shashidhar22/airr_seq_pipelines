library(LymphoSeq2, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(readxl, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(viridis, quietly = TRUE)
library(fossil, quietly = TRUE)
library(ggraph, quietly = TRUE)
library(tidygraph, quietly = TRUE)
library(optparse, quietly = TRUE)
# Set opt parse for the following variables
# sample_name
# meta_object
# network_file
option_list <- list(
            make_option(c("-s", "--sample_name"), 
                        type="character",  
                        help="Sample name for which the network diagram is being plotted",
                        dest="sample"),
            make_option(c("-m", "--meta_data"), 
                        type="character",  
                        help="RDA file name for meta data object",
                        dest="meta"),
            make_option(c("-n", "--network_file"), 
                        type="character",  
                        help="Network file name",
                        dest="network")           
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
# Load metadata and productive amino acid table
load(parser$meta)
# Read network data and set names
sample <-  parser$sample
patient <- sample_table %>% 
            dplyr::filter(repertoire_id == sample) %>%
            pull(patientID) %>%
            unique()
date <- sample_table %>% 
        dplyr::filter(repertoire_id == sample) %>%
        pull(dateTCR) %>%
        unique() %>%
        as.character() %>%
        str_replace_all("-","_")
out_path <- paste(patient, date, "network.pdf", sep = "_")
net_table <- readr::read_tsv(parser$network) %>%
             mutate(ao_len = nchar(amino_one), 
                    at_len = nchar(amino_one)) %>%
             filter((ao_len >= 10 & ao_len <= 20) & (at_len >= 10 & at_len <= 20)) %>%
             select(amino_one, amino_two, editDist)
ngraph <- tidygraph::as_tbl_graph(net_table) %>% 
          dplyr::mutate(relatedness = centrality_degree(mode = 'in'), 
                        community = group_components(), 
                        nodecolor = if_else(community <= 10L, community, 11L)) %>% 
          dplyr::group_by(community) %>% 
          dplyr::mutate(occupency = n()) %>%
          dplyr::ungroup() 
nplot <- ggraph(ngraph, layout = 'drl') + 
            geom_edge_fan(aes(alpha = stat(index)), show.legend = FALSE) + 
            geom_node_point(aes(size = relatedness, colour = nodecolor)) + 
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
date_time <- sample_table %>% 
             dplyr::filter(repertoire_id == sample) %>%
             pull(dateTCR) %>%
             unique()
sgraph <- ngraph %>%
          as_tibble() %>%
          mutate(repertoire_id = sample,
                 patientID = patient,
                 dateTCR = date_time)


write_tsv(sgraph, paste(sample, "study", "NT.tsv", sep = "_"))